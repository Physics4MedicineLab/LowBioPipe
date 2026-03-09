#!/usr/bin/env nextflow

/*
================================================================================
    LowBioPipe - Post-Taxprofiler Processing for Low Biomass Microbiome Data
================================================================================
    Processes nf-core/taxprofiler output for low biomass samples:
    - Contaminant filtering
    - Recentrifuge taxonomic refinement
    - Abundance table generation
    - Diversity analysis

    Prerequisites: Run nf-core/taxprofiler first.
================================================================================
*/

nextflow.enable.dsl = 2

include { FILTER_TAXA } from './modules/filter_taxa'
include { RECENTRIFUGE } from './modules/recentrifuge'
include { RCF_TO_ABUNDANCE } from './modules/rcf_to_abundance'
include { DIVERSITY_ANALYSIS } from './modules/diversity'

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!params.taxprofiler_results) {
        error """
        ERROR: --taxprofiler_results is required.

        LowBioPipe processes output from nf-core/taxprofiler.
        Run taxprofiler first, then provide its results directory.

        Example:
            nextflow run main.nf \\
                --taxprofiler_results /path/to/taxprofiler/results \\
                --kraken2_db_name your_db_name \\
                --samples sample1,sample2 \\
                --controls blank1,blank2 \\
                -profile docker
        """
    }

    if (!params.kraken2_db_name) {
        error "ERROR: --kraken2_db_name is required (name of Kraken2 database used in taxprofiler)."
    }

    // Validate file paths exist
    if (!file(params.taxdump).exists()) {
        error "ERROR: --taxdump path does not exist: ${params.taxdump}\nDownload from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    }

    if (!file(params.exclude_taxa).exists()) {
        error "ERROR: --exclude_taxa file does not exist: ${params.exclude_taxa}"
    }

    // Validate samples and controls are provided
    if (!params.samples || params.samples.size() == 0) {
        error "ERROR: --samples is required. Provide comma-separated sample IDs."
    }

    // Clean sample IDs: handle String vs List and filter empty entries
    def cleaned_samples = (params.samples instanceof List ? params.samples : [params.samples]).collect { it?.trim() }.findAll { it }
    if (cleaned_samples.size() < params.samples.size()) {
        log.warn "WARNING: Empty sample IDs detected and removed from --samples"
    }
    if (cleaned_samples.size() == 0) {
        error "ERROR: --samples contains no valid sample IDs after cleaning."
    }

    if (!params.controls || params.controls.size() == 0) {
        log.warn "WARNING: No --controls provided. Running without negative control integration."
    }

    // Clean control IDs
    def cleaned_controls = !params.controls ? [] : (params.controls instanceof List ? params.controls : [params.controls]).collect { it?.trim() }.findAll { it }

    // Validate groups file if provided
    if (params.groups_file && !file(params.groups_file).exists()) {
        error "ERROR: --groups_file does not exist: ${params.groups_file}"
    }

    // Validate abundance_ranks contains valid values
    def valid_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'any']
    params.abundance_ranks.each { rank ->
        if (!valid_ranks.contains(rank.toLowerCase())) {
            error "ERROR: Invalid rank '${rank}'. Must be one of: ${valid_ranks.join(', ')}"
        }
    }

    log.info """
        =====================================================
        LowBioPipe v${workflow.manifest.version}
        =====================================================
        Taxprofiler results : ${params.taxprofiler_results}
        Kraken2 database    : ${params.kraken2_db_name}
        Output directory    : ${params.outdir}
        Taxdump             : ${params.taxdump}
        Exclusion list      : ${params.exclude_taxa}
        =====================================================
        """.stripIndent()

    // Load Kraken2 classified reads from taxprofiler output
    kraken_results = Channel
        .fromPath("${params.taxprofiler_results}/kraken2/${params.kraken2_db_name}/*.classifiedreads.txt")
        .ifEmpty { error "No Kraken2 results found at ${params.taxprofiler_results}/kraken2/${params.kraken2_db_name}/" }

    // Step 1: Filter contaminants
    FILTER_TAXA(
        kraken_results,
        file(params.taxdump),
        file(params.exclude_taxa)
    )

    // Step 2: Recentrifuge analysis
    RECENTRIFUGE(
        FILTER_TAXA.out.filtered_reads.collect(),
        file(params.taxdump),
        cleaned_samples,
        cleaned_controls
    )

    // Step 3: Abundance tables at multiple ranks
    RCF_TO_ABUNDANCE(
        RECENTRIFUGE.out.rcf_files,
        Channel.of(params.abundance_ranks).flatten()
    )

    // Step 4: Diversity analysis
    groups_ch = params.groups_file ? Channel.fromPath(params.groups_file) : Channel.of(file('NO_FILE'))

    DIVERSITY_ANALYSIS(
        RCF_TO_ABUNDANCE.out.abundance_tables.collect(),
        groups_ch
    )
}

workflow.onComplete {
    log.info """
        =====================================================
        LowBioPipe completed at: ${workflow.complete}
        Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration: ${workflow.duration}
        =====================================================
        Results: ${params.outdir}
        """.stripIndent()
}

def helpMessage() {
    log.info """
    =========================================================================
    LowBioPipe v${workflow.manifest.version}
    =========================================================================

    Post-taxprofiler processing for low biomass microbiome analysis.

    Usage:
        nextflow run main.nf --taxprofiler_results <path> --kraken2_db_name <name> [options]

    Required:
        --taxprofiler_results   Path to nf-core/taxprofiler output directory
        --kraken2_db_name       Name of Kraken2 database used in taxprofiler

    Sample configuration:
        --samples               Comma-separated list of sample IDs
        --controls              Comma-separated list of negative control IDs
        --groups_file           TSV file mapping samples to groups (for PERMANOVA)

    Reference data:
        --taxdump               Path to NCBI taxdump directory (default: data/taxdump)
        --exclude_taxa          Path to contaminant TaxID file (default: config/contaminants_example.txt)

    Contaminant filtering:
        --filter_include_ancestors      Also filter ancestor taxa (default: false)
        --filter_include_descendants    Also filter descendant taxa (default: false)

    Recentrifuge:
        --recentrifuge_min_score        Minimum score threshold (default: 10)

    Abundance tables:
        --abundance_ranks           Ranks to generate (default: species,genus,phylum)
        --abundance_aggregate       Aggregate taxa to rank (default: true)
        --abundance_min_count       Minimum total count threshold (default: 1)
        --abundance_min_samples     Minimum samples threshold (default: 2)
        --abundance_relative        Generate relative abundance (default: true)
        --abundance_keep_unranked   Keep taxa without resolvable rank (default: false)

    Output:
        --outdir                Output directory (default: results)

    Other:
        --help                  Show this message
        -profile                docker, singularity, conda, or test
        -resume                 Resume previous run

    Example:
        nextflow run main.nf \\
            --taxprofiler_results taxprofiler_output \\
            --kraken2_db_name standard_db \\
            --samples S1,S2,S3,S4 \\
            --controls BLANK1,BLANK2 \\
            --groups_file metadata/groups.tsv \\
            -profile docker
    """
}
