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
    Parameters are validated against nextflow_schema.json (see --help).
================================================================================
*/

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { FILTER_TAXA } from './modules/filter_taxa'
include { RECENTRIFUGE } from './modules/recentrifuge'
include { RCF_TO_ABUNDANCE } from './modules/rcf_to_abundance'
include { DIVERSITY_ANALYSIS } from './modules/diversity'

workflow {

    // Types, required params, file existence and allowed values come from nextflow_schema.json
    validateParameters()
    log.info paramsSummaryLog(workflow)

    def cleaned_samples = parseList(params.samples)
    if (cleaned_samples.size() < 2) {
        error "ERROR: --samples needs at least 2 sample IDs (diversity analysis compares samples)."
    }

    def cleaned_controls = parseList(params.controls)
    if (!cleaned_controls) {
        log.warn "WARNING: No --controls provided. Running without negative control integration."
    }

    def all_ids = cleaned_controls + cleaned_samples
    def duplicate_ids = all_ids.countBy { id -> id }.findAll { _id, n -> n > 1 }.keySet()
    if (duplicate_ids) {
        error "ERROR: IDs listed more than once across --samples/--controls: ${duplicate_ids.join(', ')}"
    }

    // 'domain' replaced 'superkingdom' in NCBI taxdumps from 2025
    def valid_ranks = ['domain', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'any']
    def abundance_ranks = parseList(params.abundance_ranks).collect { rank -> rank.toLowerCase() }
    def invalid_ranks = abundance_ranks.findAll { rank -> !valid_ranks.contains(rank) }
    if (!abundance_ranks || invalid_ranks) {
        error "ERROR: Invalid --abundance_ranks '${params.abundance_ranks}'. Use a comma-separated list of: ${valid_ranks.join(', ')}"
    }

    // auto: Recentrifuge's control-filtered counts when controls are given, raw counts otherwise
    def counts_source = params.abundance_source == 'auto' ? (cleaned_controls ? 'ctrl' : 'raw') : params.abundance_source
    if (counts_source == 'ctrl' && !cleaned_controls) {
        error "ERROR: --abundance_source ctrl needs --controls (Recentrifuge only filters against controls)."
    }

    // Load Kraken2 classified reads from taxprofiler output, keeping only the requested IDs
    def kraken_dir = "${params.taxprofiler_results}/kraken2/${params.kraken2_db_name}"
    def kraken_files = files("${kraken_dir}/*.classifiedreads.txt")
    if (!kraken_files) {
        error "ERROR: No Kraken2 results found at ${kraken_dir}/"
    }
    kraken_results = channel.fromList(matchIdsToFiles(all_ids, kraken_files, kraken_dir))

    // Step 1: Filter contaminants (the only contaminant filter; Recentrifuge gets no exclusions)
    FILTER_TAXA(
        kraken_results,
        file(params.taxdump),
        file(params.exclude_taxa)
    )

    // Step 2: Recentrifuge analysis
    RECENTRIFUGE(
        FILTER_TAXA.out.filtered_reads.toList().map { pairs -> [pairs.collect { p -> p[0] }, pairs.collect { p -> p[1] }] },
        file(params.taxdump),
        cleaned_samples,
        cleaned_controls
    )

    // Step 3: Abundance tables at multiple ranks
    RCF_TO_ABUNDANCE(
        RECENTRIFUGE.out.rcf_tables,
        file(params.taxdump),
        channel.fromList(abundance_ranks),
        counts_source,
        cleaned_controls
    )

    // Step 4: Diversity analysis
    DIVERSITY_ANALYSIS(
        RCF_TO_ABUNDANCE.out.abundance_tables.collect(),
        params.groups_file ? file(params.groups_file) : []
    )

    // Registered here with local copies of `workflow` and `params`, which are no longer
    // resolvable when the handler runs (both with Nextflow 25.10 and 26.04)
    def run = workflow
    def outdir = params.outdir
    workflow.onComplete {
        log.info """
            =====================================================
            LowBioPipe completed at: ${run.complete}
            Status: ${run.success ? 'SUCCESS' : 'FAILED'}
            Duration: ${run.duration}
            =====================================================
            Results: ${outdir}
            """.stripIndent()
    }
}

// Split a comma-separated parameter into trimmed, non-empty items
def parseList(value) {
    return value ? value.toString().tokenize(',').collect { item -> item.trim() }.findAll { item -> item } : []
}

// Pair each ID with exactly one Kraken2 file. A file belongs to the longest ID that
// prefixes its name followed by '.' or '_', so 'S1' never claims the file of 'S1_rep2'.
def matchIdsToFiles(ids, kraken_files, kraken_dir) {
    def owner = { f ->
        ids.findAll { id -> f.name.size() > id.size() && f.name.startsWith(id) && f.name[id.size()] in ['.', '_'] }
           .max { id -> id.size() }
    }
    def files_by_id = kraken_files.groupBy(owner)
    return ids.collect { id ->
        def hits = files_by_id[id] ?: []
        if (hits.size() == 0) {
            error "ERROR: No Kraken2 file for ID '${id}' in ${kraken_dir}/ (expected a name starting with '${id}.' or '${id}_')"
        }
        if (hits.size() > 1) {
            error "ERROR: ID '${id}' matches ${hits.size()} Kraken2 files (${hits.collect { h -> h.name }.join(', ')}). " +
                  "List each ID explicitly, or run taxprofiler with --perform_runmerging to get one file per sample."
        }
        [id, hits[0]]
    }
}
