/*
================================================================================
    RCF_TO_ABUNDANCE MODULE
================================================================================
    Converts Recentrifuge TSV output to abundance matrices at specified
    taxonomic ranks. Supports filtering and aggregation modes.
================================================================================
*/

process RCF_TO_ABUNDANCE {
    tag "${rcf_data.simpleName}_${rank}"
    label 'process_low'
    publishDir "${params.outdir}/abundance", mode: 'copy'

    // Same environment as RECENTRIFUGE: it ships the pandas needed to read its own output
    conda "bioconda::recentrifuge=2.1.2"
    container 'quay.io/biocontainers/recentrifuge:2.1.2--pyhdfd78af_0'

    input:
    tuple path(rcf_data), path(rcf_stats)
    path taxdump
    each rank
    val counts_source
    val controls

    output:
    path "abundance_*_counts_${rank}.tsv", emit: abundance_tables
    path "abundance_*_relative_${rank}.tsv", optional: true, emit: relative_tables
    path "abundance_*_taxa_metadata_${rank}.tsv", emit: metadata_tables

    script:
    // rcf_to_abundance.py finds ${rcf_stats} (<base>.stat.tsv) next to the data file
    def out_prefix = "abundance_${rcf_data.simpleName}"
    def relative = params.abundance_relative ? '--relative' : ''
    def aggregate = params.abundance_aggregate ? '--aggregate' : ''
    def keep_unranked = params.abundance_keep_unranked ? '--keep-unranked' : ''
    def controls_arg = controls ? "--controls ${controls.join(' ')}" : ''
    """
    rcf_to_abundance.py \\
        ${rcf_data} \\
        --taxdump ${taxdump} \\
        --source ${counts_source} \\
        ${controls_arg} \\
        --rank ${rank} \\
        --min-count ${params.abundance_min_count} \\
        --min-samples ${params.abundance_min_samples} \\
        --out-prefix ${out_prefix} \\
        ${relative} \\
        ${aggregate} \\
        ${keep_unranked}
    """
}
