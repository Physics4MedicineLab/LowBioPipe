/*
================================================================================
    FILTER_TAXA MODULE
================================================================================
    Filters Kraken2 classified reads by removing specified contaminant taxa.
    Supports ancestor and descendant filtering for comprehensive removal.
================================================================================
*/

process FILTER_TAXA {
    tag "${classifiedreads.simpleName}"
    label 'process_medium'
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'

    input:
    path classifiedreads
    path taxdump
    path exclude_taxa

    output:
    path "*.filtered.txt", emit: filtered_reads

    script:
    def ancestors = params.filter_include_ancestors ? '--include-ancestors' : ''
    def descendants = params.filter_include_descendants ? '--include-descendants' : ''
    """
    filter_taxa.py \\
        --taxdump ${taxdump} \\
        --exclude ${exclude_taxa} \\
        --indir . \\
        --outdir . \\
        --pattern "${classifiedreads.name}" \\
        ${ancestors} \\
        ${descendants}
    """
}
