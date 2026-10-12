/*
================================================================================
    FILTER_TAXA MODULE
================================================================================
    Filters Kraken2 classified reads by removing specified contaminant taxa.
    Supports ancestor and descendant filtering for comprehensive removal.
================================================================================
*/

process FILTER_TAXA {
    tag "${id}"
    label 'process_low'
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'

    conda "conda-forge::python=3.12"
    container 'python:3.12-slim'

    input:
    tuple val(id), path(classifiedreads)
    path taxdump
    path exclude_taxa

    output:
    tuple val(id), path("*.filtered.txt"), emit: filtered_reads

    script:
    def ancestors = params.filter_include_ancestors ? "--include-ancestors --ancestors-max-rank ${params.filter_ancestors_max_rank}" : ''
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
