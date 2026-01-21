/*
================================================================================
    RCF_TO_ABUNDANCE MODULE
================================================================================
    Converts Recentrifuge Excel output to abundance matrices at specified
    taxonomic ranks. Supports filtering and aggregation modes.
================================================================================
*/

process RCF_TO_ABUNDANCE {
    tag "${rcf_file.simpleName}_${rank}"
    label 'process_low'
    publishDir "${params.outdir}/abundance", mode: 'copy'

    conda "conda-forge::python=3.9 conda-forge::pandas>=1.5 conda-forge::openpyxl>=3.0 bioconda::ete3>=3.1"
    container 'python:3.9-slim'

    input:
    path rcf_file
    each rank

    output:
    path "abundance_*_counts_${rank}.tsv", emit: abundance_tables
    path "abundance_*_relative_${rank}.tsv", optional: true, emit: relative_tables
    path "abundance_*_taxa_metadata_${rank}.tsv", emit: metadata_tables

    script:
    def out_prefix = "abundance_${rcf_file.simpleName.replace('.rcf', '')}"
    def relative = params.abundance_relative ? '--relative' : ''
    def aggregate = params.abundance_aggregate ? '--aggregate' : ''
    def keep_unranked = params.abundance_keep_unranked ? '--keep-unranked' : ''
    """
    rcf_to_abundance.py \\
        ${rcf_file} \\
        --rank ${rank} \\
        --min-count ${params.abundance_min_count} \\
        --min-samples ${params.abundance_min_samples} \\
        --out-prefix ${out_prefix} \\
        ${relative} \\
        ${aggregate} \\
        ${keep_unranked}
    """
}
