/*
================================================================================
    DIVERSITY_ANALYSIS MODULE
================================================================================
    Comprehensive diversity analysis:
    - Alpha diversity (observed OTUs, Shannon, Simpson, Chao1)
    - Beta diversity (Bray-Curtis, Jaccard)
    - PCoA ordination
    - Hierarchical clustering heatmaps
    - PERMANOVA statistical tests (if groups provided)

    Automatically detects and processes all taxonomic levels.
================================================================================
*/

process DIVERSITY_ANALYSIS {
    label 'process_medium'
    publishDir "${params.outdir}/diversity", mode: 'copy'

    conda "bioconda::scikit-bio=0.5.9"
    container 'quay.io/biocontainers/scikit-bio:0.5.9--py39h3d4b85c_0'

    input:
    path abundance_tables
    path groups_file

    output:
    path "*/alpha_metrics_*.tsv", emit: alpha_metrics
    path "*/alpha_*.png", emit: alpha_plots
    path "*/beta_*.tsv", emit: beta_matrices
    path "*/pcoa_*_coords_*.tsv", emit: pcoa_coords
    path "*/pcoa_*.png", emit: pcoa_plots
    path "*/heatmap_*.png", emit: heatmaps
    path "*/permanova_*.txt", optional: true, emit: permanova_results

    script:
    def groups_arg = groups_file.name != 'NO_FILE' ? "--groups ${groups_file}" : ''
    """
    diversity_analysis.py \\
        --indir . \\
        --outdir . \\
        ${groups_arg}
    """
}
