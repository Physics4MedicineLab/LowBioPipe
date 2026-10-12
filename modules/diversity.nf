/*
================================================================================
    DIVERSITY_ANALYSIS MODULE
================================================================================
    Comprehensive diversity analysis:
    - Alpha diversity (observed OTUs, Shannon, Simpson, Chao1)
    - Beta diversity (Bray-Curtis, Jaccard, Aitchison)
    - PCoA ordination
    - Hierarchical clustering heatmaps
    - PERMANOVA statistical tests (if groups provided)

    Automatically detects and processes all taxonomic levels.
================================================================================
*/

process DIVERSITY_ANALYSIS {
    label 'process_medium'
    publishDir "${params.outdir}/diversity", mode: 'copy'

    conda "conda-forge::python=3.12 conda-forge::scikit-bio=0.7.4 conda-forge::seaborn=0.13.2 conda-forge::pandas=3.0.6 conda-forge::procps-ng"
    // Seqera Containers image built from the conda spec above
    container 'community.wave.seqera.io/library/python_scikit-bio_seaborn_pandas_procps-ng:b9637d4d60955f5b'

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
    path "*/permdisp_*.txt", optional: true, emit: permdisp_results

    script:
    def groups_arg = groups_file ? "--groups ${groups_file}" : ''
    """
    # Single-threaded BLAS: threaded summation order breaks PERMANOVA ties with the
    # observed statistic differently each run, so p-values would not be reproducible
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

    diversity_analysis.py \\
        --indir . \\
        --outdir . \\
        --heatmap-top-n ${params.diversity_heatmap_top_n} \\
        ${groups_arg}
    """
}
