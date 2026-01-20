/*
================================================================================
    RECENTRIFUGE MODULE
================================================================================
    Runs Recentrifuge for refined taxonomic profiling with control-based
    contamination assessment.
================================================================================
*/

process RECENTRIFUGE {
    label 'process_high'
    publishDir "${params.outdir}/recentrifuge", mode: 'copy'

    input:
    path classifiedreads
    path taxdump
    val samples
    val controls

    output:
    path "*.rcf.html", emit: html_reports
    path "*.rcf.xlsx", emit: rcf_files

    script:
    def exclude_args = file(params.exclude_taxa).readLines()
        .findAll { it && !it.startsWith('#') }
        .collect { it.split()[0] }
        .collect { "-x ${it}" }
        .join(' ')

    def control_files = controls.collect { ctrl ->
        classifiedreads.find { it.name.contains(ctrl) }
    }.findAll().collect { "-k ${it}" }.join(' ')

    def sample_files = samples.collect { sample ->
        classifiedreads.find { it.name.contains(sample) }
    }.findAll().collect { "-k ${it}" }.join(' ')

    def num_controls = controls.size()
    """
    rcf -n ${taxdump} \\
        -s KRAKEN \\
        -y ${params.recentrifuge_min_score} \\
        ${exclude_args} \\
        ${control_files} \\
        ${sample_files} \\
        -c ${num_controls} \\
        -o lowbiopipe_recentrifuge.rcf.html
    """
}
