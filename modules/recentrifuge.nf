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

    conda "bioconda::recentrifuge=2.1.2"
    container 'quay.io/biocontainers/recentrifuge:2.1.2--pyhdfd78af_0'

    input:
    tuple val(ids), path(classifiedreads)
    path taxdump
    val samples
    val controls

    output:
    path "*.rcf.html", emit: html_reports
    tuple path("*.rcf.data.tsv"), path("*.rcf.stat.tsv"), emit: rcf_tables

    script:
    // Contaminants were already removed by FILTER_TAXA, so no -x exclusions here.
    // Link each file as <ID>.krk so Recentrifuge names samples exactly by their ID.
    // Controls go first: rcf treats the first -c N inputs as negative controls.
    // TSV output (-e TSV) instead of Excel, whose 16,384-column limit is hit at ~500 samples.
    def file_by_id = [ids, [classifiedreads].flatten()].transpose().collectEntries()
    def ordered_ids = controls + samples
    def links = ordered_ids.collect { id -> "ln -s ${file_by_id[id]} ${id}.krk" }.join('\n')
    def inputs = ordered_ids.collect { id -> "-k ${id}.krk" }.join(' ')

    log.info "Recentrifuge: Processing ${samples.size()} samples with ${controls.size()} controls"
    """
    ${links}

    rcf -n ${taxdump} \\
        -s KRAKEN \\
        -y ${params.recentrifuge_min_score} \\
        ${inputs} \\
        -c ${controls.size()} \\
        -e TSV \\
        -o lowbiopipe_recentrifuge.rcf.html
    """
}
