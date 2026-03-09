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

    conda "bioconda::recentrifuge>=1.12.0"
    container 'quay.io/biocontainers/recentrifuge:1.12.0--pyhdfd78af_0'

    input:
    path classifiedreads
    path taxdump
    val samples
    val controls

    output:
    path "*.rcf.html", emit: html_reports
    path "*.rcf.xlsx", emit: rcf_files

    script:
    // Parse exclusion taxa
    def exclude_args = file(params.exclude_taxa).readLines()
        .findAll { it && !it.startsWith('#') }
        .collect { it.split()[0] }
        .collect { "-x ${it}" }
        .join(' ')

    // Find control files with validation
    def control_matches = controls.collect { ctrl ->
        def match = classifiedreads.find { it.name.startsWith(ctrl) && ctrl.size() < it.name.size() && (it.name[ctrl.size()] == '.' as char || it.name[ctrl.size()] == '_' as char) }
        if (!match) {
            log.warn "WARNING: No file found matching control '${ctrl}'"
        }
        return match
    }.findAll()

    if (control_matches.size() == 0 && controls.size() > 0) {
        log.warn "WARNING: No control files matched. Check control IDs against file names."
    }

    def control_files = control_matches.collect { "-k ${it}" }.join(' ')

    // Find sample files with validation
    def sample_matches = samples.collect { sample ->
        def match = classifiedreads.find { it.name.startsWith(sample) && sample.size() < it.name.size() && (it.name[sample.size()] == '.' as char || it.name[sample.size()] == '_' as char) }
        if (!match) {
            log.warn "WARNING: No file found matching sample '${sample}'"
        }
        return match
    }.findAll()

    if (sample_matches.size() == 0) {
        error "ERROR: No sample files matched. Check sample IDs against file names in ${classifiedreads}"
    }

    def sample_files = sample_matches.collect { "-k ${it}" }.join(' ')

    def num_controls = control_matches.size()

    log.info "Recentrifuge: Processing ${sample_matches.size()} samples with ${num_controls} controls"
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
