process QUARTO_RUN_REPORT {
    tag "quarto"
    label "process_medium"

    conda "${moduleDir}/environment.yml"

    input:
    path(coverage) //coverage.tsv files
    path(resistance) //resistance.tsv files
    val(report_author) //name that you want to appear at the top of the report
    val(run_name) //run name to list at the top of the report
    val(n_samples) //number of samples in samplesheet to list on report

    output:
    path("${run_name}.html")
    tuple val("${task.process}"), val('quarto'), eval('quarto -v'), emit: versions_quarto, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${run_name}"
    """
    cp ${projectDir}/bin/run_report.qmd .

    quarto render \
        run_report.qmd \
        -P report_author:"${report_author}" \
        -P run:${run_name} \
        -P n_samples:${n_samples} \
        --output ${prefix}.html
    """

    //TODO figure out stub here, add parameters to yml file, 
}