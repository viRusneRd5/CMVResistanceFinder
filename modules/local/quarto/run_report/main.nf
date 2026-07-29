process QUARTO_RUN_REPORT {
    tag "quarto"
    label "process_medium"

    conda "conda-forge::quarto==1.9.38 conda-forge::r-quarto conda-forge::r-rmarkdown=2.29 r-tidyverse bioconda::bioconductor-genomicalignments r-gt r-knitr r-rmarkdown conda-forge::r-plotly conda-forge::r-dt"

    input:
    path(coverage) //coverage.tsv files
    path(resistance) //resistance.tsv files
    val(run_name) //run name to list at the top of the report
    val(n_samples) //number of samples in samplesheet to list on report

    output:
    path("${run_name}.html")

    script:
    def prefix = task.ext.prefix ?: "${run_name}"
    """
    cp ${projectDir}/bin/run_report.qmd .

    quarto render \
        run_report.qmd \
        -P run:${run_name} \
        -P n_samples:${n_samples} \
        --output ${prefix}.html
    """
}