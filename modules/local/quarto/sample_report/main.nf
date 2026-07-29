process QUARTO_REPORT {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::quarto==1.9.38 conda-forge::r-quarto conda-forge::r-rmarkdown=2.29 r-pacman r-cowplot r-tidyverse bioconda::bioconductor-genomicalignments r-gt r-knitr r-rmarkdown bioconda::bioconductor-rtracklayer conda-forge::r-gridpattern==1.3.1 conda-forge::r-ggpattern bioconda::bioconductor-genomeinfodb r-ggcoverage"

    input:
    tuple val(meta), path(bam), path(resistance), path(formatted) //sorted bam file, resistance and formatted output files from MERGE_FORMAT
    path(mut_bed) //this is the bed file listing the locations of the mutations for the report plot
    val(min_depth) //sequencing depth threshold for calling mutations
    val(cmv_merlin_name) //reference sequence name to list on report

    output:
    tuple val(meta), path("${meta.id}.html")

    script:
    """
    cp ${projectDir}/bin/sample_report.qmd .

    quarto render \
        sample_report.qmd \
        -P sample:${meta.id} \
        -P resistance:${resistance} \
        -P formatted:${formatted} \
        -P bam:${bam} \
        -P mut_bed:${mut_bed} \
        -P threshold:${min_depth} \
        -P reference:${cmv_merlin_name} \
        --output ${meta.id}.html
    """
}