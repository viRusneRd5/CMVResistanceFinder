process QUARTO_SAMPLE_REPORT {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::quarto==1.9.38 conda-forge::r-rmarkdown=2.29 r-pacman r-cowplot r-tidyverse bioconda::bioconductor-genomicalignments r-gt r-knitr r-rmarkdown conda-forge::r-htmltools bioconda::bioconductor-rsamtools genomedk::tinytex"

    input:
    tuple val(meta), path(bam), path(bai), path(resistance), path(formatted) //sorted bam file, resistance and formatted output files from MERGE_FORMAT
    path(mut_bed) //this is the bed file listing the locations of the mutations for the report plot
    val(protocol) //one of amplicon or wgs
    val(report_author) //name that you want to appear at the top of the report
    val(bibliography) //.bib .bibtex format
    val(database_name) //name of the database used
    val(cmv_resis_db_version) //date
    val(min_depth) //sequencing depth threshold for calling mutations
    val(cmv_merlin_name) //reference sequence name to list on report
    val(gene_positions) //this is the start/end positions of the AMR genes for plotting coverage.

    output:
    tuple val(meta), path("${meta.id}.html")

    script:
    def ul27_pos = gene_positions['ul27']
    def ul51_pos = gene_positions['ul51']
    def ul54_pos = gene_positions['ul54']
    def ul56_pos = gene_positions['ul56']
    def ul89_a_pos = gene_positions['ul89_a']
    def ul97_pos = gene_positions['ul97']

    """
    cp ${projectDir}/bin/sample_report.qmd .

    quarto render \
        sample_report.qmd \
        -P report_author:"${report_author}" \
        -P bibliography:${bibliography} \
        -P database_name:${database_name} \
        -P database_version:${cmv_resis_db_version} \
        -P sample:${meta.id} \
        -P resistance:${resistance} \
        -P formatted:${formatted} \
        -P bam:${bam} \
        -P mut_bed:${mut_bed} \
        -P protocol:${protocol} \
        -P threshold:${min_depth} \
        -P reference:${cmv_merlin_name} \
        -P ul27_pos:${ul27_pos} \
        -P ul51_pos:${ul51_pos} \
        -P ul54_pos:${ul54_pos} \
        -P ul56_pos:${ul56_pos} \
        -P ul89_a_pos:${ul89_a_pos} \
        -P ul97_pos:${ul97_pos} \
        --output ${meta.id}.html
    """
}

