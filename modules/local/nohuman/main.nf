process NOHUMAN {
    tag "$meta.id"
    label 'process_high' // kraken2 loads the whole database (~6GB for HPRC.r2) into memory

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nohuman:0.5.1--hfa8f182_1' :
        'quay.io/biocontainers/nohuman:0.5.1--hfa8f182_1' }"

    input:
    tuple val(meta), path(reads, stageAs: "input_reads/*")
    path(nohuman_db)

    output:
    tuple val(meta), path('*.nohuman.fq.gz')       , emit: fastq
    tuple val(meta), path('*.kraken2.report.txt')  , emit: report
    tuple val(meta), path('*.nohuman.log')         , emit: log
    tuple val("${task.process}"), val('nohuman'), eval("nohuman --version | sed 's/^nohuman //'"), emit: versions_nohuman, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def sorted_reads = meta.single_end ? [reads].flatten() : reads.sort { read -> read.name }
    def out_cmd      = meta.single_end ? "--out1 ${prefix}.nohuman.fq.gz" : "--out1 ${prefix}_1.nohuman.fq.gz --out2 ${prefix}_2.nohuman.fq.gz"
    """
    nohuman \\
        ${args} \\
        --threads ${task.cpus} \\
        --db ${nohuman_db} \\
        ${out_cmd} \\
        --kraken-report ${prefix}.nohuman.kraken2.report.txt \\
        ${sorted_reads.join(' ')} \\
        2> >(tee ${prefix}.nohuman.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outs   = meta.single_end ? "${prefix}.nohuman.fq.gz" : "${prefix}_1.nohuman.fq.gz ${prefix}_2.nohuman.fq.gz"
    """
    for f in ${outs}; do echo "" | gzip > \$f; done
    touch ${prefix}.nohuman.kraken2.report.txt
    touch ${prefix}.nohuman.log
    """
}
