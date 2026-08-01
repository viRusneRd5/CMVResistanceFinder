process FREEBAYES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freebayes:1.3.10--hbefcdb2_0'
        : 'quay.io/biocontainers/freebayes:1.3.10--hbefcdb2_0'}"

    input:
    tuple val(meta), path(input), path(input_index), path(target_bed)
    tuple val(meta2), path(fasta), path(fasta_fai)

    output:
    tuple val(meta), path("*_freebayes.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('freebayes'), eval('freebayes --version 2>&1 | sed "s/version:\s*v//g"'), emit: versions_freebayes, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def targets_file = target_bed  ? "--target ${target_bed}"       : ""
    """
    freebayes \\
        -f ${fasta} \\
        ${targets_file} \\
        --ploidy 1 \\
        -F 0.05 \\
        --genotype-qualities \\
        ${args} \\
        ${input} > ${prefix}_freebayes.vcf

    bgzip ${prefix}_freebayes.vcf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_freebayes.vcf.gz
    """
}
