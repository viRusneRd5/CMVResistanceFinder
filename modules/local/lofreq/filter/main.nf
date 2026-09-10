process LOFREQ_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4' :
        'quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4' }"

    input:
    tuple val(meta), path(vcf) // vcf(.gz) from LOFREQ_CALL

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('lofreq'), eval("lofreq version | sed -n '1s/^.* //p'"), emit: versions_lofreq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // lofreq call has no minimum-frequency option; lofreq filter is its own tool
    // for applying an allele-frequency (and, if wanted, strand-bias / quality)
    // cutoff. ext.args carries e.g. "--no-defaults -a <min_freq>".
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq filter \\
        ${args} \\
        -i ${vcf} \\
        -o ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
