process IVAR_BGZIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/33/33a1f2c7f36ec58339e41cbea096d121f606918778a91cfbef944b40ba7ce48b/data'
        : 'community.wave.seqera.io/library/htslib_xz:49c8c84af5c4b3b9'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_ivar.vcf.gz"), emit: vcfgz
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1! d; s/bgzip (htslib) //'"), topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bgzip -@ ${task.cpus} ${vcf}
    """

    stub:
    """
    touch "${vcf}.gz"
    """
}
