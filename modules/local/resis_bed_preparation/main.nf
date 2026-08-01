process RESIS_BED_PREP {
    tag "prepare_resis_bed"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79e3c1e22b660e6a4f3655b1aeced00469a42bdff308be6e44910f4de0210ea0/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pandas_python_pruned:46d87e2ad1f8a063' }"

    input:
    val(resis_db)
    val(gff)

    output:
    path("mutation_list.txt"), emit: mut_list
    path("resistance_sites.bed"), emit: resis_bed
    tuple val("${task.process}"), val('python'), eval("python --version | sed '1!d; s/^.*Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get_mutation_list.py \
        --input ${resis_db} \
        --output mutation_list.txt

    mutation_bed.py \
    --gff ${gff} \
    --mutations mutation_list.txt \
    --output resistance_sites.bed
    """
}