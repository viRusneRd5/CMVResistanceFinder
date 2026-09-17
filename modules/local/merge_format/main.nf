process MERGE_FORMAT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79e3c1e22b660e6a4f3655b1aeced00469a42bdff308be6e44910f4de0210ea0/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pandas_python_pruned:46d87e2ad1f8a063' }"

    input:
    tuple val(meta), path(vcfs) // channel combining all vcf files for each sample
    path(cmv_resis_db)

    output:
    tuple val(meta), path("*.merged.tsv"), emit: merged
    tuple val(meta), path("*.formatted.tsv"), emit: formatted
    tuple val(meta), path("*.resistance.tsv"), emit: resistance
    tuple val("${task.process}"), val('python'), eval("python --version | sed '1!d; s/^.*Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    combine_snpsift_output.py \\
        --sample $prefix \\
        --output ${prefix}.merged.tsv \\
        $args \\
        ${vcfs.join(' ')}

    format_variants.py \\
    --input ${prefix}.merged.tsv \
    --output ${prefix}.formatted.tsv

    db_matching.py \
    --variants ${prefix}.formatted.tsv \
    --database $cmv_resis_db \
    --output ${prefix}.resistance.tsv
    """
}
