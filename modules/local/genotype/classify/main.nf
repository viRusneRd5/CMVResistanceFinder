process GENOTYPE_CLASSIFY {
    tag "${meta.id}:${gene}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79e3c1e22b660e6a4f3655b1aeced00469a42bdff308be6e44910f4de0210ea0/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pandas_python_pruned:46d87e2ad1f8a063' }"

    input:
    tuple val(meta), path(consensus) // per-gene consensus FASTA (GENE_CONSENSUS output; may contain other genes too)
    val(gene)                        // gene name - must match a record header in consensus
    val(strand)                      // '+' or '-', the gene's strand on the reference genome
    path(reference_panel)            // FASTA of per-genotype reference sequences (header = genotype label)

    output:
    tuple val(meta), path("*.genotype.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed '1!d; s/^.*Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${gene}"
    """
    genotype_classify.py \\
        --sample ${meta.id} \\
        --gene ${gene} \\
        --strand '${strand}' \\
        --consensus ${consensus} \\
        --reference_panel ${reference_panel} \\
        --output ${prefix}.genotype.tsv \\
        ${args}
    """
}
