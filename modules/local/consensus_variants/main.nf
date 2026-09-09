process CONSENSUS_VARIANTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/79e3c1e22b660e6a4f3655b1aeced00469a42bdff308be6e44910f4de0210ea0/data' :
        'community.wave.seqera.io/library/biopython_matplotlib_pandas_python_pruned:46d87e2ad1f8a063' }"

    input:
    tuple val(meta), path(fasta) // one multi-record FASTA per sample (one record per resistance gene)
    path(genes_bed)              // resistance-genes BED (chrom, start, end, gene)
    path(reference)              // CMV reference FASTA

    output:
    tuple val(meta), path("*.sanger.vcf")          , emit: vcf
    tuple val(meta), path("*.assessed.tsv")         , emit: assessed
    tuple val(meta), path("*.gene_consensus.fasta") , emit: gene_fasta
    tuple val("${task.process}"), val('python'), eval("python --version | sed '1!d; s/^.*Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def frameshifts = params.include_frameshifts ? '--include-frameshifts' : ''
    """
    fasta_consensus_variants.py \\
        --sample ${prefix} \\
        --fasta ${fasta} \\
        --genes-bed ${genes_bed} \\
        --reference ${reference} \\
        --min-identity ${params.fasta_min_identity} \\
        ${frameshifts} \\
        --output-vcf ${prefix}.sanger.vcf \\
        --output-assessed ${prefix}.assessed.tsv \\
        --output-gene-fasta ${prefix}.gene_consensus.fasta
    """
}
