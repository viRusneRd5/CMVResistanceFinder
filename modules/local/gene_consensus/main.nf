process GENE_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(bam), path(bai) //sorted bam file and index
    path(genes_bed)                       //resistance-genes BED (chrom, start, end, gene)
    val(min_depth)                        //positions below this depth become N
    val(ambig_freq)                       //minor-allele fraction above which an IUPAC code is emitted

    output:
    tuple val(meta), path("*.consensus.fasta"), emit: fasta
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    : > ${prefix}.consensus.fasta

    while IFS=\$'\\t' read -r chrom start end gene
    do
        [ -z "\$gene" ] && continue
        samtools consensus \\
            --format fasta \\
            --ambig \\
            --min-depth ${min_depth} \\
            --het-fract ${ambig_freq} \\
            --show-del yes --show-ins no \\
            -a \\
            -r "\$chrom:\$start-\$end" \\
            "${bam}" \\
        | awk -v gene="\$gene" 'NR==1 {print ">"gene; next} {print}' \\
        >> ${prefix}.consensus.fasta
    done < ${genes_bed}
    """
}
