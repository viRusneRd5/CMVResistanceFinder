process GENOTYPE_MERGE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Just needs awk/coreutils to concatenate TSVs - reuse the same small,
    // already-verified samtools image the rest of the pipeline uses rather
    // than pull in a new, unverified container reference for this.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    path(tsvs) // one *.genotype.tsv per sample per genotyped gene, each with a header row

    output:
    path("run_genotypes.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' ${tsvs.join(' ')} > run_genotypes.tsv
    """
}
