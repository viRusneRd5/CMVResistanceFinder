//
// Turn per-sample FASTA (Sanger / assembled per-gene consensus) input into
// per-sample VCFs + a table of the region of each gene that was assessed.
//

include { CONSENSUS_VARIANTS } from '../../../modules/local/consensus_variants/main'

workflow FASTA_INPUT {

    take:
    ch_sample_fasta // channel: [ val(meta), path(fasta) ]  one multi-record FASTA per sample
    genes_bed       // path: resistance-genes BED
    reference       // path: CMV reference FASTA

    main:

    CONSENSUS_VARIANTS (
        ch_sample_fasta,
        genes_bed,
        reference
    )

    emit:
    vcf        = CONSENSUS_VARIANTS.out.vcf        // channel: [ val(meta), path(vcf) ]
    assessed   = CONSENSUS_VARIANTS.out.assessed   // channel: [ val(meta), path(tsv) ]
    gene_fasta = CONSENSUS_VARIANTS.out.gene_fasta // channel: [ val(meta), path(fasta) ]
}
