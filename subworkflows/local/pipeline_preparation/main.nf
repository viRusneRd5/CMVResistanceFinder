//
// Uncompress and prepare reference genome files
//
include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main'
include { BWAMEM3_INDEX                 } from '../../../modules/nf-core/bwamem3/index/main'

// local modules taken from nf-core/viralrecon
include { SNPEFF_BUILD                  } from '../../../modules/local/snpeff/build'

workflow PREPARE {

    take:
    fasta
    gff
   // primer_bed

    main:

    ch_versions = Channel.empty()

    //
    // Create chromosome sizes file
    //
    ch_chrom = Channel.of( [params.cmv_merlin_name, params.cmv_merlin_fasta, params.cmv_merlin_fai] )
    SAMTOOLS_FAIDX (
        ch_chrom,
        true
    )
    ch_fai         = SAMTOOLS_FAIDX.out.fai.map { it[1] }


    ch_merlin_fasta = Channel.of( [params.cmv_merlin_name, params.cmv_merlin_fasta] )
    if (params.aligner == 'bowtie2') {
        //
        // Build reference bowtie2 index
        //
        BOWTIE2_BUILD (
            ch_merlin_fasta
        )
        ch_merlin_index = BOWTIE2_BUILD.out.index
    }

    if (params.aligner == 'bwa-mem3') {
        ch_merlin_fasta = Channel.value( [ [id: params.cmv_merlin_name, single_end:false], params.cmv_merlin_fasta] )
        BWAMEM3_INDEX (
            ch_merlin_fasta
        )
        ch_merlin_index = BWAMEM3_INDEX.out.index
    }

    //
    // Make snpEff database
    //
    ch_fasta = Channel.value(file(fasta))
    ch_gff = Channel.value(file(gff))
    SNPEFF_BUILD (
        ch_fasta,
        ch_gff
    )
    ch_snpeff_db     = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config
    ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)


    emit:
    fasta                = ch_fasta                // path: genome.fasta
    gff                  = ch_gff                  // path: genome.gff
    fai                  = ch_fai                  // path: genome.fai
    snpeff_db            = ch_snpeff_db            // path: snpeff_db
    snpeff_config        = ch_snpeff_config        // path: snpeff.config
    cmv_merlin_index     = ch_merlin_index         // path: genome.bt2

    versions             = ch_versions             // channel: [ versions.yml ]
}
