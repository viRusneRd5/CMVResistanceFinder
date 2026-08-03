//
// Prepare reference genome files, index reference file, make snpeff database, prepare bed file for report.
//
include { BOWTIE2_BUILD                 } from '../../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main'
include { BWAMEM3_INDEX                 } from '../../../modules/nf-core/bwamem3/index/main'
include { MINIMAP2_INDEX                } from '../../../modules/nf-core/minimap2/index/main'


// local modules taken from nf-core/viralrecon
include { SNPEFF_BUILD                  } from '../../../modules/local/snpeff/build'

// my local module
include { RESIS_BED_PREP                } from '../../../modules/local/resis_bed_preparation/main'


workflow PIPELINE_PREPARATION {

    take:
    ref_name
    fasta
    gff
   // primer_bed

    main:

    //ch_versions = Channel.empty()

    //
    // Create chromosome sizes file
    //
    //ch_chrom = Channel.of( [params.cmv_ref_name, params.cmv_ref_fasta, params.cmv_ref_fai] )
    ch_chrom = Channel.of( [ref_name, fasta, []] )
    SAMTOOLS_FAIDX (
        ch_chrom,
        true
    )
    ch_fai         = SAMTOOLS_FAIDX.out.fai


    //ch_ref_fasta = Channel.of( [params.cmv_ref_name, params.cmv_ref_fasta] )
    ch_ref_fasta = Channel.value( [ [id: ref_name, single_end:false], fasta] )
    if (params.aligner == 'bowtie2') {
        //
        // Build reference bowtie2 index
        //
        //ch_ref_fasta = Channel.of( [ref_name, fasta] )
        BOWTIE2_BUILD (
            ch_ref_fasta
        )
        ch_ref_index = BOWTIE2_BUILD.out.index
    }

    if (params.aligner == 'bwa-mem3') {
        //ch_ref_fasta = Channel.value( [ [id: params.cmv_ref_name, single_end:false], params.cmv_ref_fasta] )
        //ch_ref_fasta = Channel.value( [ [id: ref_name, single_end:false], fasta] )
        BWAMEM3_INDEX (
            ch_ref_fasta
        )
        ch_ref_index = BWAMEM3_INDEX.out.index
    }

    if (params.aligner == 'minimap2') {
        //ch_ref_fasta = Channel.value( [ [id: ref_name, single_end:false], fasta] )
        MINIMAP2_INDEX (
            ch_ref_fasta
        )
        ch_ref_index = MINIMAP2_INDEX.out.index
    }


    //
    // Make snpEff database
    //
    //ch_fasta = Channel.value(file(fasta))
    ch_gff = Channel.value(file(gff))
    SNPEFF_BUILD (
        ch_ref_fasta,
        ch_gff
    )
    ch_snpeff_db     = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config
    //ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)


    //
    // Prepare the mutation bed file for report plots. Only needed if updating the database or using a different reference
    //
    if (params.cmv_ref_name != 'NC_006273.2' || params.cmv_genes_bed == null) {
        RESIS_BED_PREP (
            params.cmv_resis_db,
            ch_gff
        )
        ch_resis_bed = RESIS_BED_PREP.out.resis_bed
    } else {
        ch_resis_bed = params.mut_bed
    }


    emit:
    fasta                = ch_ref_fasta            // path: genome.fasta
    gff                  = ch_gff                  // path: genome.gff
    fai                  = ch_fai                  // path: genome.fai
    snpeff_db            = ch_snpeff_db            // path: snpeff_db
    snpeff_config        = ch_snpeff_config        // path: snpeff.config
    cmv_ref_index        = ch_ref_index            // path: genome.bt2 (bowtie2) .mmi (minimap2) bwamem3 directory
    resis_bed            = ch_resis_bed            // path: resistance_sites.bed
    //versions             = ch_versions             // channel: [ versions.yml ]
}
