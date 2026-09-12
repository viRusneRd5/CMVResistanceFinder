//
// Convert ivar tsv to a vcf file and bgzip for input into bcftools merge
//
include { IVAR_BGZIP        } from '../../../modules/local/ivar/bgzip/main.nf'
include { IVAR_CONVERT      } from '../../../modules/local/ivar/convert/main.nf'


workflow CONVERT_IVAR {
    take:
    ch_ivar_tsv                // channel: [ val(meta), [ ivar_tsv_out ] ]

    main:

    ivar_vcfgz         = Channel.empty()

    IVAR_CONVERT (
        ch_ivar_tsv
    )
    ch_ivar_vcf = IVAR_CONVERT.out.vcf


    IVAR_BGZIP (
        ch_ivar_vcf
    )
    ivar_vcfgz = IVAR_BGZIP.out.vcfgz

    emit:
    ivar_vcfgz  // channel: [ val(meta), [ vcf.gz ] ]
}
