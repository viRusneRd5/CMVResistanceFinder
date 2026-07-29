/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_cmvresistancefinder_pipeline'
include { INPUT_CHECK                                } from '../subworkflows/local/input_check/main'
include { PREPARE                                    } from '../subworkflows/local/pipeline_preparation/main'
include { HOSTILE_CLEAN                              } from '../modules/nf-core/hostile/clean/main'
include { KRAKEN2_KRAKEN2                            } from '../modules/nf-core/kraken2/kraken2/main'
include { FASTQC as RAW_FASTQC                       } from '../modules/nf-core/fastqc/main'
include { FASTQC as FINAL_FASTQC                     } from '../modules/nf-core/fastqc/main'

//Trimming
include { TRIMMOMATIC                                } from '../modules/nf-core/trimmomatic/main'
include { FASTP                                      } from '../modules/nf-core/fastp/main'

//Align to reference
include { BOWTIE2_ALIGN                              } from '../modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN                             } from '../modules/nf-core/minimap2/align/main'
include { BWAMEM3_MEM                                } from '../modules/nf-core/bwamem3/mem/main'
include { SAMTOOLS_SORT                              } from '../modules/nf-core/samtools/sort/main'

//Variant calling
include { IVAR_VARIANTS                              } from '../modules/nf-core/ivar/variants/main'
include { FREEBAYES                                  } from '../modules/local/freebayes/main'
include { BCFTOOLS_SORT                              } from '../modules/nf-core/bcftools/sort/main'
include { LOFREQ_CALL                                } from '../modules/nf-core/lofreq/call/main'
include { BCFTOOLS_INDEX as INDEX_FREEBAYES          } from '../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as INDEX_LOFREQ             } from '../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as IVAR_INDEX               } from '../modules/nf-core/bcftools/index/main'
include { CONVERT_IVAR                               } from '../subworkflows/local/convert_ivar/main'
include { SNPEFF_ANN     as SNPEFF_FREEBAYES         } from '../modules/local/snpeff/ann'
include { SNPEFF_ANN     as SNPEFF_LOFREQ            } from '../modules/local/snpeff/ann'
include { SNPEFF_ANN     as SNPEFF_IVAR              } from '../modules/local/snpeff/ann'
include { SNPSIFT_EXTRACTFIELDS as SNPSIFT_FREEBAYES } from '../modules/local/snpsift/extractfields/main'
include { SNPSIFT_EXTRACTFIELDS as SNPSIFT_LOFREQ    } from '../modules/local/snpsift/extractfields/main'
include { SNPSIFT_EXTRACTFIELDS as SNPSIFT_IVAR      } from '../modules/local/snpsift/extractfields/main'
include { MERGE_FORMAT                               } from '../modules/local/merge_format/main'
include { QUARTO_REPORT                              } from '../modules/local/quarto/sample_report/main'
include { SAMTOOLS_COVERAGE                          } from '../modules/local/samtools/coverage/main'
include { QUARTO_RUN_REPORT                          } from '../modules/local/quarto/run_report/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CMVRESISTANCEFINDER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = Channel.empty()
    def ch_multiqc_files = Channel.empty()
    ch_raw_fastqs = Channel.empty()
    ch_trimmed_fastq = Channel.empty()

    sample_input = file(params.input, checkIfExists:true)

    // SUBWORKFLOW: Prepare the CMV reference genome and primers
    PREPARE (
        params.cmv_merlin_fasta,
        params.cmv_merlin_gff
        //params.cmv_primers_bed
    )

    //
    // Create raw data channel
    //
    INPUT_CHECK(sample_input)
    ch_raw_fastqs = ch_raw_fastqs.mix(INPUT_CHECK.out.reads) // meta[id, reads]
    //ch_raw_fastqs.view()


    //
    // MODULE: Run Fastqc
    //
    RAW_FASTQC (
        ch_raw_fastqs
    )


    if (params.trimmer == "fastp") {
        //
        // MODULE: Run fastp
        //
        FASTP (
            ch_raw_fastqs.map { meta, reads -> tuple(meta, reads, '') },
            false,
            false,
            false
        )
        ch_trimmed_fastq = ch_trimmed_fastq.mix(TRIMMOMATIC.out.trimmed_reads)
    }

    if (params.trimmer == "trimmomatic") {
        //
        // MODULE: Run trimmomatic
        //
        TRIMMOMATIC (
                ch_raw_fastqs
            )
        ch_trimmed_fastq = ch_trimmed_fastq.mix(TRIMMOMATIC.out.trimmed_reads)
    }

    if (!params.skip_host_removal && !params.skip_kraken2) {
        //
        // MODULE: Run Hostile to remove host reads
        //
        ch_hostile_reference = Channel.value( [params.hostile_index, params.hostile_reference_path] )
        HOSTILE_CLEAN (
            ch_trimmed_fastq,
            ch_hostile_reference
        )
        ch_final_reads = HOSTILE_CLEAN.out.fastq

        //
        // MODULE: Kraken2 to identify the number of reads remaining that belong to CMV
        //
        KRAKEN2_KRAKEN2 (
            ch_final_reads,
            params.kraken2_db,
            params.kraken2_output,
            params.kraken2_output
        )
        ch_multiqc_files =  ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    }


    if (params.skip_host_removal && !params.skip_kraken2) {
        //
        // MODULE: Kraken2 to identify the number of reads remaining that belong to CMV
        //
        KRAKEN2_KRAKEN2 (
            ch_trimmed_fastq,
            params.kraken2_db,
            params.kraken2_output,
            params.kraken2_output
        )
        ch_multiqc_files =  ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))

        ch_final_reads = ch_trimmed_fastq
    } else {
        ch_final_reads = ch_trimmed_fastq
    }

    //
    // MODULE: Run Fastqc
    //
    FINAL_FASTQC (
        ch_final_reads
    )


    if (params.aligner == 'bowtie2') {
        //
        // MODULE: Alignment with Bowtie2
        //
        ch_ref = Channel.value( [ [id: params.cmv_merlin_name, single_end:false], params.cmv_merlin_fasta] )
        cmv_merlin_index = PREPARE.out.cmv_merlin_index
        cmv_merlin_index = cmv_merlin_index.map { ref_name, index_directory ->
        [ [id: ref_name, single_end:false], index_directory]
        }
        //cmv_merlin_index = cmv_merlin_index.first()
        ch_ref = ch_ref.first()
        
        BOWTIE2_ALIGN (
            ch_trimmed_fastq,
            cmv_merlin_index,
            ch_ref,
            false, //save unaligned
            true //samtools sort
        )
        ch_bam           = BOWTIE2_ALIGN.out.bam
        ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    }


    if (params.aligner == 'minimap2') {
        //
        // MODULE: Alignment with minimap2
        //
        ch_ref = Channel.value( [ [id: params.cmv_merlin_name, single_end:false], params.cmv_merlin_fasta] )
        cmv_merlin_index = PREPARE.out.cmv_merlin_index
        cmv_merlin_index = cmv_merlin_index.map { ref_name, index_directory ->
        [ [id: ref_name, single_end:false], index_directory]
        }
        //cmv_merlin_index = cmv_merlin_index.first()
        ch_ref = ch_ref.first()
        MINIMAP2_ALIGN (
            ch_trimmed_fastq,
            ch_ref,
            true,
            'bai',
            false,
            false
        )
        ch_bam           = MINIMAP2_ALIGN.out.bam
    }


    if (params.aligner == 'bwa-mem3') {
        //
        // MODULE: Alignment with bwa-mem3
        //
        ch_ref = Channel.value( [ [id: params.cmv_merlin_name, single_end:false], params.cmv_merlin_fasta] )
        cmv_merlin_index = PREPARE.out.cmv_merlin_index
        //cmv_merlin_index = cmv_merlin_index.first()
        BWAMEM3_MEM (
            ch_trimmed_fastq,
            cmv_merlin_index,
            ch_ref,
            true
        )
        ch_bam           = BWAMEM3_MEM.out.aligned
    }


    //
    // MODULE: Sort bam file and index
    //

    ch_chrom = Channel.value( [params.cmv_merlin_name, params.cmv_merlin_fasta, params.cmv_merlin_fai] )
    SAMTOOLS_SORT(
        ch_bam,
        ch_chrom,
        'bai'
    )
    ch_sorted_bam = SAMTOOLS_SORT.out.bam
    ch_sort_bam_index = SAMTOOLS_SORT.out.index


    // To-do: get stats on samples that passed minimum mapped reads


    //
    // MODULE: Call variants with IVar
    //
    ch_vcf                    = Channel.empty()
    ch_tbi                    = Channel.empty()
    ch_snpsift_txt            = Channel.empty()
    if (params.variant_caller == "all" || params.variant_caller == "ivar") {
        IVAR_VARIANTS (
            ch_sorted_bam,
            params.cmv_merlin_fasta,
            params.cmv_merlin_fai,
            params.cmv_merlin_gff,
            true
        )
        ch_ivar_tsv           = IVAR_VARIANTS.out.tsv
        ch_multiqc_files = ch_multiqc_files.mix(IVAR_VARIANTS.out.tsv.collect{it[1]}.ifEmpty([]))

        CONVERT_IVAR (
            ch_ivar_tsv
        )
        ch_ivar_vcf = CONVERT_IVAR.out.ivar_vcfgz

        IVAR_INDEX (
            ch_ivar_vcf
        )
        ch_ivar_index = IVAR_INDEX.out.index

        SNPEFF_IVAR (
            ch_ivar_vcf,
            PREPARE.out.snpeff_db,
            PREPARE.out.snpeff_config,
            params.cmv_merlin_fasta
        )

        SNPSIFT_IVAR (
            SNPEFF_IVAR.out.vcf
        )
        ch_snpsift_ivar = SNPSIFT_IVAR.out.txt
    }

    if (params.variant_caller == "all" || params.variant_caller == "freebayes") {
        ch_freebayes = ch_sorted_bam.join(ch_sort_bam_index)
        ch_freebayes = ch_freebayes.map { meta, bam, bai -> 
        tuple(meta, bam, bai, params.cmv_genes_freebayes) 
        }
        FREEBAYES (
            ch_freebayes,
            ch_chrom
        )
        ch_freebayes_vcf = FREEBAYES.out.vcf

        BCFTOOLS_SORT (
            ch_freebayes_vcf
        )
        ch_freebayes_vcf = BCFTOOLS_SORT.out.vcf

        INDEX_FREEBAYES (
            ch_freebayes_vcf
        )
        ch_freebayes_index = INDEX_FREEBAYES.out.index

        SNPEFF_FREEBAYES (
            ch_freebayes_vcf,
            PREPARE.out.snpeff_db,
            PREPARE.out.snpeff_config,
            params.cmv_merlin_fasta
        )

        SNPSIFT_FREEBAYES (
            SNPEFF_FREEBAYES.out.vcf
        )
        ch_snpsift_freebayes = SNPSIFT_FREEBAYES.out.txt
    }

    if (params.variant_caller == "all" || params.variant_caller == "lofreq") {
        ch_lofreq = ch_sorted_bam.map { meta, bam -> 
        tuple(meta, bam, params.cmv_genes_freebayes)
        }
        LOFREQ_CALL (
            ch_lofreq,
            params.cmv_merlin_fasta
        )
        ch_lofreq_vcf = LOFREQ_CALL.out.vcf

        INDEX_LOFREQ (
            ch_lofreq_vcf
        )
        ch_lofreq_index = INDEX_LOFREQ.out.index

        SNPEFF_LOFREQ (
            ch_lofreq_vcf,
            PREPARE.out.snpeff_db,
            PREPARE.out.snpeff_config,
            params.cmv_merlin_fasta
        )

        SNPSIFT_LOFREQ (
            SNPEFF_LOFREQ.out.vcf
        )
        ch_snpsift_lofreq = SNPSIFT_LOFREQ.out.txt
    }

    if (params.variant_caller == "all") {
        ch_merged_vcf = ch_snpsift_ivar.join(ch_snpsift_freebayes)
        ch_merged_vcf = ch_merged_vcf.join(ch_snpsift_lofreq).map { meta, ivar, freebayes, lofreq ->
            tuple(meta, [ivar, freebayes, lofreq])
        }

        MERGE_FORMAT (
            ch_merged_vcf,
            params.cmv_resis_db
        )
        ch_resistance = MERGE_FORMAT.out.resistance
        ch_formatted = MERGE_FORMAT.out.formatted
    }

    ch_quarto = ch_sorted_bam.join(ch_resistance)
    ch_quarto = ch_quarto.join(ch_formatted)//.map { meta, bam, resistance, formatted ->
    //    tuple(meta, [bam, resistance, formatted])
    //}
    QUARTO_REPORT (
        ch_quarto,
        params.mut_bed,
        params.min_depth,
        params.cmv_merlin_name
    )

    ch_bam_bai = ch_sorted_bam.join(ch_sort_bam_index)
    SAMTOOLS_COVERAGE (
        ch_bam_bai,
        params.cmv_genes_bed,
        params.min_depth
    )

    ch_coverages = SAMTOOLS_COVERAGE.out.coverage.map { meta, coverage -> coverage }.collect()
    ch_resistances = ch_resistance.map { meta, resistance -> resistance }.collect()
    ch_total_count = ch_raw_fastqs.count()
    QUARTO_RUN_REPORT (
        ch_coverages,
        ch_resistances,
        params.run_name,
        ch_total_count
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'cmvresistancefinder_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'cmvresistancefinder'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
