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
include { FASTA_INPUT                                } from '../subworkflows/local/fasta_input/main'
include { PIPELINE_PREPARATION                       } from '../subworkflows/local/pipeline_preparation/main'
include { HOSTILE_CLEAN                              } from '../modules/nf-core/hostile/clean/main'
include { KRAKEN2_KRAKEN2                            } from '../modules/nf-core/kraken2/kraken2/main'
include { FASTQC as RAW_FASTQC                       } from '../modules/nf-core/fastqc/main'
include { FASTQC as FINAL_FASTQC                     } from '../modules/nf-core/fastqc/main'
include { SAMTOOLS_STATS                             } from '../modules/nf-core/samtools/stats/main'

//Trimming
include { TRIMMOMATIC                                } from '../modules/nf-core/trimmomatic/main'
include { FASTP                                      } from '../modules/nf-core/fastp/main'
include { FASTPLONG                                  } from '../modules/nf-core/fastplong/main'

//Align to reference
include { BOWTIE2_ALIGN                              } from '../modules/nf-core/bowtie2/align/main'
include { MINIMAP2_ALIGN                             } from '../modules/nf-core/minimap2/align/main'
include { BWAMEM3_MEM                                } from '../modules/nf-core/bwamem3/mem/main'
include { SAMTOOLS_SORT                              } from '../modules/nf-core/samtools/sort/main'

//Variant calling
include { IVAR_VARIANTS                              } from '../modules/local/ivar/variants/main'
include { FREEBAYES                                  } from '../modules/local/freebayes/main'
include { BCFTOOLS_SORT                              } from '../modules/nf-core/bcftools/sort/main'
include { LOFREQ_CALL                                } from '../modules/nf-core/lofreq/call/main'
include { LOFREQ_FILTER                              } from '../modules/local/lofreq/filter/main'
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
include { SNPEFF_ANN            as SNPEFF_SANGER     } from '../modules/local/snpeff/ann'
include { SNPSIFT_EXTRACTFIELDS as SNPSIFT_SANGER    } from '../modules/local/snpsift/extractfields/main'
include { MERGE_FORMAT                               } from '../modules/local/merge_format/main'
include { QUARTO_SAMPLE_REPORT                       } from '../modules/local/quarto/sample_report/main'
include { SAMTOOLS_COVERAGE                          } from '../modules/local/samtools/coverage/main'
include { GENE_CONSENSUS                             } from '../modules/local/gene_consensus/main'
include { GENE_CONSENSUS as GENE_CONSENSUS_GENOTYPE  } from '../modules/local/gene_consensus/main'
include { GENOTYPE_CLASSIFY as GENOTYPE_CLASSIFY_GB  } from '../modules/local/genotype/classify/main'
include { GENOTYPE_MERGE                              } from '../modules/local/genotype/merge/main'
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
    PIPELINE_PREPARATION (
        params.cmv_ref_name,
        params.cmv_ref_fasta,
        params.cmv_ref_gff
    )
    ch_resis_mut_bed = PIPELINE_PREPARATION.out.resis_bed
    ch_cmv_ref_fasta = PIPELINE_PREPARATION.out.fasta
    ch_cmv_ref_fai = PIPELINE_PREPARATION.out.fai
    ch_cmv_ref_fai = ch_cmv_ref_fai.map { id, fai ->
            tuple(
                [id: id, single_end: false],
                fai
            )
    }
    ch_fasta_fai = ch_cmv_ref_fasta.join(ch_cmv_ref_fai).first()
    ch_cmv_aligner_index = PIPELINE_PREPARATION.out.cmv_ref_index
    ch_cmv_ref_gff = PIPELINE_PREPARATION.out.gff


    ch_genes_bed = Channel.fromPath(params.cmv_genes_bed).first()
    ch_genotype_genes_bed = Channel.fromPath(params.genotype_genes_bed).first()
    ch_gb_genotype_references = Channel.fromPath(params.gb_genotype_references).first()

    ch_snpsift_all = Channel.empty()
    ch_fasta_assessed = Channel.empty()

    //
    // SUBWORKFLOW: FASTA (Sanger / assembled per-gene consensus) input
    // -> per-sample VCF -> shared snpEff / SnpSift, then MERGE_FORMAT below.
    //
    if (params.input_type == 'fasta') {
        FASTA_INPUT (
            ch_samplesheet,                        // [ meta, fasta ] from PIPELINE_INITIALISATION
            ch_genes_bed,
            ch_cmv_ref_fasta.map { meta, fasta -> fasta }
        )
        ch_fasta_assessed = FASTA_INPUT.out.assessed

        SNPEFF_SANGER (
            FASTA_INPUT.out.vcf,
            PIPELINE_PREPARATION.out.snpeff_db,
            PIPELINE_PREPARATION.out.snpeff_config,
            ch_cmv_ref_fasta
        )
        SNPSIFT_SANGER (
            SNPEFF_SANGER.out.vcf
        )
        ch_snpsift_all = ch_snpsift_all.mix(SNPSIFT_SANGER.out.txt)
    }

    //
    // Read input (short-read / long-read): reads -> QC/trimming -> host removal -> alignment -> variant callers.
    // (everything down to `ch_merged_vcf` below is reads-only; skipped for --input_type fasta)
    //
    if (params.input_type != 'fasta') {

    //
    // SUBWORKFLOW: Create raw data channel
    //
    INPUT_CHECK(sample_input)
    ch_raw_fastqs = ch_raw_fastqs.mix(INPUT_CHECK.out.reads) // meta[id, reads]


    //
    // MODULE: Run Fastqc on raw fastqs.
    //
    RAW_FASTQC (
        ch_raw_fastqs
    )
    ch_multiqc_files = ch_multiqc_files.mix(RAW_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))


    //
    // MODULE: Run fastp
    //
    if (params.trimmer == "fastp" && params.input_type == "short-read") {
        FASTP (
            ch_raw_fastqs.map { meta, reads -> tuple(meta, reads, []) },
            false,
            false,
            false
        )
        ch_trimmed_fastq = ch_trimmed_fastq.mix(FASTP.out.reads)
    }

    //
    // MODULE: Run trimmomatic
    //
    if (params.trimmer == "trimmomatic" && params.adapter_file != null && params.input_type == "short-read") {
        TRIMMOMATIC (
                ch_raw_fastqs
            )
        ch_trimmed_fastq = ch_trimmed_fastq.mix(TRIMMOMATIC.out.trimmed_reads)
    } else if (params.trimmer == "trimmomatic" && params.adapter_file == null && params.input_type == "short-read"){
        println "Comment: You must provide an adapter fasta file for trimming and specify --adapter_file <path/to/adapter/file"
    } else if (params.trimmer == "trimmomatic" && params.adapter_file != null && params.input_type == "long-read"){
        println "Comment: This pipeline is not set up to perform trimming using trimmomatic on long-read sequencing data. Please use fastplong for trimming."
    }

    if (params.input_type == "long-read") {
        FASTPLONG (
            ch_raw_fastqs,
            [], //adapter file not needed as can auto detect
            false, //discard trimmed passed and only save reports
            false //save trimmed failed
        )
        ch_trimmed_fastq = ch_trimmed_fastq.mix(FASTPLONG.out.reads)
    }

    //
    // MODULE: Host-read removal (Hostile) and/or read classification (Kraken2).
    // `ch_final_reads` is what goes on to FastQC and alignment: host-removed reads
    // if host removal ran, otherwise the trimmed reads.
    //
    ch_final_reads = ch_trimmed_fastq

    if (!params.skip_host_removal) {
        ch_hostile_reference = Channel.value( [params.hostile_index, params.hostile_reference_path] )
        HOSTILE_CLEAN (
            ch_trimmed_fastq,
            ch_hostile_reference
        )
        ch_final_reads = HOSTILE_CLEAN.out.fastq
    }

    if (!params.skip_kraken2) {
        KRAKEN2_KRAKEN2 (
            ch_final_reads,
            params.kraken2_db,
            params.kraken2_output,
            params.kraken2_output
        )
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.collect{it[1]}.ifEmpty([]))
    }

    //
    // MODULE: Run Fastqc for trimmed and host-removed reads.
    //
    FINAL_FASTQC (
        ch_final_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FINAL_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))


    //
    // MODULE: Alignment with Bowtie2
    //
    if (params.aligner == 'bowtie2' && params.input_type == "short-read") {
        BOWTIE2_ALIGN (
            ch_final_reads,
            ch_cmv_aligner_index,
            ch_cmv_ref_fasta,
            false, //save unaligned
            true //samtools sort
        )
        ch_bam           = BOWTIE2_ALIGN.out.bam
    }

    //
    // MODULE: Alignment with minimap2
    //
    if (params.aligner == 'minimap2') {
        MINIMAP2_ALIGN (
            ch_final_reads,
            ch_cmv_aligner_index,
            true,
            'bai',
            false,
            false
        )
        ch_bam           = MINIMAP2_ALIGN.out.bam
    }

    //
    // MODULE: Alignment with bwa-mem3
    //
    if (params.aligner == 'bwa-mem3' && params.input_type == "short-read") {
        BWAMEM3_MEM (
            ch_final_reads,
            ch_cmv_aligner_index,
            ch_cmv_ref_fasta,
            true
        )
        ch_bam           = BWAMEM3_MEM.out.aligned
    }


    //
    // MODULE: Sort bam file and index
    //
    SAMTOOLS_SORT(
        ch_bam,
        ch_fasta_fai,
        'bai'
    )
    ch_sorted_bam = SAMTOOLS_SORT.out.bam
    ch_sort_bam_index = SAMTOOLS_SORT.out.index

    //
    // MODULE: Get bam file stats
    //
    ch_bam_bai = ch_sorted_bam.join(ch_sort_bam_index)
    SAMTOOLS_STATS (
        ch_bam_bai,
        ch_fasta_fai
    )
    multiqc_input_ch = SAMTOOLS_STATS.out.stats.map { meta, stats_file -> stats_file }.collect()

    //
    // MODULE: Optionally export a per-gene consensus FASTA per sample
    //
    if (params.export_gene_consensus) {
        GENE_CONSENSUS (
            ch_bam_bai,
            ch_genes_bed,
            params.min_depth,
            params.consensus_ambig_freq
        )
    }

    //
    // MODULE: Optional UL55 (gB) genotyping - independent of --export_gene_consensus,
    // this builds its own consensus over the genotyping gene(s) and classifies each
    // against a curated reference panel by nearest-neighbour identity. gN (UL73)
    // follows the same pattern once its reference panel is ready.
    //
    ch_genotype_tsvs = Channel.empty()
    if (params.genotype_gb) {
        GENE_CONSENSUS_GENOTYPE (
            ch_bam_bai,
            ch_genotype_genes_bed,
            params.min_depth,
            params.consensus_ambig_freq
        )
        GENOTYPE_CLASSIFY_GB (
            GENE_CONSENSUS_GENOTYPE.out.fasta,
            'UL55',
            '-',
            ch_gb_genotype_references
        )
        ch_genotype_tsvs = ch_genotype_tsvs.mix(GENOTYPE_CLASSIFY_GB.out.tsv)
    }


    // TODO: get stats on samples that passed minimum mapped reads


    //
    // Work out which variant caller(s) to run.
    // --variant_caller accepts "all", or a comma-separated subset of ivar, freebayes, lofreq
    // (e.g. "ivar,freebayes"). At least one valid caller must be selected.
    //
    def valid_callers = ['ivar', 'freebayes', 'lofreq']
    def selected_callers = params.variant_caller == 'all'
        ? valid_callers
        : params.variant_caller.tokenize(',').collect { it.trim().toLowerCase() }

    def unknown_callers = selected_callers - valid_callers
    if (unknown_callers) {
        error("Unknown --variant_caller value(s): ${unknown_callers.join(', ')}. " +
              "Valid options are 'all', or a comma-separated list from: ${valid_callers.join(', ')}.")
    }
    if (!selected_callers) {
        error("--variant_caller must select at least one caller ('all', or a comma-separated list from: ${valid_callers.join(', ')}).")
    }

    //
    // MODULE: Call variants with IVar, convert tsv to vcf, index, annotate with snpeff and extract with snpsift
    //
    if (params.input_type != 'fasta' && 'ivar' in selected_callers) {
        IVAR_VARIANTS (
            ch_sorted_bam,
            ch_fasta_fai,
            ch_cmv_ref_gff,
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
            PIPELINE_PREPARATION.out.snpeff_db,
            PIPELINE_PREPARATION.out.snpeff_config,
            ch_cmv_ref_fasta
        )

        SNPSIFT_IVAR (
            SNPEFF_IVAR.out.vcf
        )
        ch_snpsift_ivar = SNPSIFT_IVAR.out.txt
        ch_snpsift_all = ch_snpsift_all.mix(ch_snpsift_ivar)
    }


    //
    // MODULE: Call variants with Freebayes, index, annotate with snpeff and extract with snpsift
    //
    if (params.input_type != 'fasta' && 'freebayes' in selected_callers) {
        // amplicon: restrict calling to the six resistance genes (faster; the
        // reports only use those anyway). wgs: call across the whole genome.
        ch_freebayes = ch_sorted_bam.join(ch_sort_bam_index)
        ch_freebayes = params.protocol == 'wgs'
            ? ch_freebayes.map { meta, bam, bai -> tuple(meta, bam, bai, []) }
            : ch_freebayes.combine(ch_genes_bed).map { meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }

        FREEBAYES (
            ch_freebayes,
            ch_fasta_fai
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
            PIPELINE_PREPARATION.out.snpeff_db,
            PIPELINE_PREPARATION.out.snpeff_config,
            ch_cmv_ref_fasta
        )

        SNPSIFT_FREEBAYES (
            SNPEFF_FREEBAYES.out.vcf
        )
        ch_snpsift_freebayes = SNPSIFT_FREEBAYES.out.txt
        ch_snpsift_all = ch_snpsift_all.mix(ch_snpsift_freebayes)
    }

    //
    // MODULE: Call variants with Lofreq, index, annotate with snpeff and extract with snpsift
    //
    if (params.input_type != 'fasta' && 'lofreq' in selected_callers) {
        // amplicon: restrict calling to the six resistance genes. wgs: whole genome.
        ch_lofreq = params.protocol == 'wgs'
            ? ch_sorted_bam.map { meta, bam -> tuple(meta, bam, []) }
            : ch_sorted_bam.combine(ch_genes_bed).map { meta, bam, bed -> tuple(meta, bam, bed) }
        LOFREQ_CALL (
            ch_lofreq,
            params.cmv_ref_fasta
        )

        // lofreq call has no minimum-frequency option, so apply --min_freq here
        // (keeps the published lofreq vcf consistent with ivar / freebayes)
        LOFREQ_FILTER (
            LOFREQ_CALL.out.vcf
        )
        ch_lofreq_vcf = LOFREQ_FILTER.out.vcf

        INDEX_LOFREQ (
            ch_lofreq_vcf
        )
        ch_lofreq_index = INDEX_LOFREQ.out.index

        SNPEFF_LOFREQ (
            ch_lofreq_vcf,
            PIPELINE_PREPARATION.out.snpeff_db,
            PIPELINE_PREPARATION.out.snpeff_config,
            ch_cmv_ref_fasta
        )

        SNPSIFT_LOFREQ (
            SNPEFF_LOFREQ.out.vcf
        )
        ch_snpsift_lofreq = SNPSIFT_LOFREQ.out.txt
        ch_snpsift_all = ch_snpsift_all.mix(ch_snpsift_lofreq)
    }

    } // end of FASTQ-only reads/callers section

    //
    // MODULE: Merge vcf files, format output files, extract matches in database
    //
    // Groups the snpsift output of whichever caller(s) were selected by sample,
    // so this runs the same way whether 1, 2, or all 3 callers were chosen
    // (and the same way for a single "sanger" pseudo-caller from FASTA input).
    //
    ch_merged_vcf = ch_snpsift_all.groupTuple()

    MERGE_FORMAT (
        ch_merged_vcf,
        params.cmv_resis_db
    )
    ch_resistance = MERGE_FORMAT.out.resistance
    ch_formatted = MERGE_FORMAT.out.formatted


    //
    // MODULE: Reports
    //
    ch_db_name = file(params.cmv_resis_db).baseName
    gene_positions_ch = ch_genes_bed
        .map { bed_file ->
            def gene_positions = [:]
            bed_file.eachLine { line ->
                def (chr, start, end, gene) = line.tokenize('\t')
                gene_positions[gene.toLowerCase()] = "${start}-${end}"
            }
            gene_positions
        }

    if (params.input_type != 'fasta') {

        //
        // MODULE: Output html report for each sample using quarto
        //
        if (!params.skip_sample_report){
            // sixth tuple slot (assessed.tsv) is FASTA-only -> empty here.
            // seventh slot (genotype tsv) is only populated when genotyping is
            // enabled - joining against ch_genotype_tsvs when it's empty would
            // drop every sample, so branch instead of joining unconditionally.
            // One row per sample while only gB is wired up; revisit if/when gN
            // needs combining into the same per-sample file.
            if (params.genotype_gb) {
                ch_quarto = ch_sorted_bam.join(ch_sort_bam_index).join(ch_resistance).join(ch_formatted).join(ch_genotype_tsvs)
                    .map { meta, bam, bai, resistance, formatted, genotype -> tuple(meta, bam, bai, resistance, formatted, [], genotype) }
            } else {
                ch_quarto = ch_sorted_bam.join(ch_sort_bam_index).join(ch_resistance).join(ch_formatted)
                    .map { meta, bam, bai, resistance, formatted -> tuple(meta, bam, bai, resistance, formatted, [], []) }
            }
            QUARTO_SAMPLE_REPORT (
                ch_quarto,
                ch_resis_mut_bed,
                params.protocol,
                params.input_type,
                params.report_author,
                params.bibliography,
                ch_db_name,
                params.cmv_resis_db_version,
                params.min_depth,
                params.cmv_ref_name,
                gene_positions_ch
            )
        }

        //
        // MODULE: Get coverage data for bam file to output in report
        //
        //TODO: Figure out wgs vs amplicon report output.
        ch_bam_bai = ch_sorted_bam.join(ch_sort_bam_index)
        SAMTOOLS_COVERAGE (
            ch_bam_bai,
            ch_genes_bed,
            params.min_depth
        )
        ch_coverages = SAMTOOLS_COVERAGE.out.coverage.map { meta, coverage -> coverage }.collect()

        //
        // MODULE: Summarise coverage/depth and mutation results for whole run and output run summary report using quarto.
        //
        ch_resistances = ch_resistance.map { meta, resistance -> resistance }.collect()
        ch_total_count = ch_raw_fastqs.count()

        // Run-level genotype summary (one row per sample per genotyped gene) -
        // empty list when no genotyping was enabled, so the report can skip
        // the whole Genotyping section.
        ch_run_genotypes = Channel.value([])
        if (params.genotype_gb) {
            GENOTYPE_MERGE ( ch_genotype_tsvs.map { meta, tsv -> tsv }.collect() )
            ch_run_genotypes = GENOTYPE_MERGE.out.tsv
        }

        QUARTO_RUN_REPORT (
            ch_coverages,
            ch_resistances,
            params.report_author,
            params.run_name,
            ch_total_count,
            params.input_type,
            ch_run_genotypes
        )
    }

    //
    // MODULE: Reports for FASTA (Sanger consensus) input - no coverage/depth,
    // adds a "regions assessed" summary from the consensus differ instead.
    //
    if (params.input_type == 'fasta') {

        if (!params.skip_sample_report){
            // seventh tuple slot (genotype tsv(s)) is short/long-read only -> empty here
            ch_quarto_fa = ch_resistance.join(ch_formatted).join(ch_fasta_assessed)
                .map { meta, resistance, formatted, assessed -> tuple(meta, [], [], resistance, formatted, assessed, []) }
            QUARTO_SAMPLE_REPORT (
                ch_quarto_fa,
                ch_resis_mut_bed,
                params.protocol,
                params.input_type,
                params.report_author,
                params.bibliography,
                ch_db_name,
                params.cmv_resis_db_version,
                params.min_depth,
                params.cmv_ref_name,
                gene_positions_ch
            )
        }

        ch_resistances = ch_resistance.map { meta, resistance -> resistance }.collect()
        ch_assessed_all = ch_fasta_assessed.map { meta, assessed -> assessed }.collect()
        ch_total_count = ch_samplesheet.count()
        QUARTO_RUN_REPORT (
            ch_assessed_all,
            ch_resistances,
            params.report_author,
            params.run_name,
            ch_total_count,
            params.input_type,
            []
        )
    }

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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'cmvresistancefinder_software_'  + 'mqc_'  + 'versions.yml',
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
