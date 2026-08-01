# viRusneRd5/cmvresistancefinder

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/viRusneRd5/cmvresistancefinder)
[![GitHub Actions CI Status](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/nf-test.yml/badge.svg)](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/linting.yml/badge.svg)](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/viRusneRd5/cmvresistancefinder)

## Introduction

**viRusneRd5/cmvresistancefinder** is a bioinformatics pipeline that identifies variants from short-read sequencing data and compares to a database of known clinically relevant mutations. This pipeline is for research use only and should not be used for informing clinical practice.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.



First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

*Hostile Reference*

If you want to perform host-removal, you will need to download the hostile reference and index using bowtie2 build, as per below:

```bash
wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz

bowtie2-build /path/to/hostile_reference/human-t2t-hla.fa.gz /path/to/hostile_bowtie_index
```
If you do not wish to perform host-removal, add ```bash --skip_host_removal false ``` in place of ```bash --hostile_reference_path <path/to/hostile_reference>``` below.

*Kraken2 database*

If you wish to run Kraken2 to perform taxonomic classification on your reads prior to alignment to the reference, you will need to download the Kraken2 viral database. You can you this by:

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz
tar -xzvf /path/to/k2_viral_20260226.tar.gz
```



Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile <singularity/conda> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --protocol amplicon \
   --hostile_reference_path <path/to/hostile_reference> \
   --kraken2_db <path/to/kraken2_virus_db \
   --report_author <your_name/institute> \
   --run_name <run_name>

```
> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).




## FAQ's

**Why are there two different bed files?**

The first bed file (```bash ch_genes_bed``` or ```bash assets/CMV_NC_006273.2_resis_genes.bed```) describes the start and end positions of the six genes involved in CMV antiviral resistance (UL27, UL51, UL54, UL56, UL89 and UL97). Note that UL89 has a splice site in the middle, so there are two coding positions for this gene. This bed file should look like this:

```tsv
NC_006273.2	33182	35008	UL27
NC_006273.2	74576	75049	UL51
NC_006273.2	78194	81922	UL54
NC_006273.2	84752	87304	UL56
NC_006273.2	138810	139697	UL89_A
NC_006273.2	133770	134906	UL89_B
NC_006273.2	141798	143921	UL97
```
If you are not using the merlin reference sequence NC_006273.2, then this will need to match the header name and start and end locations of the genes for your selected reference sequence.

The second bed file (```bash ch_resis_mut_bed``` or ```bash assets/resistance_sites.bed```) describes the start and end locations of known CMV antiviral resistant mutations. This is only used in the quarto report file to list and plot the locations with insufficient sequencing depth to call a mutation. This file must be recreated if you are not using the merlin reference sequence NC_006273.2 or if the phenotyping_database is updated.

**What is the difference between *amplicon* and *wgs* protocol?**

This primary purpose of this pipeline is to identify mutations in the six genes associated with antiviral resistance and provide characterisation data to contextualise mutations where available. However, I do know that people may want to perform whole genome sequencing and get variants across all CDS and look at coverage across the genome. For this reason, I made the pipeline with two protocols. The amplicon protocol is not exclusively for use on amplicon sequencing data, this protocol can be used on whole genome sequencing data. The only difference between these protocols is that the reports will only output variants and coverage data from the six resistance genes, ignoring variants elsewhere in the genome and making analysis slightly faster. The wgs protocol can also be used on amplicon sequencing data, and could prove useful to identify non-specific amplification, but this analysis will take a little longer.

**Why doesn't the pipeline remove duplicate reads using picard?**

Sequencing CMV usually requires enrichment, whether this is probe hybridisation with PCR or just PCR. This means we need to duplicate the reads to be able to detect it with sequencing, especially when dealing with host DNA as well. Yes, we are assuming that all reads are amplified equally to get accurate variant frequencies, but if we were to remove duplicates, we may inadvertantly enhance minority variants and overstate their importance.

**Why is the default threshold's minimum 40 read depth, read quality 20, 0.5 (5%) variant frequency etc.?**

I don't really have a good reason for you. We want good quality reads so we can be more confident in the variants we are calling. 40 seemed like a good depth threshold for decent confidence in a variant and less than 5% frequency is likely not clinically relevant, so this is what I set. If you need to adjust the depth or variant frequency thresholds, you can change these parameters in the command line when you run the pipeline, or make a custom parameter/config file. I haven't added in the parameter to adjust the read quality, but if people ask for it, I might be able to.


**What do I need to provide if I want to run this pipeline with my own reference sequence?**

You will need to provide the following:
1. fasta file
2. gff file 
      must have extension .gff, if it ends in .gff3 you will need to rename it
3. the name of your reference
      the CHROM name in the header of the fasta file - these must match! 
      e.g. "NC_006273.2" for merlin or "X17403.1" for AD169
4. bed file
      the locations of the six resistance genes as described above. Yes you need this even if using wgs protocol.

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile <singularity/conda> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --protocol amplicon \
   --hostile_reference_path <path/to/hostile_reference> \
   --kraken2_db <path/to/kraken2_virus_db \
   --report_author <your_name/institute> \
   --run_name <run_name> \
   --cmv_ref_fasta <path/to/reference/fasta> \
   --cmv_ref_gff <path/to/gff/file> \
   --cmv_ref_name "<ref_name>" \
   --cmv_genes_bed <path/to/bed>
```

**Can I run this pipeline with a fasta file input from Sanger sequencing?**

If the demand is high enough, I might be able to make this work. 

**What if I want to provide my own database?**

This is complicated. The database used in this pipeline was compiled from a huge systematic review. Honestly, it is a bit of a mess, but if you start changing the database it will break the report writing. If you are aware of some mutations that are missing from the database, please get in touch and I can add it into the next update.




## Credits

viRusneRd5/cmvresistancefinder was originally written by Jocelyn Hume.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use viRusneRd5/cmvresistancefinder for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
