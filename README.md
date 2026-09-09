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

**viRusneRd5/cmvresistancefinder** is a bioinformatics pipeline that identifies CMV antiviral-resistance mutations and compares them to a database of known clinically relevant mutations. It can start from sequencing reads (`--input_type fastq`, short- or long-read, aligned and variant-called) or from an already-assembled per-gene consensus (`--input_type fasta`, e.g. Sanger sequencing, aligned to the reference and diffed). **This pipeline is for research use only and should not be used for informing clinical practice.**

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the major pipeline sections and the types of output it produces. You're giving an overview to someone new to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test,conda_singularity` before running the workflow on actual data.

> [!IMPORTANT]
> This pipeline currently needs **both Conda and Apptainer/Singularity** available on the machine that runs it, and is run with `-profile conda_singularity`. Most tools run from containers, but a few steps (the Quarto reports) do not yet have a usable container and fall back to Conda. We hope to make it run on either engine alone in a future update.



First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```
Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

The `bin/create_samplesheet.py` script can help you make the samplesheet from fastq files within a directory. You can run this script as follows:

```
python create_samplesheet.py /path/to/fastq_dir samplesheet.csv
```

To start from a per-gene consensus instead (e.g. Sanger sequencing), add `--input_type fasta` and use a two-column samplesheet (`sample,fasta`) where each `fasta` is one multi-record FASTA per sample, one record per resistance gene with the gene name at the start of the header. See [docs/usage.md](docs/usage.md#samplesheet-input-fasta--sanger-consensus).

*Hostile Reference*

If you want to perform host-removal, you will need to download the hostile reference and index using bowtie2 build, as per below:

```bash
wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz

bowtie2-build /path/to/hostile_reference/human-t2t-hla.fa.gz /path/to/hostile_reference/hostile_bowtie_index
```
If you do not wish to perform host-removal, add `--skip_host_removal true` in place of `--hostile_reference_path </path/to/hostile_reference>` below.

*Kraken2 database*

If you wish to run Kraken2 to perform taxonomic classification on your reads prior to alignment to the reference, you will need to download the Kraken2 viral database. You can you this by:

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz
tar -xzvf /path/to/kraken_database/k2_viral_20260226.tar.gz
```
If you do not wish to perform kraken2 classification, add `--skip_kraken2 true` in place of `--kraken2_db </path/to/kraken2_virus_db>` below.

You can also use the Kraken2 standard database (up to date datases [found here](https://benlangmead.github.io/aws-indexes/k2)), but this database is larger and may take longer to run. Note that the pluspf database will not identify viruses, so should not be used here.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile conda_singularity \
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

This pipeline is still a work in progress and any feedback or suggestions are welcome!

## Credits

viRusneRd5/cmvresistancefinder was originally written by Jocelyn Hume.

We thank the following people for their extensive assistance in the development of this pipeline:

Dr Emma Sweeney - the University of Queensland

Dr Kym Lowry - the University of Queensland
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
