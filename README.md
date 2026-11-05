# CMVResistanceFinder

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/viRusneRd5/cmvresistancefinder)
[![GitHub Actions CI Status](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/nf-test.yml/badge.svg)](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/linting.yml/badge.svg)](https://github.com/viRusneRd5/cmvresistancefinder/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**CMVResistanceFinder** is a bioinformatics pipeline that identifies CMV antiviral-resistance mutations and compares them to a database of known clinically relevant mutations. It can start from sequencing reads (`--input_type short-read` or `long-read`, aligned and variant-called) or from an already-assembled per-gene consensus (`--input_type fasta`, e.g. Sanger sequencing, aligned to the reference and diffed). **This pipeline is for research use only and should not be used for informing clinical practice.**

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->

Read input:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter/quality trimming ([`fastp`](https://github.com/OpenGene/fastp), [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) or, for long reads, [`fastplong`](https://github.com/OpenGene/fastplong)), plus optional amplicon primer trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/), `--primer_sequences`)
3. Optional host-read removal ([`nohuman`](https://github.com/mbhall88/nohuman) by default, or [`Hostile`](https://github.com/bede/hostile)) and taxonomic classification ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
4. Alignment to the CMV reference ([`minimap2`](https://github.com/lh3/minimap2), [`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/) or [`BWA-MEM3`](https://github.com/bwa-mem2/bwa-mem2))
5. Variant calling with one or more of [`iVar`](https://andersen-lab.github.io/ivar/html/), [`FreeBayes`](https://github.com/freebayes/freebayes) and [`LoFreq`](https://csb5.github.io/lofreq/)

FASTA (consensus) input:

1. Align each per-gene consensus record to its reference gene region and call the differences (`bin/fasta_consensus_variants.py`, [`Biopython`](https://biopython.org/))

Both:

2. Annotate variants ([`SnpEff`](https://pcingola.github.io/SnpEff/) / [`SnpSift`](https://pcingola.github.io/SnpEff/snpsift/introduction/)) and match against the phenotype database
3. Per-sample and per-run resistance reports ([`Quarto`](https://quarto.org/), HTML + PDF)
4. Aggregate QC report ([`MultiQC`](http://multiqc.info/))

<img width="2000" height="1512" alt="Pipeline Tube Map-selection" src="https://github.com/user-attachments/assets/1d514f4e-d91d-406f-8e65-16de6a4ed8f1" />

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test,conda_singularity` before running the workflow on actual data.

> [!IMPORTANT]
> This pipeline currently needs **both Conda and a container engine** (Singularity or Docker) available on the machine that runs it, and is run with `-profile conda_singularity` (or `conda_docker`). Most tools run from containers, but a few steps (the Quarto reports) do not yet have a usable container and fall back to Conda. We hope to make it run on either engine alone in a future update.

This section is a quick start. See [docs/usage.md](docs/usage.md) for input types, the full parameter list, profiles and troubleshooting.

### Test the installation

Before running on your own data, check that Nextflow, Conda and your container engine work together using the bundled test dataset — a small simulated CMV amplicon dataset with known spiked-in resistance mutations:

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile test,conda_singularity \
   --outdir <TEST_OUTDIR>
```

This runs the whole workflow end to end and needs no reference downloads (host removal and Kraken2 are skipped). The per-sample report written to `<TEST_OUTDIR>/report_files/` should list the six spiked mutations: UL97 A594V, UL54 A809V, UL27 DEL218, UL51 P91S, UL56 C25F and UL89 N329S. Swap `conda_singularity` for `conda_docker` to test a Docker setup.

### Running on your own data

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row is a single fastq file (single-end) or a pair of fastq files (paired-end). The `bin/create_samplesheet.py` helper builds one from a directory of fastq files: `python bin/create_samplesheet.py /path/to/fastq_dir samplesheet.csv`.

To start from a per-gene consensus instead (e.g. Sanger sequencing), add `--input_type fasta` and use a two-column samplesheet (`sample,fasta`) where each `fasta` is one multi-record FASTA per sample, one record per resistance gene with the gene name at the start of the header. See [docs/usage.md](docs/usage.md#samplesheet-input-fasta--sanger-consensus).

Host-read removal ([`nohuman`](https://github.com/mbhall88/nohuman) by default, or [`Hostile`](https://github.com/bede/hostile) with `--host_removal_tool hostile`) and taxonomic classification ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/)) each need a database downloaded once — see the [`--nohuman_db`](docs/usage.md#--nohuman_db) (or [`--hostile_reference_path`](docs/usage.md#--hostile_reference_path)) and [`--kraken2_db`](docs/usage.md#--kraken2_db) sections of `docs/usage.md` for the download commands. To skip either step, pass `--skip_host_removal` and/or `--skip_kraken2` instead of the corresponding path.

Then run the pipeline using:

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile conda_singularity \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --protocol amplicon \
   --nohuman_db <path/to/nohuman_db/HPRC.r2> \
   --kraken2_db <path/to/kraken2_virus_db> \
   --report_author <your_name/institute> \
   --run_name <run_name>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

This pipeline is still a work in progress and any feedback or suggestions are welcome!

## Credits

CMVResistanceFinder was originally written by Jocelyn Hume.

We thank the following people for their assistance in the development of this pipeline:

Dr Emma Sweeney - the University of Queensland

Dr Kym Lowry - the University of Queensland

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
