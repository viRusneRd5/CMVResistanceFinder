# viRusneRd5/cmvresistancefinder: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Generating a samplesheet automatically

If your FASTQ files are all in one directory, `bin/create_samplesheet.py` can build the samplesheet for you by pairing up R1/R2 files based on their filenames:

```bash
python bin/create_samplesheet.py /path/to/fastq_dir samplesheet.csv
```

It recognises the common Illumina naming pattern (`SAMPLE_S1_L001_R1_001.fastq.gz`), simpler `SAMPLE_R1/R2.fastq.gz` and `SAMPLE_1/2.fastq.gz` patterns, and treats any FASTQ with no R1/R2 marker as single-end. It warns about anything it couldn't confidently pair up rather than guessing - always check the output before running the pipeline with it.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile singularity,conda \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --protocol amplicon \
   --hostile_reference_path <path/to/hostile_reference> \
   --kraken2_db <path/to/kraken2_virus_db \
   --report_author <your_name/institute> \
   --run_name <run_name>
```

This will launch the pipeline with the `singularity,conda` configuration profile, which will run singularity or if not available then conda (working to fix so that it can run on either). See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run viRusneRd5/cmvresistancefinder -profile singularity,conda -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
protocol: 'amplicon'
hostile_reference_path: '<path/to/hostile_reference>'
kraken2_db: '<path/to/kraken2_virus_db'
run_name: '<run_name>'
report_author: '<report_author/institution>'
<...>
```
You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull viRusneRd5/cmvresistancefinder
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [viRusneRd5/cmvresistancefinder releases page](https://github.com/viRusneRd5/cmvresistancefinder/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### `--seq_type`

Specify if the samples in the run are either `'long-read'` (ONT) or `'short-read'` (Illumina). Default is `'short-read'`. If running long-read sequencing, the `--aligner` must be `minimap2`. 

**Note:** The long-read function has only been tested on two publicly available long-read samples, and the results were not super promising, so use with caution! If you have long-read sequencing data you would like to share with me to optimise this function, please reach out!

### `--trimmer`

Specify adaptor/quality trimming tool to be used. One of `'trimmomatic'` or `'fastp'`. Fastp is the default as it can auto-detect the adaptors to trim. If using trimmomatic, you will need to provide `--adapter_file <path/to/adapter.fasta>`. Note that if running on long-read sequences, this will default to fastplong for adaptor/quality trimming.

### `--skip_host_removal`

This allows you to skip performing host read removal using hostile. I would only recommend this if your pipeline already performs host-removal prior to running this pipeline, or if you are running public sequencing data where host-removal has presumably already occurred, but it doesn't hurt to leave this on.

### `--hostile_index`

The name of the hostile reference index used, e.g. `'human-t2t-hla'`.

### `--hostile_reference_path`

If you want to perform host-removal, you will need to download the hostile reference and index using bowtie2 build, as per below:

```bash
wget https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/human-t2t-hla.fa.gz

bowtie2-build /path/to/hostile_reference/human-t2t-hla.fa.gz /path/to/hostile_reference/hostile_bowtie_index
```
Then provide the `--hostile_reference_path /path/to/hostile_reference/hostile_bowtie_index` when running your pipeline. Ensure that the `--hostile_index` matches to the one you downloaded.

### `--kraken2_db`

The path to your Kraken2 database. If you wish to run Kraken2 to perform taxonomic classification on your reads prior to alignment to the reference, you will need to download the Kraken2 viral database. You can you this by:

```bash
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20260226.tar.gz
tar -xzvf /path/to/kraken_database/k2_viral_20260226.tar.gz
```
If you do not wish to perform kraken2 classification, add `--skip_kraken2 true` in place of `--kraken2_db </path/to/kraken2_virus_db>` below.

You can also use the Kraken2 standard database (up to date datases [found here](https://benlangmead.github.io/aws-indexes/k2)), but this database is larger and may take longer to run. Note that the pluspf database will not identify viruses, so should not be used here.

### `--kraken2_output`

Boolean to specify whether you want kraken2 to save the output fastq and reads assignment file. Default is `false`.

### `--skip_kraken2`

Boolean to skip Kraken2 read classification (default is `false`). If you do skip, it might save a little bit of time and doesn't effect the downstream of which reads are aligned to the reference sequence, will only change the multiqc output.

### `--variant_caller`

Currently, this must be set to `'all'` so that the pipeline will perform variant calling using freebayes, ivar and lofreq. If you change this setting to only run a single variant caller, this will break the compilation of variants later in the pipeline. Working to fix this.
<!-- TODO variant caller: Once this is fixed, I will need to change this. -->

### `--aligner`

You can choose between three aligning tools - bwa-mem3, bowtie2 and minimap2. The default is minimap2 as I found this gave me the best coverage, but you may have one that you prefer. If running longread sequencing, this must be set to minimap2.

### `--min_depth`

This sets the minimum read depth at the codon location required for calling a variant (default is `40`). This parameter sets the minimum depth for all three variant callers.

### `--min_freq`

This sets the minimum variant frequency required to call a variant (default is `0.05`, i.e. 5%). This parameter sets the minimum frequency flag for ivar and freebayes, while lofreq is filtered later in the pipeline when merging the data from all variant callers.

### `--skip_sample_report`

Boolean to skip the generation of the individual sample reports. Default is `false`.

## Extra custom configuration (for advanced users)

This pipeline includes the merlin CMV reference genome for aligning your reads to and calling variants. If you wish to use your own reference sequence, the pipeline can accommodate this, however you will need to prepare a few files. More info on this can be found in the [FAQ](faq.md).

### `--cmv_ref_fasta`

BYO reference sequence fasta. Read [FAQ](faq.md) for more detail.

### `--cmv_ref_gff`

BYO reference sequence gff. Read [FAQ](faq.md) for more detail.

### `--cmv_ref_name`

BYO reference sequence name (same has fasta file header). Read [FAQ](faq.md) for more detail.

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
