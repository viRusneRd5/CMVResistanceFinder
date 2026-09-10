# viRusneRd5/cmvresistancefinder: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC (read input only)
- [Per-caller variant files](#per-caller-variant-files) - raw and SnpEff/SnpSift-annotated VCFs from each variant caller
- [Variant and resistance files](#variant-and-resistance-files) - per-sample merged variant, resistance and (fasta input) assessed-region tables
- [Reports](#reports) - per-sample and per-run HTML/PDF resistance reports
- [Gene consensus](#gene-consensus) - optional per-gene consensus FASTA export
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Per-caller variant files

<details markdown="1">
<summary>Output files (read input)</summary>

- `vcf_files/ivar/`
  - `<sample>_ivar.tsv`: raw iVar variant table.
  - `<sample>_ivar.vcf.gz` (+ `.csi`): the same calls converted to VCF.
  - `snpeff/<sample>_ivar.snpeff.vcf`: iVar calls annotated by SnpEff, plus `.snpeff.csv`, `.snpeff.genes.txt` and `.snpeff.summary.html` effect summaries.
  - `snpsift/<sample>_ivar.snpsift.txt`: the annotated fields pulled out into a flat table by SnpSift (the input to `merged.tsv`).
- `vcf_files/freebayes/`
  - `<sample>.vcf.gz` (+ `.csi`): coordinate-sorted FreeBayes calls.
  - `snpeff/<sample>_freebayes.snpeff.vcf` (+ `.csv`, `.genes.txt`, `.summary.html`) and `snpsift/<sample>_freebayes.snpsift.txt`.
- `vcf_files/lofreq/`
  - `<sample>_lofreq.vcf.gz` (+ `.csi`): LoFreq calls.
  - `snpeff/<sample>_lofreq.snpeff.vcf` (+ `.csv`, `.genes.txt`, `.summary.html`) and `snpsift/<sample>_lofreq.snpsift.txt`.

</details>

<details markdown="1">
<summary>Output files (fasta input)</summary>

- `vcf_files/sanger/`
  - `snpeff/<sample>_sanger.snpeff.vcf` (+ `.csv`, `.genes.txt`, `.summary.html`): the consensus-vs-reference differences, annotated by SnpEff.
  - `snpsift/<sample>_sanger.snpsift.txt`: the annotated fields as a flat table.

</details>

Which callers run depends on `--variant_caller` (default `all` = iVar + FreeBayes + LoFreq). Each caller's `snpsift.txt` is what feeds the merged variant table below; the raw and annotated VCFs are kept so a specific call can be traced back to the caller and evidence that produced it.

### Variant and resistance files

<details markdown="1">
<summary>Output files</summary>

- `variant_files/`
  - `<sample>.merged.tsv`: annotated variants merged across the selected variant callers (read input) or from the consensus differ (fasta input). One row per variant, with a `detected_by` column listing which caller(s) found it.
  - `<sample>.formatted.tsv`: the above with amino-acid changes normalised and synonymous variants dropped.
  - `<sample>.resistance.tsv`: the subset of variants that match a mutation in the phenotype database, one row per matching drug (and, for combination-mutation database entries, per entry). Carries the EC50 fold-change ratio, level of resistance and supporting references.
  - `<sample>.assessed.tsv` (**fasta input only**): one row per submitted gene record - the genome region it covered, alignment identity, whether it was reverse-complemented, and how many variants it yielded.

</details>

For fasta input `resistance.tsv` and `formatted.tsv` also carry `mixed_population` / `mixed_with_reference` columns, derived from IUPAC ambiguity codes in the consensus.

### Reports

<details markdown="1">
<summary>Output files</summary>

- `report_files/`
  - `<sample>.html` / `<sample>.pdf`: per-sample resistance report.
  - `<run_name>.html`: cross-sample run summary.

</details>

The per-sample report has the following sections:

- **Detected Mutations - Predicted Resistance** and **Key Mutations Characterised in the Literature** - the variants from `resistance.tsv`: a coloured drug-by-mutation grid of the predicted resistance level, then a table with the EC50 fold-change ratio, level of resistance and the supporting references for each mutation/drug combination.
- **Mutations of Unknown Significance** - non-synonymous variants found in the CMV resistance genes (`UL27`, `UL51`, `UL54`, `UL56`, `UL89`, `UL97`; or genome-wide with `--protocol wgs`) that are **not** listed in the phenotype database. These have no established resistance interpretation - they may be natural strain polymorphisms, sequencing artefacts, or genuinely novel mutations that have not been phenotypically characterised. The pipeline cannot classify them, so it just lists them (position, change, gene, amino-acid change, predicted effect, and - for read input - frequency, depth and which callers found them) for a reviewer to weigh up.
- **AMR Gene Coverage** (read input only) - every characterised resistance mutation has a known codon position (from `assets/resistance_sites.bed`). For each, the pipeline takes the mean read depth from the aligned BAM; where that is below `--min_depth` (default 40) the codon **could not be assessed** - a resistance mutation could be present there but there were not enough reads to call it. These codons are shaded red on the per-gene coverage plots and counted per gene in the "Codons below threshold" table. A report with no red means every known resistance locus had adequate coverage, so a "no resistance mutation detected" result can be trusted for those positions.

The run report drops the per-sample detail and shows a cross-sample resistance summary plus (read input) coverage/depth heatmaps or (fasta input) a per-sample, per-gene "records assessed" overview.

Fasta-input reports replace the coverage/depth sections with a "Regions Assessed" table (from `assessed.tsv`) - listing which gene regions a consensus record covered, so an absent gene reads as "not assessed" rather than "no mutation" - collapse same-codon mixed calls to a single notation (`M460V/I`, `A594V/A`) and add a "Mixed population" column.

### Gene consensus

<details markdown="1">
<summary>Output files</summary>

- `gene_consensus/` (only with `--export_gene_consensus`)
  - `<sample>.consensus.fasta`: one record per resistance gene, in reference coordinates and forward orientation, N-padded to the full gene region.

</details>

For read input this is built with `samtools consensus` (N below `--min_depth`, IUPAC codes above `--consensus_ambig_freq`); for `fasta` input it is the aligned, re-oriented submitted consensus.

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

Read QC, trimming, host removal, alignment and per-caller variant calling only run for read input (`--input_type short-read` or `long-read`), not for `--input_type fasta`.

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `software_versions.yml`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
