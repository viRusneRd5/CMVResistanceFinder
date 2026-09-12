# viRusneRd5/cmvresistancefinder: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of viRusneRd5/cmvresistancefinder, built with the [nf-core](https://nf-co.re/) pipeline template.

### `Added`

- Read (short-read / long-read) workflow: FastQC, adapter/quality trimming (`fastp`, `fastplong` or `Trimmomatic`), optional host-read removal (`Hostile`) and taxonomic classification (`Kraken2`), alignment (`minimap2`, `Bowtie2` or `BWA-MEM3`) and variant calling with any combination of `iVar`, `FreeBayes` and `LoFreq`.
- FASTA (Sanger / per-gene consensus) input workflow (`--input_type fasta`): aligns each per-gene consensus record to the reference and calls the differences, skipping all read processing and coverage-based report sections.
- Variant annotation with `SnpEff` / `SnpSift` and matching against a long-format meta-analysis phenotype database.
- Per-sample and per-run resistance reports (`Quarto`, HTML + PDF) and an aggregate QC report (`MultiQC`).
- `--protocol amplicon | wgs`: `amplicon` confines `FreeBayes`/`LoFreq` and the reports to the six resistance genes; `wgs` calls genome-wide with all callers.
- `--min_mapq` (default 20) applied consistently across `iVar`, `FreeBayes` and `LoFreq`; `LoFreq` minimum variant frequency applied via `lofreq filter`.
- Bundled simulated CMV amplicon test dataset with known spiked-in resistance mutations, exercised by `-profile test` and nf-test.
- `conda_singularity` and `conda_docker` profiles.
- Optional strain genotyping (`--genotype_gb`, `--genotype_gn`): classifies `UL55` (gB) and/or `UL73` (gN) by nearest-neighbour percent identity against a curated reference panel, from either a full-gene (WGS) or partial (amplicon) consensus. Entirely separate from resistance-mutation calling. One combined `<outdir>/genotype/run_genotypes.tsv` per run (no per-sample files), plus a "Genotyping" section in the sample and run reports.

### `Fixed`

### `Dependencies`

### `Deprecated`
