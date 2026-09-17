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
- `--min_base_quality` (default 20): minimum base quality to count a base towards a variant call, applied consistently across `iVar` (its own `-q` plus the feeding `samtools mpileup -Q`), `FreeBayes` (`--min-base-quality`) and `LoFreq` (`-q`/`--min-alt-bq`) - previously hardcoded to `20` in each caller with no way to change it.
- Bundled simulated CMV amplicon test dataset with known spiked-in resistance mutations, exercised by `-profile test` and nf-test.
- `conda_singularity` and `conda_docker` profiles.
- Optional strain genotyping (`--genotype_gb`, `--genotype_gn`): classifies `UL55` (gB) and/or `UL73` (gN) by nearest-neighbour percent identity against a curated reference panel, from either a full-gene (WGS) or partial (amplicon) consensus. Entirely separate from resistance-mutation calling. One combined `<outdir>/genotype/run_genotypes.tsv` per run (no per-sample files), plus a "Genotyping" section in the sample and run reports.
- `--genotype_min_depth` (default `10`): depth below which a genotyping-consensus position becomes `N`, independent of and much lower than `--min_depth` (used for resistance-mutation calling) - genotyping only needs one representative base per position, so low-titer/limit-of-detection samples that can't reach `--min_depth` can often still be genotyped.
- `confidence` column (`High`/`Medium`/`Low`, `N/A` for fasta input) on every variant in `merged.tsv`/`formatted.tsv`/`resistance.tsv`: flags calls detected by only one caller and/or with few estimated supporting reads, or sitting right at the `--min_depth`/`--min_freq` floor. Purely informational - nothing is filtered out, so sensitivity is unaffected. Highlighted (with a footnote-marked column header explaining the tiers) in the sample report's "Detected Mutations - Predicted Phenotype", "Key Mutations Characterised in the Literature" and "Mutations of Unknown Significance" tables, the equivalent run-report summary table, and a new run-report "Confidence overview" section (per-sample Low/Medium/High counts).
- Sample report now also shows `--min_freq`, `--min_base_quality` and `--min_mapq` alongside the existing depth threshold at the top of the report; run report shows `--min_depth`/`--min_base_quality` next to the AMR Gene Coverage plot.
- `--lofreq_max_depth` (default `10000`): caps LoFreq's own `-d`/`--max-depth`, which otherwise defaults to effectively unbounded (1000000) and can make an extremely deep site (e.g. a limit-of-detection dilution panel) take hours, since LoFreq's per-base significance test scales with depth. 10000x is comfortably above what's needed to call anything down to a few percent frequency, so this costs no practical sensitivity at typical `--min_freq` thresholds.

### `Fixed`

- gB genotyping's `gB3` reference (`M60930.1`) was a partial CDS (1002bp vs ~2721bp for gB1/gB2/gB4) covering a region that turned out to be nearly invariant between gB1/gB2/gB3 - unable to actually discriminate gB3 from the others, and a direct cause of unstable/spuriously-ambiguous or outright-wrong gB calls (confirmed against real production data spanning many repeat runs of the same sample). Replaced with `M85228.1` (complete CDS, strain C076A), the full-length gB3 representative used in a peer-reviewed real-time PCR gB-genotyping assay (Görzer et al.), verified to be genuinely divergent from gB1/gB2/gB4 (93-94% identity, consistent with their divergence from each other). Note: this also changes the Merlin reference genome's own UL55 classification from `gB3` to the now-correct `gB1`.
- Sample report's AMR Gene Coverage table listed and counted "affected codons" once per catalogued mutation rather than once per position, double-counting/duplicating codons where more than one mutation is catalogued at the same site (e.g. A594V and A594T).

### `Dependencies`

### `Deprecated`
