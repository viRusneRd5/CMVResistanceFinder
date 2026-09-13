## FAQ's

**Why are there two different bed files?**

The first bed file (`ch_genes_bed` or `assets/CMV_NC_006273.2_resis_genes.bed`) describes the start and end positions of the six genes involved in CMV antiviral resistance (UL27, UL51, UL54, UL56, UL89 and UL97). Note that UL89 has a splice site in the middle, so there are two coding positions for this gene. This bed file should look like this:

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

The second bed file (`ch_resis_mut_bed` or `assets/resistance_sites.bed`) describes the start and end locations of known CMV antiviral resistant mutations. This is only used in the quarto report file to list and plot the locations with insufficient sequencing depth to call a mutation. This file must be recreated if you are not using the merlin reference sequence NC_006273.2 or if the phenotype database (`--cmv_resis_db`) is updated. When you run with a non-default reference the pipeline regenerates it for you from the database and GFF.

**What is the difference between *amplicon* and *wgs* protocol?**

The primary purpose of this pipeline is to identify mutations in the six genes associated with antiviral resistance, and provide characterisation data to contextualise them where available. However, some people will do whole-genome sequencing and want variants and coverage across the whole genome, so there are two protocols.

- **`amplicon`** (default): FreeBayes and LoFreq only call variants within the six resistance genes, and the reports only show variants and coverage for those genes. Faster, and appropriate whichever way the sample was prepared - most whole-genome runs still only need the resistance genes assessed.
- **`wgs`**: all three callers run across the whole genome, and the reports show genome-wide variants and a whole-genome coverage plot. Slower. Useful on whole-genome data, or on amplicon data to check for non-specific amplification.

`amplicon` can be used on whole-genome data and `wgs` can be used on amplicon data - the choice only controls how much of the genome is looked at, not how the reads were generated.

**Why doesn't the pipeline remove duplicate reads using picard?**

Picard `MarkDuplicates` identifies duplicates by fragment start/end coordinates. In amplicon data every read from an amplicon shares the same coordinates by design, so Picard would flag almost the whole amplicon as duplicate and collapse it to a single read, destroying the depth we need for variant calling.

Removing genuine PCR duplicates from amplicon data requires unique molecular identifiers (UMIs), which my assay doesn't use. Without them the pipeline assumes that alleles are amplified roughly equally, which is the assumption behind treating read counts as variant frequencies. One consequence is that a PCR error introduced in an early cycle is copied into many reads and can appear as a stable low-frequency variant - another reason to interpret calls near the frequency threshold with caution.

(For hybridisation-capture whole-genome data, where fragmentation is random and duplicates behave normally, deduplication would be appropriate. The `wgs` protocol does not currently do this.)

**Why is the default threshold's minimum 40 read depth, read quality 20, 0.05 (5%) variant frequency etc.?**

These are pragmatic defaults rather than values from a formal analysis.

**5% variant frequency** is a reporting threshold near the lower end of what's commonly used for CMV NGS resistance testing; there is no field consensus. Note that at 40x a 5% call is supported by only ~2 reads, so treat near-threshold, low-depth calls with caution - and be aware that low-frequency resistant variants can still be clinically significant, as they can expand under drug selection.

**40x depth** gives reasonable confidence for majority variant calls. **Base quality 20** (≈1% error) was also set pragmatically, as a permissive floor; a stricter value such as 30 would give more confidence in low-frequency calls but is not currently exposed as a parameter.

**Minimum mapping quality** (`--min_mapq`) defaults to 20 and is passed to all three callers, so reads that aligned ambiguously are dropped before calling. The six resistance genes sit in the unique-long region of the genome rather than in CMV's repeat regions, so this has little effect on coverage there; set `--min_mapq 0` to disable it.

You can change `--min_depth`, `--min_freq` and `--min_mapq` on the command line or in a params/config file.


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
      You don't have to build this by hand - `bin/create_resis_genes_bed.py` will generate it from your GFF file, as long as the genes are annotated with standard `gene=` attributes on their CDS features:
      ```bash
      python bin/create_resis_genes_bed.py <path/to/gff/file> --ref-name "<ref_name>"
      ```
      This writes `CMV_<ref_name>_resis_genes.bed` for you to pass to `--cmv_genes_bed`. It handles spliced genes (like UL89 in the Merlin annotation) automatically, splitting them into `_A`/`_B` entries, and will warn you if a gene is missing or looks unexpectedly split.

```bash
nextflow run viRusneRd5/cmvresistancefinder \
   -profile conda_singularity \
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

Yes. Run with `--input_type fasta` and give a two-column samplesheet (`sample,fasta`) where each `fasta` is one multi-record FASTA per sample, with one record per resistance gene and the gene name at the start of each header (e.g. `>UL54`). Each record is aligned to the reference and diffed; read-processing steps are skipped and the reports drop the coverage sections in favour of a "regions assessed" table. See the [usage docs](usage.md#samplesheet-input-fasta--sanger-consensus) for the full format and the FASTA-specific options.

**What if I want to provide my own database?**

If you want to use your own database, then as long as you follow the same format as the existing database then you shouldn't have any problems. The database used in this pipeline was compiled from a huge systematic review [hopefully can cite soon!]. If you are aware of some mutations that are missing from the database, please get in touch and I can add it into the next update.

**How long does the pipeline take to run?**

This depends on how large your fastq files are, and whether you are running this on a HPC or a laptop. Running 33 paired end Illumina fastq files sequenced on a MiSeq i100 5M 300-cycle kit took less than 10 minutes to run with 48GB memory and 24 threads and using cached conda/singularity containers.

**Why doesn't this pipeline perform primer trimming using IVAR trim?**

This pipeline is meant to be usable across many sequencing assays - different amplicon panels, hybridisation capture, whole genome - without asking each user to supply an assay-specific primer BED file and tune trimming settings. `ivar trim` needs exactly that, so it isn't built in.

For our own amplicon assay this is also fine scientifically: the primers were designed so that no primer-binding site overlaps (or sits within a read length of) a codon that is interrogated for resistance. Primer-derived sequence therefore can't mask a resistance variant or bias its allele frequency, because the reads covering a resistance codon are always template-derived.

If your primer scheme does place primers within a read length of a resistance codon, trim them externally (e.g. `ivar trim`, `samtools ampliconclip`, fgbio) before running this pipeline - note that this workflow has not been tested on externally primer-trimmed input.

**Can I run this pipeline with Ion Torrent (or other non-Illumina short-read) data?**

Structurally, yes - a samplesheet can freely mix paired-end and single-end rows (`--input_type short-read`), and each sample is auto-detected and processed accordingly. Single-end samples are flagged in their sample report ("Sequencing mode: Single-end") as a heads-up, since single-end short-read data in a samplesheet otherwise dominated by paired-end Illumina reads is often a sign of a different platform (e.g. the Vela Sentosa SQ CMV assay, which runs on Ion Torrent).

Scientifically, treat results from non-Illumina platforms with more caution, particularly for **indel calls** (deletions/insertions, reported as `DEL`/`INS` mutations - e.g. `UL27` deletions, `UL54` insertions): the default aligner (`bwa-mem3`) and all three variant callers (`iVar`, `FreeBayes`, `LoFreq`) were tuned and tested against Illumina (substitution-dominated error) and, to a lesser extent, Nanopore data. Ion Torrent's dominant error mode is homopolymer-associated indel slippage - exactly the variant class this pipeline reports as resistance mutations - and `LoFreq` in particular assumes an Illumina-style quality-score error model. Of the three callers, `FreeBayes` (haplotype-based) is the most platform-robust choice. We'd recommend manually checking any DEL/INS call from non-Illumina data against the aligned reads (e.g. in IGV) before relying on it, and are happy to hear from anyone who has validated this pipeline against Ion Torrent data.

**What is UL55/UL73 genotyping (`--genotype_gb` / `--genotype_gn`), and how is it different from the resistance mutations?**

It's strain typing, not resistance: gB (`UL55`) and gN (`UL73`) genotypes (e.g. `gB1`-`gB4`, `gN1`/`gN2`/`gN3a`/`gN3b`/`gN4a`-`gN4c`) are classic CMV epidemiological/strain markers used to track transmission and distinguish co-infecting strains - they don't indicate antiviral resistance, and the result is kept entirely separate from the resistance database and report.

The pipeline builds a consensus over the relevant gene (the same mechanism as `--export_gene_consensus`, just independent of that flag) and classifies it by nearest-neighbour percent identity against a small curated reference panel built from published reference sequences for each genotype. A worked example from our own reference: the Merlin genome's own (unmodified) `UL55`/`UL73` sequences classify confidently as `gB1` and `gN4d` respectively against these panels - Merlin is a real clinical/lab strain like any other, not a "neutral" or central reference for either gene, so this isn't surprising, it's just worth knowing when interpreting a WGS sample that happens to look similar to the reference.

Caveats: the reference panel is only as good as the genotypes it includes - a genuinely novel or unrepresented lineage will come back as `unresolved` rather than a wrong confident answer, which is the intended behaviour, not a bug. `gN4d` in particular is a rare lineage (previously reported mainly in the Merlin strain itself) rather than one of the "standard" seven gN genotypes most typing panels report - we included it because it's the closest match for Merlin-like sequences, but treat a `gN4d` call with the same "rare, check it" mindset you'd apply to any minority finding.
