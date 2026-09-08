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

The second bed file (`ch_resis_mut_bed` or `assets/resistance_sites.bed`) describes the start and end locations of known CMV antiviral resistant mutations. This is only used in the quarto report file to list and plot the locations with insufficient sequencing depth to call a mutation. This file must be recreated if you are not using the merlin reference sequence NC_006273.2 or if the phenotyping_database is updated.

**What is the difference between *amplicon* and *wgs* protocol?**

This primary purpose of this pipeline is to identify mutations in the six genes associated with antiviral resistance and provide characterisation data to contextualise mutations where available. However, I do know that people may want to perform whole genome sequencing and get variants across all CDS and look at coverage across the genome. For this reason, I made the pipeline with two protocols. The amplicon protocol is not exclusively for use on amplicon sequencing data, this protocol can be used on whole genome sequencing data. The only difference between these protocols is that the reports will only output variants and coverage data from the six resistance genes, ignoring variants elsewhere in the genome and making analysis slightly faster. The wgs protocol can also be used on amplicon sequencing data, and could prove useful to identify non-specific amplification, but this analysis will take a little longer.

**Why doesn't the pipeline remove duplicate reads using picard?**

Sequencing CMV usually requires enrichment, whether this is probe hybridisation with PCR or just PCR. This means we need to duplicate the reads to be able to detect it with sequencing, especially when dealing with host DNA as well. Yes, we are assuming that all reads are amplified equally to get accurate variant frequencies, but if we were to remove duplicates, we may inadvertantly enhance minority variants and overstate their importance.

**Why is the default threshold's minimum 40 read depth, read quality 20, 0.05 (5%) variant frequency etc.?**

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
      You don't have to build this by hand - `bin/create_resis_genes_bed.py` will generate it from your GFF file, as long as the genes are annotated with standard `gene=` attributes on their CDS features:
      ```bash
      python bin/create_resis_genes_bed.py <path/to/gff/file> --ref-name "<ref_name>"
      ```
      This writes `CMV_<ref_name>_resis_genes.bed` for you to pass to `--cmv_genes_bed`. It handles spliced genes (like UL89 in the Merlin annotation) automatically, splitting them into `_A`/`_B` entries, and will warn you if a gene is missing or looks unexpectedly split.

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

If the demand is high enough, I might be able to make this work, but currently this is not possible.

**What if I want to provide my own database?**

This is complicated. The database used in this pipeline was compiled from a huge systematic review. Honestly, it is a bit of a mess, but if you start changing the database it will break the report writing. If you are aware of some mutations that are missing from the database, please get in touch and I can add it into the next update.

**How long does the pipeline take to run?**

This depends on how large your fastq files are, and whether you are running this on a HPC or a laptop. Running 33 paired end Illumina fastq files sequenced on a MiSeq i100 5M 300-cycle kit took less than 10 minutes to run with 48GB memory and 24 threads and using cached conda/singularity containers.

**Why do I have to use all three variant callers, why can't I choose which one to use?**

When developing this pipeline, I wanted to use multiple variant callers to be more confident in the variants that are called, and to make sure we aren't missing mutations. For this reason, the pipeline performs variant calling using all three variant callers and then combines the outputs of these.

**Why doesn't this pipeline perform primer trimming using IVAR trim?**

I looked into this option originally, considering my first focus was amplicon sequencing analysis, however I found that this removed too many reads and considerably dropped my sequencing coverage. My amplicon assay was designed outside of any regions containing AMR mutations and so it was decided not to include IVAR trim in the pipeline.