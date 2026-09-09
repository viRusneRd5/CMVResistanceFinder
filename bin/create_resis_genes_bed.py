#!/usr/bin/env python
"""
Build a CMV_<accession>_resis_genes.bed file from a reference GFF.

The pipeline needs a BED file giving the start/end of the six genes involved
in CMV antiviral resistance (UL27, UL51, UL54, UL56, UL89, UL97 by default -
see docs/faq.md). If you're using a reference other than the bundled Merlin
sequence (NC_006273.2) you'd otherwise have to build that BED file by hand.
This script builds it for you from the reference's GFF file, so long as the
GFF annotates those genes with standard `gene=` attributes on CDS features.

UL89 is spliced (an intron in the middle) in the Merlin annotation, so it
gets two BED entries, UL89_A and UL89_B, one per CDS segment, in the order
they appear in the GFF. Any other target gene with more than one CDS segment
is treated the same way (gene_A, gene_B, ...) and flagged with a warning, in
case that's unexpected for your reference/annotation.

Example:
    python create_resis_genes_bed.py my_reference.gff --ref-name X17403.1
    -> writes CMV_X17403.1_resis_genes.bed
"""

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger()

DEFAULT_GENES = ["UL27", "UL51", "UL54", "UL56", "UL89", "UL97"]


def parse_attributes(field):
    """Parse a GFF3 column-9 attribute string into a dict, e.g. {'gene': 'UL97', ...}."""
    attrs = {}
    for part in field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, value = part.partition("=")
        attrs[key.strip()] = value.strip()
    return attrs


def extract_cds_by_gene(gff_file, target_genes):
    """
    Walk a GFF3 file and collect CDS features for the requested genes.

    Returns:
        dict: gene name (as requested, case-preserved) -> list of
            (seqid, start, end) tuples, in file order.
        set: seqids seen across the collected features (to sanity-check
            that everything comes from a single reference sequence).
    """
    wanted = {g.upper(): g for g in target_genes}
    hits = defaultdict(list)
    seqids = set()

    with gff_file.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _source, feature_type, start, end, _score, _strand, _frame, attr_field = fields[:9]
            if feature_type != "CDS":
                continue
            attrs = parse_attributes(attr_field)
            gene = attrs.get("gene")
            if not gene or gene.upper() not in wanted:
                continue
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                logger.warning(f"Line {line_no}: couldn't parse start/end ({start}, {end}) - skipping.")
                continue
            hits[wanted[gene.upper()]].append((seqid, start_i, end_i))
            seqids.add(seqid)

    return hits, seqids


def build_bed_rows(hits, target_genes):
    rows = []
    missing = []
    for gene in target_genes:
        segments = hits.get(gene, [])
        if not segments:
            missing.append(gene)
            continue
        if len(segments) == 1:
            seqid, start, end = segments[0]
            rows.append((seqid, start, end, gene))
        else:
            if len(segments) > 2:
                logger.warning(
                    f"{gene}: found {len(segments)} CDS segments (expected 1, or 2 for a spliced gene "
                    f"like UL89 in the Merlin annotation) - double-check these are all genuinely part "
                    f"of the same coding sequence."
                )
            suffixes = "ABCDEFGH"
            for i, (seqid, start, end) in enumerate(segments):
                rows.append((seqid, start, end, f"{gene}_{suffixes[i]}"))
    return rows, missing


def write_bed(rows, file_out):
    file_out.parent.mkdir(parents=True, exist_ok=True)
    with file_out.open("w") as handle:
        for seqid, start, end, name in rows:
            handle.write(f"{seqid}\t{start}\t{end}\t{name}\n")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("gff_file", type=Path, help="Reference GFF3 file (must have gene= attributes on CDS features).")
    parser.add_argument(
        "--ref-name",
        help="Reference/accession name to use in the output filename and BED chrom column check "
        "(must match the CHROM name in your reference fasta header). If omitted, it's inferred "
        "from the GFF, as long as all matched genes are on a single sequence.",
    )
    parser.add_argument(
        "--genes",
        default=",".join(DEFAULT_GENES),
        help=f"Comma-separated list of gene names to extract (default: {','.join(DEFAULT_GENES)}).",
    )
    parser.add_argument("--outdir", type=Path, default=Path("."), help="Directory to write the BED file to (default: current directory).")
    parser.add_argument(
        "-l", "--log-level", default="INFO", choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG")
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.gff_file.is_file():
        logger.critical(f"{args.gff_file} not found.")
        return 2

    target_genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    hits, seqids = extract_cds_by_gene(args.gff_file, target_genes)

    if not hits:
        logger.critical(
            f"None of the requested genes ({', '.join(target_genes)}) were found as CDS features with a "
            f"'gene=' attribute in {args.gff_file}. Check the GFF uses standard gene= tags."
        )
        return 1

    ref_name = args.ref_name
    if ref_name is None:
        if len(seqids) != 1:
            logger.critical(
                f"Matched genes span multiple sequences ({', '.join(sorted(seqids))}) - pass --ref-name "
                f"explicitly so the output can be built unambiguously."
            )
            return 1
        ref_name = next(iter(seqids))
    elif seqids and ref_name not in seqids:
        logger.warning(
            f"--ref-name '{ref_name}' doesn't match the seqid(s) the matched genes were found on "
            f"({', '.join(sorted(seqids))}). Make sure this matches the CHROM name in your reference fasta header."
        )

    rows, missing = build_bed_rows(hits, target_genes)
    if missing:
        logger.warning(f"No CDS found for: {', '.join(missing)}. They will be missing from the BED file.")

    file_out = args.outdir / f"CMV_{ref_name}_resis_genes.bed"
    write_bed(rows, file_out)
    logger.info(f"Wrote {file_out} with {len(rows)} entries for {len(target_genes) - len(missing)}/{len(target_genes)} genes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
