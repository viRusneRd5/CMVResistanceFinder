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
gets two BED entries, UL89_A and UL89_B. For UL89 specifically, the CDS
segment with the larger genomic coordinates is labelled UL89_A and the
segment with the smaller coordinates is labelled UL89_B, regardless of GFF
ordering.

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
    """Parse a GFF3 column-9 attribute string into a dict."""
    attrs = {}

    for part in field.strip().split(";"):
        part = part.strip()

        if not part or "=" not in part:
            continue

        key, _, value = part.partition("=")
        attrs[key.strip()] = value.strip()

    return attrs


def normalise_gene_name(gene):
    """
    Normalise CMV gene names across annotation conventions.

    Examples
    --------
    UL97      -> UL97
    HCMVUL97  -> UL97
    hcmvul54  -> UL54
    """
    gene = gene.upper()

    if gene.startswith("HCMVUL"):
        gene = gene.replace("HCMV", "", 1)

    return gene


def extract_cds_by_gene(gff_file, target_genes):
    """
    Walk a GFF3 file and collect CDS features for the requested genes.

    Returns
    -------
    dict
        gene name -> list of (seqid, start, end) tuples.
    set
        seqids represented by matched CDS features.
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

            (
                seqid,
                _source,
                feature_type,
                start,
                end,
                _score,
                _strand,
                _frame,
                attr_field,
            ) = fields[:9]

            if feature_type != "CDS":
                continue

            attrs = parse_attributes(attr_field)

            gene = attrs.get("gene")

            if not gene:
                continue

            norm_gene = normalise_gene_name(gene)

            if norm_gene not in wanted:
                continue

            try:
                start_i = int(start)
                end_i = int(end)

            except ValueError:
                logger.warning(
                    f"Line {line_no}: couldn't parse start/end "
                    f"({start}, {end}) - skipping."
                )
                continue

            hits[wanted[norm_gene]].append((seqid, start_i, end_i))
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

            if gene == "UL89" and len(segments) == 2:

                # UL89 is on the negative strand.
                # Assign UL89_A to the CDS segment with the larger
                # genomic coordinates and UL89_B to the smaller one.
                segments = sorted(
                    segments,
                    key=lambda x: max(x[1], x[2]),
                    reverse=True,
                )

                rows.append(
                    (segments[0][0], segments[0][1], segments[0][2], "UL89_A")
                )
                rows.append(
                    (segments[1][0], segments[1][1], segments[1][2], "UL89_B")
                )

            else:

                if len(segments) > 2:
                    logger.warning(
                        f"{gene}: found {len(segments)} CDS segments "
                        f"(expected 1, or 2 for a spliced gene like UL89) "
                        f"- double-check these are all genuinely part of the "
                        f"same coding sequence."
                    )

                suffixes = "ABCDEFGH"

                for i, (seqid, start, end) in enumerate(segments):
                    rows.append(
                        (seqid, start, end, f"{gene}_{suffixes[i]}")
                    )

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

    parser.add_argument(
        "gff_file",
        type=Path,
        help="Reference GFF3 file (must have gene= attributes on CDS features).",
    )

    parser.add_argument(
        "--ref-name",
        help=(
            "Reference/accession name to use in the output filename and BED "
            "chrom column check (must match the CHROM name in your reference "
            "fasta header). If omitted, it's inferred from the GFF, as long "
            "as all matched genes are on a single sequence."
        ),
    )

    parser.add_argument(
        "--genes",
        default=",".join(DEFAULT_GENES),
        help=f"Comma-separated list of gene names to extract "
             f"(default: {','.join(DEFAULT_GENES)}).",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("."),
        help="Directory to write the BED file to (default: current directory).",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        default="INFO",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
    )

    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)

    logging.basicConfig(
        level=args.log_level,
        format="[%(levelname)s] %(message)s",
    )

    if not args.gff_file.is_file():
        logger.critical(f"{args.gff_file} not found.")
        return 2

    target_genes = [
        g.strip()
        for g in args.genes.split(",")
        if g.strip()
    ]

    hits, seqids = extract_cds_by_gene(
        args.gff_file,
        target_genes,
    )

    if not hits:
        logger.critical(
            f"None of the requested genes ({', '.join(target_genes)}) "
            f"were found as CDS features with a 'gene=' attribute in "
            f"{args.gff_file}. Check the GFF uses standard gene= tags."
        )
        return 1

    ref_name = args.ref_name

    if ref_name is None:

        if len(seqids) != 1:
            logger.critical(
                f"Matched genes span multiple sequences "
                f"({', '.join(sorted(seqids))}) - pass --ref-name "
                f"explicitly so the output can be built unambiguously."
            )
            return 1

        ref_name = next(iter(seqids))

    elif seqids and ref_name not in seqids:
        logger.warning(
            f"--ref-name '{ref_name}' doesn't match the seqid(s) the "
            f"matched genes were found on ({', '.join(sorted(seqids))}). "
            f"Make sure this matches the CHROM name in your reference "
            f"fasta header."
        )

    rows, missing = build_bed_rows(hits, target_genes)

    if missing:
        logger.warning(
            f"No CDS found for: {', '.join(missing)}. "
            f"They will be missing from the BED file."
        )

    file_out = args.outdir / f"CMV_{ref_name}_resis_genes.bed"

    write_bed(rows, file_out)

    logger.info(
        f"Wrote {file_out} with {len(rows)} entries for "
        f"{len(target_genes) - len(missing)}/{len(target_genes)} genes."
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
