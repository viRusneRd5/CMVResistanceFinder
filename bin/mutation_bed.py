#!/usr/bin/env python3

import re
import argparse


def parse_gff(gff_file):
    """
    Parse CDS features from a GFF3 file.

    Returns:
        dict:
            gene -> {
                "chrom": chromosome,
                "strand": "+"|"-",
                "cds_coords": [genomic positions in CDS order]
            }
    """

    cds_by_gene = {}

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")

            if len(fields) != 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attrs = fields

            if feature_type != "CDS":
                continue

            start = int(start)
            end = int(end)

            attr_dict = {}

            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k] = v

            gene = (
                attr_dict.get("gene")
                or attr_dict.get("gene_name")
                or attr_dict.get("Name")
                or attr_dict.get("locus_tag")
                or attr_dict.get("Parent")
            )

            if gene is None:
                continue

            # Normalize names such as HCMVUL97 -> UL97
            m = re.search(r'(UL\d+)', gene, re.IGNORECASE)
            if m:
                gene = m.group(1).upper()


            if gene not in cds_by_gene:
                cds_by_gene[gene] = {
                    "chrom": chrom,
                    "strand": strand,
                    "segments": []
                }

            cds_by_gene[gene]["segments"].append((start, end))

    genes = {}

    for gene, info in cds_by_gene.items():

        segments = info["segments"]

        if info["strand"] == "+":
            segments.sort(key=lambda x: x[0])

            cds_coords = []

            for start, end in segments:
                cds_coords.extend(range(start, end + 1))

        else:
            segments.sort(key=lambda x: x[0], reverse=True)

            cds_coords = []

            for start, end in segments:
                cds_coords.extend(range(end, start - 1, -1))

        genes[gene] = {
            "chrom": info["chrom"],
            "strand": info["strand"],
            "cds_coords": cds_coords
        }

    return genes


def aa_to_nt_range(gene_info, aa_position):
    """
    Convert amino-acid position to genomic nucleotide range.

    Returns:
        (start_nt, end_nt)
    """

    cds = gene_info["cds_coords"]

    codon_start = (aa_position - 1) * 3

    if codon_start + 2 >= len(cds):
        raise ValueError(
            f"AA position {aa_position} exceeds CDS length"
        )

    nts = cds[codon_start:codon_start + 3]

    return min(nts), max(nts)


def mutation_to_bed(genes, mutation):

    m = re.match(r"(.+)_DEL(\d+)-(\d+)$", mutation)

    if m:
        gene, aa1, aa2 = m.groups()

        aa1 = int(aa1)
        aa2 = int(aa2)

        start1, end1 = aa_to_nt_range(genes[gene], aa1)
        start2, end2 = aa_to_nt_range(genes[gene], aa2)

        start = min(start1, start2)
        end = max(end1, end2)

    elif m := re.match(r"(.+)_DEL(\d+)$", mutation):

        gene, aa = m.groups()

        start, end = aa_to_nt_range(
            genes[gene],
            int(aa)
        )

    elif m := re.match(
        r"(.+)_([A-Z\*])(\d+)([A-Z\*])$",
        mutation
    ):

        gene, ref, aa, alt = m.groups()

        start, end = aa_to_nt_range(
            genes[gene],
            int(aa)
        )

    else:
        raise ValueError(
            f"Unrecognized mutation format: {mutation}"
        )

    chrom = genes[gene]["chrom"]

    return (
        chrom,
        start - 1,  # BED 0-based
        end,
        mutation
    )


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--gff",
        required=True
    )

    parser.add_argument(
        "--mutations",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )

    args = parser.parse_args()

    genes = parse_gff(args.gff)

    with open(args.mutations) as fin, \
            open(args.output, "w") as fout:

        for line in fin:

            mutation = line.strip()

            if not mutation:
                continue

            try:
                chrom, start, end, name = (
                    mutation_to_bed(
                        genes,
                        mutation
                    )
                )

                fout.write(
                    f"{chrom}\t{start}\t{end}\t{name}\n"
                )

            except Exception as e:
                print(
                    f"Skipping {mutation}: {e}"
                )


if __name__ == "__main__":
    main()
