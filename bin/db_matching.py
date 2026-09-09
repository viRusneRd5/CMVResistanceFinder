#!/usr/bin/env python3

import argparse
import re

import pandas as pd

# Drug columns in the phenotype database, in report display order
DRUGS = ["Ganciclovir", "Foscarnet", "Cidofovir", "Maribavir", "Letermovir"]

# Columns copied straight from the variant (formatted.tsv) into the output
VARIANT_COLS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "DP[*]",
    "AF[*]",
    "detected_by",
    "n_callers",
    "mixed_population",
    "mixed_with_reference",
]

# Columns produced from the database for every variant/drug match
DB_OUTPUT_COLS = [
    "Gene",
    "Mutation",
    "Gene_Mutation",
    "match_type",
    "Drug",
    "EC50_ratio",
    "Level_of_resistance",
    "References",
]


def extract_gene(raw_gene):
    """Return the UL gene name (e.g. 'UL54') if present, else the raw value upper-cased."""
    raw_gene = str(raw_gene)
    match = re.search(r"(UL\d+)", raw_gene, re.IGNORECASE)
    return match.group(1).upper() if match else raw_gene.upper()


def split_components(gene_mutation):
    """
    Split a database 'Gene_Mutation' string into a set of (gene, mutation) tuples.

    Single mutation:   'UL54_A594T'                  -> {('UL54', 'A594T')}
    Compound mutation: 'UL54_A594T/UL54_V715M'       -> {('UL54', 'A594T'), ('UL54', 'V715M')}
    Cross-gene:        'UL27_R233S/UL97_L337M'       -> {('UL27', 'R233S'), ('UL97', 'L337M')}
    """
    components = set()
    for part in str(gene_mutation).split("/"):
        part = part.strip()
        if "_" not in part:
            continue
        gene, mutation = part.split("_", 1)
        components.add((gene.upper(), mutation))
    return components


def variant_components(gene, aa_change):
    """
    (gene, aa_change) from a variant -> set of (gene, mutation) tuples.

    A compound AA change such as 'V812L/T813S' is broken into its parts so that it
    can match individual database entries.
    """
    components = set()
    for part in str(aa_change).split("/"):
        part = part.strip()
        if part:
            components.add((gene, part))
    return components


def round_ec50(value):
    """Round a numeric EC50 ratio to 2 dp; pass non-numeric values (e.g. 'Susceptible') through."""
    try:
        return round(float(value), 2)
    except (TypeError, ValueError):
        text = str(value).strip()
        return "" if text.upper() in ("", "NA", "NAN") else text


def collect_references(record):
    """Merge bibtex_list_meta and bibtex_list_geo into a unique, ', '-joined key list."""
    keys = []
    for column in ("bibtex_list_meta", "bibtex_list_geo"):
        raw = record.get(column)
        if raw is None or (isinstance(raw, float) and pd.isna(raw)):
            continue
        for key in str(raw).split(","):
            key = key.strip()
            if key and key.upper() != "NA" and key not in keys:
                keys.append(key)
    return ", ".join(keys)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True, help="formatted.tsv file")
    parser.add_argument("--database", required=True, help="phenotype database CSV")
    parser.add_argument("--output", required=True, help="output TSV")
    args = parser.parse_args()

    variants = pd.read_csv(args.variants, sep="\t")
    database = pd.read_csv(args.database, dtype=str)

    # FASTA (Sanger consensus) input carries MIXED / MIXED_REF flags from the
    # consensus differ: MIXED means the base was an IUPAC ambiguity code with a
    # non-reference allele; MIXED_REF means the reference base was also part of
    # that ambiguity set (i.e. wild-type still present in the population).
    # Surface these as friendly Yes/No columns so the report can render mixed
    # calls (e.g. M460V/I, A594V/A) and flag them.
    def _yes_no(series):
        return series.astype(str).str.strip().str.lower().map(
            {"true": "Yes", "false": "No"}
        ).fillna("No")

    if "MIXED" in variants.columns:
        variants["mixed_population"] = _yes_no(variants["MIXED"])
    else:
        variants["mixed_population"] = "No"

    if "MIXED_REF" in variants.columns:
        variants["mixed_with_reference"] = _yes_no(variants["MIXED_REF"])
    else:
        variants["mixed_with_reference"] = "No"

    # Pre-parse the database into (row, {(gene, mutation), ...}) once
    db_records = [
        (record, split_components(record["Gene_Mutation"]))
        for _, record in database.iterrows()
    ]

    available_variant_cols = [c for c in VARIANT_COLS if c in variants.columns]
    output_cols = available_variant_cols + DB_OUTPUT_COLS

    matches = []
    for _, variant in variants.iterrows():
        aa_change = str(variant.get("AA_change", "")).strip()
        if aa_change in ("", ".", "nan"):
            continue

        gene = extract_gene(variant["ANN[*].GENE"])
        query = variant_components(gene, aa_change)

        for record, db_components in db_records:
            if query.isdisjoint(db_components):
                continue

            row = {col: variant[col] for col in available_variant_cols}
            row["Gene"] = variant["ANN[*].GENE"]
            row["Mutation"] = aa_change
            row["Gene_Mutation"] = record["Gene_Mutation"]
            row["match_type"] = "exact" if len(db_components) == 1 else "compound"
            row["Drug"] = record["Drug"]
            row["EC50_ratio"] = round_ec50(record["Overall_EC50_Ratio"])
            row["Level_of_resistance"] = record["level_of_resistance"]
            row["References"] = collect_references(record)
            matches.append(row)

    if matches:
        result = pd.DataFrame(matches)[output_cols].drop_duplicates()
    else:
        result = pd.DataFrame(columns=output_cols)

    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
