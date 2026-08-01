#!/usr/bin/env python3

import pandas as pd
import argparse


def mutation_match(query_mutation, database_mutation):
    """
    Returns True if:
    - exact match
    - query mutation is one component of a compound mutation
    """
    if query_mutation == database_mutation:
        return True

    db_mutations = str(database_mutation).split("/")
    return query_mutation in db_mutations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True, help="formatted.tsv file")
    parser.add_argument("--database", required=True, help="database.csv file")
    parser.add_argument("--output", required=True, help="output TSV")
    args = parser.parse_args()

    variants = pd.read_csv(args.variants, sep="\t")
    database = pd.read_csv(args.database)

    database_cols = [
        "Gene",
        "Mutation",
        "Ganciclovir_EC50_ratio",
        "ganciclovir_resistance",
        "Foscarnet_EC50_ratio",
        "foscarnet_resistance",
        "Cidofovir_EC50_ratio",
        "cidofovir_resistance",
        "Maribavir_EC50_ratio",
        "maribavir_resistance",
        "Letermovir_EC50_ratio",
        "letermovir_resistance",
        "TITLE",
        "Authors",
        "Published_Year",
        "DOI",
        "BIBTEXKEY"
    ]

    variant_cols = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "DP[*]",
        "AF[*]",
        "detected_by",
        "n_callers"
    ]

    available_variant_cols = [c for c in variant_cols if c in variants.columns]

    matches = []

    for _, variant in variants.iterrows():
        gene = variant["ANN[*].GENE"]
        mutation = variant["AA_change"]

        gene_db = database[database["Gene"] == gene]

        for _, record in gene_db.iterrows():
            if mutation_match(mutation, record["Mutation"]):
                row = record[database_cols].to_dict()

                for col in available_variant_cols:
                    row[col] = variant[col]

                matches.append(row)

    output_cols = available_variant_cols + database_cols

    if matches:
        result = pd.DataFrame(matches)
        result = result.drop_duplicates()
    else:
        result = pd.DataFrame(columns=output_cols)

    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
