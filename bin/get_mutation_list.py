#!/usr/bin/env python3

import argparse

import pandas as pd

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Extract CMV resistance mutations from the phenotype database."
)
parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input CMV phenotype database CSV"
)
parser.add_argument(
    "-o", "--output",
    default="resistant_mutations.txt",
    help="Output text file (default: resistant_mutations.txt)"
)
args = parser.parse_args()

# Read the database. It is in long format with one row per Gene_Mutation / Drug
# combination and a single Level_of_resistance value per row.
df = pd.read_csv(args.input, dtype=str)

# Drop combination (compound) mutations - these can't be placed on the genome
# as a single site.
df = df[~df["Gene_Mutation"].str.contains("/", na=False)]

# Keep only mutations with reported resistance to at least one drug.
resistance_levels = {"Low", "Medium", "High"}
df = df[df["level_of_resistance"].isin(resistance_levels)]

# Gene_Mutation is already in "GENE_MUTATION" form (e.g. UL54_A594T).
mutation_list = df["Gene_Mutation"].drop_duplicates()

# Write output
mutation_list.to_csv(args.output, index=False, header=False)

print(f"Wrote {len(mutation_list)} resistant mutations to {args.output}")
