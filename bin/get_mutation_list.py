#!/usr/bin/env python3

import argparse
import pandas as pd

# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Extract CMV resistance mutations from a database."
)
parser.add_argument(
    "-i", "--input",
    required=True,
    help="Input CMV mutation database CSV"
)
parser.add_argument(
    "-o", "--output",
    default="resistant_mutations.txt",
    help="Output text file (default: resistant_mutations.txt)"
)
args = parser.parse_args()

# Read the database
df = pd.read_csv(args.input)

# Automatically find resistance columns
resistance_cols = [c for c in df.columns if c.endswith("_resistance")]

# Remove combination mutations
df = df[
    ~df["Gene"].str.contains("/", na=False) &
    ~df["Mutation"].str.contains("/", na=False)
]

# Keep only mutations with detected resistance
resistance_levels = {"Low", "Medium", "High"}

has_resistance = df[resistance_cols].isin(resistance_levels).any(axis=1)
df = df[has_resistance]

# Create gene_mutation strings
mutation_list = (
    df["Gene"] + "_" + df["Mutation"]
).drop_duplicates()

# Write output
mutation_list.to_csv(args.output, index=False, header=False)

print(f"Wrote {len(mutation_list)} resistant mutations to {args.output}")