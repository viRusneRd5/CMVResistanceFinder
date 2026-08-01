#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path
import re

# ------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Merge SnpSift outputs from multiple variant callers"
)

parser.add_argument(
    "--sample",
    required=True,
    help="Sample ID supplied by Nextflow"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output TSV filename"
)

parser.add_argument(
    "--min_freq",
    required=True,
    type=float,
    help="Minimum variant frequency"
)

parser.add_argument(
    "files",
    nargs="+",
    help="Input SnpSift files"
)

args = parser.parse_args()

sampleID = args.sample
output_file = args.output
input_files = args.files
min_freq = args.min_freq

# ------------------------------------------------------------------
# Read and process files
# ------------------------------------------------------------------

dfs = []

for file in input_files:

    df = pd.read_csv(file, sep="\t", dtype=str)

    # Equivalent to separate_longer_delim() on columns 5:ncol(df)
    explode_cols = list(df.columns[6:])

    for col in explode_cols:
        df[col] = (
            df[col]
            .fillna("")
            .astype(str)
            .str.split(",")
        )

    df = df.explode(explode_cols, ignore_index=True)

    # Extract variant caller from filename
    filename = Path(file).name

    m = re.search(r"_([^_]+)\.snpsift\.txt$", filename)

    if not m:
        raise ValueError(
            f"Could not determine variant caller from filename: {filename}"
        )

    variant_caller = m.group(1)

    df["sampleID"] = sampleID
    df["variant_caller"] = variant_caller

    dfs.append(df)

# ------------------------------------------------------------------
# Combine all callers
# ------------------------------------------------------------------

combined = pd.concat(dfs, ignore_index=True)

# Remove exact duplicate rows
combined = combined.drop_duplicates()

# ----------------------------------------------------
# Filter LoFreq variants with AF < 0.05
# ----------------------------------------------------

combined["AF[*]"] = pd.to_numeric(
    combined["AF[*]"],
    errors="coerce"
)

combined = combined[
    ~(
    (combined["variant_caller"].str.lower() == "lofreq") &
    (combined["AF[*]"] < min_freq)
    )
]

# ------------------------------------------------------------------
# Define columns that uniquely identify a variant
# ------------------------------------------------------------------
#
# Adjust these column names if needed.
#
# Common VCF-style columns:
# CHROM, POS, REF, ALT
#
# Example:
# variant_cols = ["sampleID", "CHROM", "POS", "REF", "ALT"]
#
# ------------------------------------------------------------------

variant_cols = ["sampleID", "CHROM", "POS", "REF", "ALT"]

missing = [x for x in variant_cols if x not in combined.columns]

if missing:
    raise ValueError(
        f"Missing required variant columns: {', '.join(missing)}\n"
        f"Available columns: {', '.join(combined.columns)}"
    )

# ------------------------------------------------------------------
# Summarise callers per variant
# ------------------------------------------------------------------

agg_dict = {
    "variant_caller": lambda x: ",".join(sorted(set(x)))
}

def collapse_caller_values(group, value_col):
    """
    Return caller:value pairs preserving caller order
    and removing duplicates caused by annotation explosion.
    """
    seen = set()
    values = []

    for caller, value in zip(group["variant_caller"], group[value_col]):

        if pd.isna(value):
            continue

        pair = (caller, str(value))

        if pair not in seen:
            seen.add(pair)
            values.append(f"{caller}:{value}")

    return ",".join(values)

# Standard aggregation
agg_dict = {
    "variant_caller": lambda x: ",".join(pd.unique(x))
}

for col in combined.columns:
    if col not in variant_cols + ["variant_caller", "DP[*]", "AF[*]"]:
        agg_dict[col] = "first"

combined.loc[
    combined["variant_caller"].str.lower() == "freebayes",
    "AF[*]"
] = pd.NA

merged = (
    combined
    .groupby(variant_cols, as_index=False)
    .agg(agg_dict)
)

# Create caller-specific DP and AF fields
dp_values = (
    combined.groupby(variant_cols, group_keys=False)
    .agg({"DP[*]": lambda x: collapse_caller_values(
        combined.loc[x.index],
        "DP[*]"
    )})
    .reset_index()
)

af_values = (
    combined.groupby(variant_cols, group_keys=False)
    .agg({"AF[*]": lambda x: collapse_caller_values(
        combined.loc[x.index],
        "AF[*]"
    )})
    .reset_index()
)

merged = merged.merge(dp_values, on=variant_cols)
merged = merged.merge(af_values, on=variant_cols)

merged = merged.rename(
    columns={"variant_caller": "detected_by"}
)

# ------------------------------------------------------------------
# Optional: count supporting callers
# ------------------------------------------------------------------

merged["n_callers"] = (
    merged["detected_by"]
    .str.split(",")
    .str.len()
)

# ------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------

merged.to_csv(
    output_file,
    sep="\t",
    index=False
)