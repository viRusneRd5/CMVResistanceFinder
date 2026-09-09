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

def calculate_freebayes_af(ao, ro):
    """
    Calculate AF values for FreeBayes records.

    Handles:
      AO=12, RO=80          -> 0.1304
      AO=12,8, RO=80,80     -> 0.1304,0.0909
    """

    if pd.isna(ao) or pd.isna(ro):
        return None

    ao_vals = str(ao).split(",")
    ro_vals = str(ro).split(",")

    # If RO contains a single value, use it for all ALT alleles
    if len(ro_vals) == 1 and len(ao_vals) > 1:
        ro_vals = ro_vals * len(ao_vals)

    af_vals = []

    for ao_v, ro_v in zip(ao_vals, ro_vals):
        try:
            ao_num = float(ao_v)
            ro_num = float(ro_v)

            denom = ao_num + ro_num

            if denom == 0:
                af_vals.append("0")
            else:
                af_vals.append(f"{ao_num / denom:.6f}")

        except (ValueError, TypeError):
            af_vals.append("")

    return ",".join(af_vals)

def max_af_value(af):
    if pd.isna(af):
        return 0.0

    try:
        values = [
            float(x)
            for x in str(af).split(",")
            if x != ""
        ]
        return max(values) if values else 0.0
    except ValueError:
        return 0.0

for file in input_files:

    # Skip anything before the real header row. Concurrent JVM tasks (snpEff,
    # SnpSift) sharing /tmp can prepend a
    #   [warning][perf,memops] Cannot use file /tmp/hsperfdata_<user>/<pid> ...
    # line to the redirected output.
    with open(file) as handle:
        skiprows = next(
            (i for i, line in enumerate(handle) if line.startswith("CHROM\t")),
            0,
        )

    df = pd.read_csv(file, sep="\t", dtype=str, skiprows=skiprows)

    # Equivalent to separate_longer_delim() on columns 5:ncol(df)
    explode_cols = [c for c in df.columns if "[*]." in c]

    for col in explode_cols:
        df[col] = (
            df[col]
            .fillna("")
            .astype(str)
            .str.split(",")
        )

    if explode_cols:
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

## ------------------------------------------------------------------
## Combine all callers
## ------------------------------------------------------------------

combined = pd.concat(dfs, ignore_index=True)

## Remove exact duplicate rows

combined = combined.drop_duplicates()

# Recalculate AF for FreeBayes using AO/(AO+RO)

freebayes_mask = (
    combined["variant_caller"].str.lower() == "freebayes"
)

# ----------------------------------------------------
# Recalculate AF for FreeBayes
# AF = AO / (AO + RO)
# ----------------------------------------------------
#
# Only applies when FreeBayes rows are present. For inputs without FreeBayes
# (e.g. Sanger/FASTA consensus, which has no AO/RO columns) this is a no-op.

if freebayes_mask.any() and {"AO[*]", "RO[*]"}.issubset(combined.columns):
    combined.loc[freebayes_mask, "AF[*]"] = (
        combined.loc[freebayes_mask]
        .apply(
            lambda row: calculate_freebayes_af(
                row["AO[*]"],
                row["RO[*]"]
            ),
            axis=1
        )
    )

## ----------------------------------------------------
## Filter LoFreq variants with AF < 0.05
## ----------------------------------------------------

combined = combined[
    ~(
        (combined["variant_caller"].str.lower() == "lofreq") &
        (combined["AF[*]"].apply(max_af_value) < min_freq)
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