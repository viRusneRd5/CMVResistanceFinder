#!/usr/bin/env python
"""
Build a viRusneRd5/cmvresistancefinder samplesheet.csv from a directory of
FASTQ and/or FASTA files.

Scans a directory for *.fastq.gz / *.fq.gz files, pairs up R1/R2 files that
belong to the same sample, and writes out a samplesheet in the
`sample,fastq_1,fastq_2` format the pipeline expects for `--input_type
short-read`/`long-read`. Single-end files (e.g. long-read data, or
already-merged fastqs) are written with an empty fastq_2 column.

It also scans for *.fa / *.fasta / *.fna files (per-gene consensus FASTA,
e.g. Sanger sequencing) and writes a separate `sample,fasta` samplesheet for
`--input_type fasta` - one row per file, sample name = filename stem. FASTQ
and FASTA samples can't be mixed into a single pipeline run, so if a
directory contains both, two samplesheets are written (see --file-out-fasta).

Recognised FASTQ naming conventions (checked in this order for each file):
    SAMPLE_S1_L001_R1_001.fastq.gz   (Illumina bcl2fastq / basespace)
    SAMPLE_R1.fastq.gz               SAMPLE_R2.fastq.gz
    SAMPLE_1.fastq.gz                SAMPLE_2.fastq.gz
    SAMPLE.fastq.gz                  (no read marker -> treated as single-end)

Example:
    python create_samplesheet.py /path/to/fastq_dir samplesheet.csv
    python create_samplesheet.py /path/to/fastq_dir samplesheet.csv --recursive
    python create_samplesheet.py /path/to/fasta_dir samplesheet_fasta.csv
"""

import argparse
import csv
import logging
import re
import sys
from pathlib import Path

logger = logging.getLogger()

VALID_SUFFIXES = (".fastq.gz", ".fq.gz")
FASTA_SUFFIXES = (".fasta", ".fna", ".fa")

# Tried in order; the first pattern that matches a filename wins.
# Each must capture "sample" and, for paired files, "read" as "1" or "2".
READ_PATTERNS = [
    # SAMPLE_S1_L001_R1_001.fastq.gz / SAMPLE_S1_L001_R2_001.fastq.gz
    re.compile(r"^(?P<sample>.+?)_S\d+_L\d+_R(?P<read>[12])_\d+$"),
    # SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz  (optionally with a trailing _001 etc.)
    re.compile(r"^(?P<sample>.+?)_R(?P<read>[12])(_\d+)?$"),
    # SAMPLE_1.fastq.gz / SAMPLE_2.fastq.gz
    re.compile(r"^(?P<sample>.+?)[._-](?P<read>[12])$"),
]


def strip_fastq_suffix(filename):
    """Return (stem, matched_suffix) for a fastq filename, or (None, None) if unrecognised."""
    for suffix in VALID_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)], suffix
    return None, None


def classify(stem):
    """Return (sample, read_number_or_None) for a fastq filename stem."""
    for pattern in READ_PATTERNS:
        match = pattern.match(stem)
        if match:
            return match.group("sample"), match.group("read")
    # No recognised R1/R2 marker -> treat the whole stem as the sample name,
    # to be written as a single-end entry.
    return stem, None


def find_fastqs(input_dir, recursive):
    globber = input_dir.rglob if recursive else input_dir.glob
    files = sorted(p for p in globber("*") if p.is_file() and p.name.endswith(VALID_SUFFIXES))
    return files


def strip_fasta_suffix(filename):
    """Return (stem, matched_suffix) for a fasta filename, or (None, None) if unrecognised."""
    for suffix in FASTA_SUFFIXES:
        if filename.endswith(suffix):
            return filename[: -len(suffix)], suffix
    return None, None


def find_fastas(input_dir, recursive):
    globber = input_dir.rglob if recursive else input_dir.glob
    files = sorted(p for p in globber("*") if p.is_file() and p.name.endswith(FASTA_SUFFIXES))
    return files


def build_fasta_samplesheet(fastas):
    """One row per FASTA file - sample name = filename stem. Returns (rows, warnings)."""
    samples = {}
    warnings = []
    for path in fastas:
        stem, _ = strip_fasta_suffix(path.name)
        sample = stem.replace(" ", "_")
        if sample in samples:
            warnings.append(
                f"{sample}: {path.name} also stems to this sample name - keeping "
                f"{samples[sample].name}, NOT adding {path.name} - rename one of them if both are real samples."
            )
            continue
        samples[sample] = path

    rows = [{"sample": sample, "fasta": str(path.resolve())} for sample, path in sorted(samples.items())]
    return rows, warnings


def build_samplesheet(fastqs):
    """Group fastq files into samples. Returns (rows, warnings)."""
    # sample -> {"1": path, "2": path, "unpaired": [paths]}
    samples = {}
    unrecognised = []
    for path in fastqs:
        stem, suffix = strip_fastq_suffix(path.name)
        if stem is None:
            unrecognised.append(path)
            continue
        sample, read = classify(stem)
        sample = sample.replace(" ", "_")
        entry = samples.setdefault(sample, {"1": None, "2": None, "unpaired": []})
        if read == "1":
            if entry["1"] is not None:
                entry["unpaired"].append(path)
            else:
                entry["1"] = path
        elif read == "2":
            if entry["2"] is not None:
                entry["unpaired"].append(path)
            else:
                entry["2"] = path
        else:
            entry["unpaired"].append(path)

    rows = []
    warnings = []
    for sample in sorted(samples):
        entry = samples[sample]
        r1, r2 = entry["1"], entry["2"]

        if r1 is None and r2 is not None:
            # Only an R2 was found on its own - treat it as a single-end file
            # rather than silently dropping it.
            r1, r2 = r2, None
            warnings.append(f"{sample}: found an R2 file with no matching R1; using it as a single-end read.")

        for extra in entry["unpaired"]:
            if r1 is None:
                r1 = extra
            else:
                warnings.append(
                    f"{sample}: extra file {extra.name} did not pair with anything and was NOT added "
                    f"to the samplesheet - if this is a real sample, add it manually."
                )

        rows.append(
            {
                "sample": sample,
                "fastq_1": str(r1.resolve()) if r1 else "",
                "fastq_2": str(r2.resolve()) if r2 else "",
            }
        )

    for path in unrecognised:
        warnings.append(f"Skipped {path.name}: doesn't end in {' or '.join(VALID_SUFFIXES)}.")

    return rows, warnings


def write_samplesheet(rows, file_out):
    file_out.parent.mkdir(parents=True, exist_ok=True)
    with file_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_1", "fastq_2"])
        writer.writeheader()
        writer.writerows(rows)


def write_fasta_samplesheet(rows, file_out):
    file_out.parent.mkdir(parents=True, exist_ok=True)
    with file_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fasta"])
        writer.writeheader()
        writer.writerows(rows)


def default_fasta_out(file_out):
    """samplesheet.csv -> samplesheet_fasta.csv (preserves any other name/suffix)."""
    return file_out.with_name(f"{file_out.stem}_fasta{file_out.suffix}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing the FASTQ and/or FASTA files.",
    )
    parser.add_argument(
        "file_out",
        type=Path,
        nargs="?",
        default=Path("samplesheet.csv"),
        help="Where to write the FASTQ samplesheet (default: samplesheet.csv in the current directory).",
    )
    parser.add_argument(
        "--file-out-fasta",
        type=Path,
        default=None,
        help="Where to write the FASTA samplesheet, if any *.fa/*.fasta/*.fna files are found "
        "(default: file_out with '_fasta' inserted, e.g. samplesheet_fasta.csv).",
    )
    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="Also look in subdirectories of input_dir.",
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
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.input_dir.is_dir():
        logger.critical(f"{args.input_dir} is not a directory.")
        return 2

    fastqs = find_fastqs(args.input_dir, args.recursive)
    fastas = find_fastas(args.input_dir, args.recursive)

    if not fastqs and not fastas:
        logger.critical(
            f"No *.fastq.gz / *.fq.gz / *.fa / *.fasta / *.fna files found in {args.input_dir}"
            + (" (recursively)" if args.recursive else "")
        )
        return 1

    if fastqs:
        rows, warnings = build_samplesheet(fastqs)
        write_samplesheet(rows, args.file_out)
        for w in warnings:
            logger.warning(w)
        paired = sum(1 for r in rows if r["fastq_2"])
        single = len(rows) - paired
        logger.info(f"Wrote {args.file_out} with {len(rows)} FASTQ sample(s): {paired} paired-end, {single} single-end.")

    if fastas:
        fasta_out = args.file_out_fasta or default_fasta_out(args.file_out)
        rows, warnings = build_fasta_samplesheet(fastas)
        write_fasta_samplesheet(rows, fasta_out)
        for w in warnings:
            logger.warning(w)
        logger.info(f"Wrote {fasta_out} with {len(rows)} FASTA sample(s) (--input_type fasta).")

    if fastqs and fastas:
        logger.warning(
            f"Found both FASTQ and FASTA files in {args.input_dir} - these can't be mixed into one pipeline "
            f"run, so two separate samplesheets were written above. Run the pipeline once per samplesheet, "
            f"with --input_type short-read/long-read for {args.file_out} and --input_type fasta for {fasta_out}."
        )

    logger.info("Double-check the samplesheet(s) before using them - filename-based sample detection can get it wrong.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
