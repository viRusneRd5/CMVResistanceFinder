#!/usr/bin/env python3
"""
genotype_classify.py

Classify a sample's consensus sequence for one gene against a curated
reference panel of known genotypes/subtypes, by nearest-neighbour percent
identity. Gene-agnostic: works for any gene given a consensus FASTA record,
its strand, and a reference panel FASTA (built for UL55/gB initially; UL73/gN
follows the same interface).

Input consensus FASTA is expected to come from the pipeline's GENE_CONSENSUS
process (samtools consensus over a gene's genomic region): forward-strand
relative to the reference genome, '-' for below-depth positions, and IUPAC
ambiguity codes for minor alleles above --consensus_ambig_freq. For a
minus-strand gene this script reverse-complements before comparing to the
(mRNA-sense) reference panel.

Outputs one TSV row (with header) for this sample/gene, to be concatenated
across samples into the run-level genotype summary - this script never writes
more than one gene's result, and the caller (Nextflow module) runs it once
per enabled gene per sample.
"""

import argparse
import re
import sys

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# IUPAC ambiguity code -> the set of bases it represents. A comparison
# between a sample base and a reference base counts as a match if the
# reference's base is a member of the sample's ambiguity set (and vice
# versa for a fully-ambiguous reference, though the panel is expected to be
# unambiguous).
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT",
}


def read_one_record(fasta_path, record_name):
    """Return the sequence (str, uppercase) of the record named record_name, or None."""
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if rec.id == record_name or rec.description.split()[0] == record_name:
            return str(rec.seq).upper()
    return None


def read_panel(fasta_path):
    """Return {genotype_label: sequence} from a reference panel FASTA (header = genotype)."""
    panel = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        panel[rec.id] = str(rec.seq).upper()
    return panel


def bases_match(sample_base, ref_base):
    """True if two (possibly IUPAC-ambiguous) bases are consistent with each other."""
    s = IUPAC.get(sample_base, "")
    r = IUPAC.get(ref_base, "")
    return bool(set(s) & set(r))


def region_coverage(seq):
    """Fraction and count of non-gap, non-N (informative) bases in a consensus sequence."""
    informative = sum(1 for b in seq if b not in ("N", "-", "*"))
    return informative, (informative / len(seq) if seq else 0.0)


def compare_to_reference(sample_seq, ref_seq, aligner):
    """
    Align sample_seq to ref_seq and return (pct_identity, n_compared) over
    positions where both sides have real (non-N/non-gap) data. IUPAC
    ambiguity codes in the sample are treated as a match if consistent with
    the reference's (unambiguous) base.
    """
    alignment = aligner.align(sample_seq, ref_seq)[0]
    aligned_sample, aligned_ref = str(alignment[0]), str(alignment[1])

    n_compared = 0
    n_match = 0
    for s_base, r_base in zip(aligned_sample, aligned_ref):
        if s_base in ("N", "-", "*") or r_base in ("-", "*"):
            continue
        n_compared += 1
        if bases_match(s_base, r_base):
            n_match += 1

    pct_identity = (100.0 * n_match / n_compared) if n_compared else 0.0
    return pct_identity, n_compared


def detect_mixed_signal(sample_seq, ref_seq, aligner, best_label, panel, best_seq):
    """
    Look for IUPAC ambiguity codes in the sample that are consistent with a
    *different* panel genotype's base at the same position (not just noise) -
    a signature of a mixed infection rather than plain sequencing ambiguity.
    Returns a sorted list of other genotype labels with such support.
    """
    alignment = aligner.align(sample_seq, best_seq)[0]
    aligned_sample, aligned_best = str(alignment[0]), str(alignment[1])

    # Map alignment columns back to best_seq positions so we can look up the
    # same columns in the other panel members via their own alignment to best_seq.
    other_support = {}
    for label, ref_seq_other in panel.items():
        if label == best_label:
            continue
        alignment2 = aligner.align(best_seq, ref_seq_other)[0]
        aligned_best2, aligned_other = str(alignment2[0]), str(alignment2[1])
        # Build best_seq-column -> other-genotype-base lookup
        other_at_best_col = {}
        best_col = 0
        for b_base, o_base in zip(aligned_best2, aligned_other):
            if b_base != "-":
                other_at_best_col[best_col] = o_base
                best_col += 1

        best_col = 0
        hits = 0
        for s_base, b_base in zip(aligned_sample, aligned_best):
            if b_base != "-":
                if s_base not in ("A", "C", "G", "T", "N", "-", "*"):
                    other_base = other_at_best_col.get(best_col)
                    if other_base and other_base != b_base and bases_match(s_base, other_base):
                        hits += 1
                best_col += 1
        if hits >= 4:  # several supporting ambiguous sites, not just one coincidental hit
            other_support[label] = hits

    # Only the single best-supported alternate - reporting every panel member
    # that shares a handful of divergent alleles with the true minor strain
    # (which is expected between related genotypes) is noise, not signal.
    if not other_support:
        return []
    return [max(other_support, key=other_support.get)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--gene", required=True, help="Gene name (must match the record header in --consensus)")
    parser.add_argument("--strand", required=True, choices=["+", "-"], help="Gene strand on the reference genome")
    parser.add_argument("--consensus", required=True, help="Per-gene consensus FASTA (GENE_CONSENSUS output)")
    parser.add_argument("--reference_panel", required=True, help="FASTA of per-genotype reference sequences (header = genotype label)")
    parser.add_argument("--min_region_covered", type=float, default=0.05, help="Minimum fraction of the gene with non-N consensus before even attempting a call (rejects an essentially-empty consensus; NOT a 'whole gene' requirement - an amplicon assay covering a small, informative fraction of a large gene is expected and fine)")
    parser.add_argument("--min_compared_bases", type=int, default=300, help="Minimum informative bases actually aligned against the best-matching reference to trust an identity score - this, not --min_region_covered, is the real quality gate for partial (e.g. amplicon-derived) coverage")
    parser.add_argument("--ambiguous_margin", type=float, default=1.0, help="Report 'ambiguous' if best and second-best identity are within this many percentage points")
    parser.add_argument("--unresolved_floor", type=float, default=90.0, help="Report 'unresolved' if the best identity is below this percentage")
    parser.add_argument("--output", required=True, help="Output TSV path (one row)")
    args = parser.parse_args()

    columns = ["sample", "gene", "genotype", "region_covered", "match_confidence", "notes"]

    def write_row(genotype, region_covered, match_confidence, notes):
        with open(args.output, "w") as fh:
            fh.write("\t".join(columns) + "\n")
            fh.write("\t".join([
                args.sample, args.gene, genotype,
                region_covered, match_confidence, notes,
            ]) + "\n")

    sample_seq = read_one_record(args.consensus, args.gene)
    if sample_seq is None:
        write_row("insufficient data", "0/0 (0.0%)", "", f"no '{args.gene}' record found in consensus FASTA")
        return

    if args.strand == "-":
        sample_seq = str(Seq(sample_seq).reverse_complement())

    gene_length = len(sample_seq)
    n_informative, frac_covered = region_coverage(sample_seq)
    region_covered_str = f"{n_informative}/{gene_length} ({frac_covered * 100:.1f}%)"

    if frac_covered < args.min_region_covered:
        write_row("insufficient data", region_covered_str, "", "below --min_region_covered")
        return

    panel = read_panel(args.reference_panel)
    if not panel:
        sys.exit(f"Reference panel {args.reference_panel} contains no sequences")

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    scores = {}
    n_compared_by_label = {}
    for label, ref_seq in panel.items():
        pct_identity, n_compared = compare_to_reference(sample_seq, ref_seq, aligner)
        scores[label] = pct_identity
        n_compared_by_label[label] = n_compared

    ranked = sorted(scores, key=scores.get, reverse=True)
    best, second = ranked[0], (ranked[1] if len(ranked) > 1 else None)
    best_score = scores[best]
    best_n_compared = n_compared_by_label[best]

    if best_n_compared < args.min_compared_bases:
        write_row("insufficient data", region_covered_str, f"{best_score:.1f}%",
                  f"only {best_n_compared} informative bases aligned to any reference (< --min_compared_bases)")
        return

    confidence = f"{best_score:.1f}%"
    if second is not None:
        confidence += f" (vs {scores[second]:.1f}% for {second})"

    notes = []

    if best_score < args.unresolved_floor:
        genotype = "unresolved"
        notes.append(f"best match {best} at only {best_score:.1f}% identity, below --unresolved_floor ({args.unresolved_floor}%)")
    elif second is not None and (best_score - scores[second]) < args.ambiguous_margin:
        genotype = f"ambiguous ({best}/{second})"
        notes.append(f"top two matches within --ambiguous_margin ({args.ambiguous_margin} pp)")
    else:
        genotype = best
        mixed_with = detect_mixed_signal(sample_seq, panel[best], aligner, best, panel, panel[best])
        if mixed_with:
            notes.append(f"possible mixed infection with {'/'.join(mixed_with)} signal present")

    write_row(genotype, region_covered_str, confidence, "; ".join(notes))


if __name__ == "__main__":
    main()
