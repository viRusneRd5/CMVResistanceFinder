#!/usr/bin/env python3
"""
Generate the small amplicon-style FASTQ test dataset used by `-profile test`.

Simulates paired-end Illumina reads tiled across the six CMV resistance genes of
the Merlin reference (NC_006273.2), with a fixed set of resistance mutations
spiked in homozygously so the pipeline's variant calling and phenotype matching
can be asserted end to end.

Also tiles reads (unmodified - no spiked edits) across UL55 (gB) and UL73 (gN)
so --genotype_gb/--genotype_gn have something to classify in CI. No sequence
editing needed: Merlin's own UL55/UL73 sequences already classify confidently
against the curated reference panels (gB3 and gN4d respectively - verified with
bin/genotype_classify.py directly against the reference panels before relying
on it here), so plain wild-type coverage gives a real, deterministic expected
genotype for nf-test to assert against.

Deterministic (fixed seed). Re-run after changing the spiked mutations or the
reference:

    python assets/test_data/make_mock_reads.py

Outputs (next to this script):
    reads_1.fastq.gz, reads_2.fastq.gz   - the test reads
    expected_mutations.tsv               - the spiked protein changes and the
                                           genomic edits that produced them
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF_FASTA = HERE.parent / "NC_006273.2.fasta"
REF_GFF = HERE.parent / "NC_006273.2.gff"
CHROM = "NC_006273.2"

SEED = 20260910
FRAG_LEN = 340
READ_LEN = 150
FRAG_STEP = 4           # spacing between fragment start positions -> ~75x depth
ERROR_RATE = 0.002      # per-base substitution error
BASE_QUAL = 37

# gene -> list of CDS segments (start, end, strand) in the order they are
# transcribed (5' segment first). 1-based inclusive, matching the GFF.
GENES: dict[str, list[tuple[int, int, str]]] = {
    "UL27": [(33182, 35008, "-")],
    "UL51": [(74576, 75049, "-")],
    "UL54": [(78194, 81922, "-")],
    "UL56": [(84752, 87304, "-")],
    "UL89": [(138810, 139697, "-"), (133770, 134906, "-")],
    "UL97": [(141798, 143921, "+")],
}

# Genotyping-only regions: tiled for read coverage the same way as GENES above,
# but with no spiked edits and no CDS/strand structure needed (nothing is
# translated for these - genotype_classify.py handles strand/orientation
# itself). 1-based inclusive, matching assets/genotype_genes.bed.
GENOTYPE_REGIONS: dict[str, tuple[int, int]] = {
    "UL55": (82066, 84789),
    "UL73": (107051, 107458),
}

# The resistance mutations to spike, as (gene, 1-based protein position, wild-type
# AA, mutant AA). "del" as the mutant AA means an in-frame deletion of that codon.
SPIKE = [
    ("UL97", 594, "A", "V"),
    ("UL54", 809, "A", "V"),
    ("UL27", 218, "?", "del"),
    ("UL51", 91, "P", "S"),
    ("UL56", 25, "C", "F"),
    ("UL89", 329, "N", "S"),
]

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_reference() -> str:
    seq = []
    with REF_FASTA.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def codon_genome_positions(segments, codon_index):
    """0-based genome positions (on the + strand) of the three bases of
    `codon_index` (0-based) of a CDS made of `segments`, in CDS reading order."""
    strand = segments[0][2]
    # CDS bases as (genome_pos_0based) in transcription order
    cds_positions = []
    for start, end, _ in segments:
        rng = range(start - 1, end) if strand == "+" else range(end - 1, start - 2, -1)
        cds_positions.extend(rng)
    i = codon_index * 3
    return cds_positions[i:i + 3], strand


def base_in_cds_orientation(ref, gpos, strand):
    b = ref[gpos]
    return b if strand == "+" else b.translate(COMP)


def apply_edits(ref: str, edits: list[tuple[int, str, str]]) -> str:
    """edits: (0-based genome pos, ref_bases, alt_bases) on the + strand,
    applied right-to-left so positions stay valid."""
    out = ref
    for pos, rb, ab in sorted(edits, key=lambda e: -e[0]):
        assert out[pos:pos + len(rb)] == rb, (pos, out[pos:pos + len(rb)], rb)
        out = out[:pos] + ab + out[pos + len(rb):]
    return out


def build_edits(ref: str):
    edits = []
    records = []
    for gene, aa_pos, wt, mut in SPIKE:
        segs = GENES[gene]
        gpos3, strand = codon_genome_positions(segs, aa_pos - 1)
        codon = "".join(base_in_cds_orientation(ref, p, strand) for p in gpos3)
        obs_aa = CODON_TABLE[codon]

        if mut == "del":
            # in-frame deletion of the whole codon: anchor on the preceding + base
            gmin, gmax = min(gpos3), max(gpos3)
            anchor = gmin - 1
            rb = ref[anchor:gmax + 1]
            ab = ref[anchor]
            edits.append((anchor, rb, ab))
            records.append((gene, f"{obs_aa}{aa_pos}del", CHROM, anchor + 1, rb, ab, strand))
            continue

        # pick an alt codon that codes `mut` with the fewest base changes
        best = None
        for alt_codon, alt_aa in CODON_TABLE.items():
            if alt_aa != mut:
                continue
            ndiff = sum(a != b for a, b in zip(codon, alt_codon))
            if best is None or ndiff < best[1]:
                best = (alt_codon, ndiff)
        alt_codon = best[0]

        for k, (cb, ab) in enumerate(zip(codon, alt_codon)):
            if cb == ab:
                continue
            gpos = gpos3[k]
            if strand == "+":
                rb, alt = cb, ab
            else:
                rb, alt = cb.translate(COMP), ab.translate(COMP)
            assert ref[gpos] == rb, (gene, gpos, ref[gpos], rb)
            edits.append((gpos, rb, alt))
        records.append((gene, f"{obs_aa}{aa_pos}{mut}", CHROM,
                        min(gpos3) + 1, codon, alt_codon, strand))
    return edits, records


def simulate(mut_ref: str, rng: random.Random):
    r1, r2 = [], []
    n = 0
    regions = {
        gene: (min(s for s, _, _ in segs) - 1, max(e for _, e, _ in segs))
        for gene, segs in GENES.items()
    }
    regions.update({
        gene: (start - 1, end) for gene, (start, end) in GENOTYPE_REGIONS.items()
    })
    for gene, (g_start, g_end) in regions.items():
        # centre fragments across the whole gene (fragments spill ~half their
        # length into the flanks) so every gene base, including the ends, gets
        # full paired-read depth.
        for centre in range(g_start, g_end, FRAG_STEP):
            start = centre - FRAG_LEN // 2
            if start < 0 or start + FRAG_LEN > len(mut_ref):
                continue
            frag = mut_ref[start:start + FRAG_LEN]
            if "N" in frag:
                continue
            fwd = list(frag[:READ_LEN])
            rev = list(revcomp(frag[-READ_LEN:]))
            for read in (fwd, rev):
                for i in range(len(read)):
                    if rng.random() < ERROR_RATE:
                        read[i] = rng.choice([b for b in "ACGT" if b != read[i]])
            name = f"mock:{gene}:{n}"
            q = chr(33 + BASE_QUAL) * READ_LEN
            r1.append(f"@{name}/1\n{''.join(fwd)}\n+\n{q}\n")
            r2.append(f"@{name}/2\n{''.join(rev)}\n+\n{q}\n")
            n += 1
    return r1, r2, n


def main():
    rng = random.Random(SEED)
    ref = read_reference()
    edits, records = build_edits(ref)
    mut_ref = apply_edits(ref, edits)

    r1, r2, n = simulate(mut_ref, rng)

    with gzip.open(HERE / "reads_1.fastq.gz", "wt") as fh:
        fh.write("".join(r1))
    with gzip.open(HERE / "reads_2.fastq.gz", "wt") as fh:
        fh.write("".join(r2))

    with (HERE / "expected_mutations.tsv").open("w") as fh:
        fh.write("gene\tprotein_change\tchrom\tgenome_pos\tref_codon_or_bases\talt_codon_or_bases\tstrand\n")
        for rec in records:
            fh.write("\t".join(str(x) for x in rec) + "\n")

    print(f"wrote {n} read pairs -> reads_1.fastq.gz / reads_2.fastq.gz")
    for rec in records:
        print("  spiked", rec[0], rec[1], "at", rec[2], rec[3])


if __name__ == "__main__":
    main()
