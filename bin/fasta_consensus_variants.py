#!/usr/bin/env python3
"""
Turn a per-sample multi-record FASTA (Sanger / assembled per-gene consensus
sequences) into a VCF of variants relative to the CMV reference, plus a table of
which region of each gene was assessed.

One record per gene: the header must start with the gene name (`>UL54`,
`>UL89_A`, or with a suffix such as `>UL54_patientX_2026`). Accepted gene names:
UL27, UL51, UL54, UL56, UL89_A, UL89_B, UL97.

Each record is locally aligned (nucleotide, both orientations, better score kept)
to that gene's region of the reference. The alignment is walked to emit
substitutions and in-frame indels; frameshift indels are dropped unless
--include-frameshifts. IUPAC ambiguity codes that include a non-reference base
are emitted as one VCF record per non-reference allele, flagged `MIXED` (and
`MIXED_REF` when the reference base is also part of the code).
"""

import argparse
import sys

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

GENES = ["UL27", "UL51", "UL54", "UL56", "UL89_A", "UL89_B", "UL97"]

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def parse_gene(header):
    tok = header.split()[0].upper()
    for g in sorted(GENES, key=len, reverse=True):
        if tok == g or tok.startswith(g + "_"):
            return g
    return None


def read_genes_bed(path):
    """gene -> (chrom, start, end) from the resistance-genes BED (chrom, start, end, name).

    NOTE: this pipeline's genes BED carries 1-based, inclusive coordinates (it is
    used that way by the report's coverage plots), not standard 0-based BED.
    """
    regions = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            regions[name] = (chrom, int(start), int(end))
    return regions


def make_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -1
    a.open_gap_score = -3
    a.extend_gap_score = -1
    return a


def best_alignment(aligner, query, ref_region):
    """Align query and its reverse complement to ref_region; return (alignment, is_rc)."""
    fwd = aligner.align(ref_region, query)[0]
    rc = str(Seq(query).reverse_complement())
    rev = aligner.align(ref_region, rc)[0]
    if rev.score > fwd.score:
        return rev, rc, True
    return fwd, query, False


def emit_variants(gene, chrom, gene_start, ref_region, alignment):
    """Walk one local alignment; yield VCF rows and return the assessed span + identity."""
    target_aln = str(alignment[0])   # aligned reference (with '-')
    query_aln = str(alignment[1])    # aligned consensus (with '-')
    t_offset = alignment.aligned[0][0][0]  # start of the alignment within ref_region

    rows = []
    ref_pos = t_offset               # 0-based offset into ref_region, points at next ref base
    aligned_cols = 0
    matches = 0
    covered = []                     # ref offsets with an informative consensus base

    i = 0
    n = len(target_aln)
    while i < n:
        r = target_aln[i]
        q = query_aln[i]

        if r != "-" and q != "-":
            aligned_cols += 1
            informative = q != "N"
            if informative:
                covered.append(ref_pos)
            if q == r:
                matches += 1
            elif informative:
                if q in "ACGT":
                    rows.append(_row(chrom, gene_start + ref_pos, r, q, mixed=False, mixed_ref=False))
                else:
                    bases = IUPAC.get(q, "")
                    alts = [b for b in bases if b != r]
                    mixed_ref = r in bases
                    for alt in alts:
                        rows.append(_row(chrom, gene_start + ref_pos, r, alt,
                                         mixed=True, mixed_ref=mixed_ref))
            ref_pos += 1
            i += 1

        elif q == "-":
            # deletion in the consensus: bases present in ref, absent in query
            start = i
            while i < n and query_aln[i] == "-":
                i += 1
            del_len = i - start
            del_bases = target_aln[start:i].replace("-", "")
            anchor_off = ref_pos - 1
            if anchor_off >= 0 and del_bases:
                anchor = ref_region[anchor_off]
                rows.append(_row(chrom, gene_start + anchor_off,
                                 anchor + del_bases, anchor,
                                 mixed=False, mixed_ref=False,
                                 frameshift=(del_len % 3 != 0)))
            ref_pos += del_len

        else:
            # insertion in the consensus: bases present in query, absent in ref
            start = i
            while i < n and target_aln[i] == "-":
                i += 1
            ins_bases = query_aln[start:i].replace("-", "").replace("N", "")
            anchor_off = ref_pos - 1
            if anchor_off >= 0 and ins_bases:
                anchor = ref_region[anchor_off]
                rows.append(_row(chrom, gene_start + anchor_off,
                                 anchor, anchor + ins_bases,
                                 mixed=False, mixed_ref=False,
                                 frameshift=(len(ins_bases) % 3 != 0)))

    identity = matches / aligned_cols if aligned_cols else 0.0
    if covered:
        span = (gene_start + min(covered), gene_start + max(covered))
    else:
        span = (0, 0)
    return rows, span, identity, aligned_cols


def _row(chrom, pos, ref, alt, *, mixed, mixed_ref, frameshift=False):
    info = ["SANGER"]
    if mixed:
        info.append("MIXED")
    if mixed_ref:
        info.append("MIXED_REF")
    return {
        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
        "filter": "frameshift" if frameshift else "PASS",
        "info": ";".join(info),
        "frameshift": frameshift,
    }


def write_vcf(path, sample, contig, rows):
    rows = sorted({(r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in rows}.values(),
                  key=lambda r: (r["chrom"], r["pos"], r["ref"], r["alt"]))
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write('##INFO=<ID=SANGER,Number=0,Type=Flag,Description="Called from a consensus (FASTA) sequence, not from reads">\n')
        fh.write('##INFO=<ID=MIXED,Number=0,Type=Flag,Description="Allele observed as an IUPAC ambiguity code in the consensus (mixed population)">\n')
        fh.write('##INFO=<ID=MIXED_REF,Number=0,Type=Flag,Description="The reference allele is also part of the ambiguity code (heterozygous with reference)">\n')
        fh.write('##FILTER=<ID=frameshift,Description="Frameshift indel from a consensus alignment (low confidence)">\n')
        fh.write(f"##contig=<ID={contig}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\t{r['filter']}\t{r['info']}\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    p.add_argument("--fasta", required=True, help="per-sample multi-record FASTA")
    p.add_argument("--genes-bed", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--output-vcf", required=True)
    p.add_argument("--output-assessed", required=True)
    p.add_argument("--output-gene-fasta", required=True,
                   help="the input records, re-oriented to reference sense")
    p.add_argument("--min-identity", type=float, default=0.9,
                   help="warn below this alignment identity")
    p.add_argument("--min-identity-hard", type=float, default=0.7,
                   help="fail if a record aligns to its gene below this identity")
    p.add_argument("--include-frameshifts", action="store_true")
    args = p.parse_args()

    reference = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.reference, "fasta")}
    regions = read_genes_bed(args.genes_bed)
    aligner = make_aligner()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    if not records:
        sys.exit(f"ERROR: {args.fasta} contains no FASTA records")

    all_rows = []
    assessed = []
    reoriented = []
    contig = None

    for rec in records:
        gene = parse_gene(rec.description)
        if gene is None:
            sys.exit(
                f"ERROR: sample {args.sample}: cannot map FASTA header '{rec.description}' "
                f"to a resistance gene. Accepted names (as a prefix): {', '.join(GENES)}"
            )
        if gene not in regions:
            sys.exit(f"ERROR: gene {gene} is not in {args.genes_bed}")

        chrom, gene_start, gene_end = regions[gene]
        contig = chrom
        ref_region = reference[chrom][gene_start - 1:gene_end]
        query = str(rec.seq).upper().replace("U", "T")

        alignment, oriented_query, is_rc = best_alignment(aligner, query, ref_region)
        rows, span, identity, aligned_cols = emit_variants(gene, chrom, gene_start, ref_region, alignment)

        if identity < args.min_identity_hard:
            sys.exit(
                f"ERROR: sample {args.sample}: the '{gene}' record aligns to the {gene} "
                f"reference region at only {identity:.0%} identity in both orientations "
                f"- it does not look like {gene}."
            )
        if identity < args.min_identity:
            sys.stderr.write(
                f"WARNING: sample {args.sample}: '{gene}' record aligns at {identity:.0%} "
                f"identity ({aligned_cols} aligned bp){' (reverse-complemented)' if is_rc else ''}\n"
            )

        if not args.include_frameshifts:
            rows = [r for r in rows if not r["frameshift"]]
        all_rows.extend(rows)

        assessed.append({
            "sample": args.sample, "gene": gene, "chrom": chrom,
            "assessed_start": span[0], "assessed_end": span[1],
            "aligned_bp": aligned_cols, "pct_identity": round(identity * 100, 1),
            "reverse_complemented": "yes" if is_rc else "no",
            "n_variants": len(rows),
        })
        reoriented.append((gene, oriented_query))

    write_vcf(args.output_vcf, args.sample, contig or "NC_006273.2", all_rows)

    cols = ["sample", "gene", "chrom", "assessed_start", "assessed_end",
            "aligned_bp", "pct_identity", "reverse_complemented", "n_variants"]
    with open(args.output_assessed, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for a in assessed:
            fh.write("\t".join(str(a[c]) for c in cols) + "\n")

    with open(args.output_gene_fasta, "w") as fh:
        for gene, seq in reoriented:
            fh.write(f">{gene}\n{seq}\n")


if __name__ == "__main__":
    main()
