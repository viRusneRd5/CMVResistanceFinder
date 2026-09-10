#!/usr/bin/env python3

import argparse
import pandas as pd
import re

# Three-letter to one-letter amino acid conversion
AA = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C",
    "Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
    "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P",
    "Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
    "Ter":"*","Stop":"*"
}


def convert_aa(hgvs):
    """
    Convert HGVS protein notation into abbreviated form.
    """

    if pd.isna(hgvs):
        return ""

    hgvs = hgvs.replace("p.", "")

    # -------------------------
    # Simple substitution
    # Val367Ile -> V367I
    # -------------------------
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|\*)", hgvs)
    if m:
        ref, pos, alt = m.groups()
        ref = AA.get(ref, ref)
        alt = "*" if alt == "*" else AA.get(alt, alt)
        return f"{ref}{pos}{alt}"

    # -------------------------
    # Stop codon
    # Ser857*
    # -------------------------
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)\*", hgvs)
    if m:
        ref, pos = m.groups()
        return f"{AA[ref]}{pos}*"

    # -------------------------
    # Frameshift deletion
    # Trp362fs -> DEL362
    # -------------------------
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)fs", hgvs)
    if m:
        _, pos = m.groups()
        return f"DEL{pos}"

    # -------------------------
    # In-frame deletion (SnpEff HGVS.p)
    # Leu218del          -> DEL218
    # Leu218_Ala221del   -> DEL218-221
    # -------------------------
    m = re.fullmatch(r"[A-Z][a-z]{2}(\d+)(?:_[A-Z][a-z]{2}(\d+))?del", hgvs)
    if m:
        start, end = m.groups()
        return f"DEL{start}" if end is None else f"DEL{start}-{end}"

    # -------------------------
    # In-frame insertion / duplication (SnpEff HGVS.p), keeping the inserted
    # residue(s) so it matches the database's INS<pos>_<residues> notation:
    #   Ser884dup                            -> INS884_S
    #   Ala884_Ser885insSer                  -> INS884_S
    #   Ser1100_Gly1101insLeuThrAlaProGlyVal -> INS1100_LTAPGV
    # -------------------------
    m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)dup", hgvs)
    if m:
        aa3, pos = m.groups()
        return f"INS{pos}_{AA.get(aa3, aa3)}"

    m = re.fullmatch(r"[A-Z][a-z]{2}(\d+)_[A-Z][a-z]{2}\d+ins((?:[A-Z][a-z]{2})+)", hgvs)
    if m:
        pos, inserted = m.groups()
        residues = "".join(
            AA.get(inserted[i:i + 3], "?") for i in range(0, len(inserted), 3)
        )
        return f"INS{pos}_{residues}"

    # range duplication or an insertion we can't fully resolve - keep the position
    m = re.fullmatch(r"[A-Z][a-z]{2}(\d+)(?:_[A-Z][a-z]{2}\d+)?(?:ins|dup).*", hgvs)
    if m:
        return f"INS{m.group(1)}"

    # -------------------------
    # Two amino acid substitution
    # AspLeu981GlyVal -> D981G/L982V
    # -------------------------
    m = re.fullmatch(
        r"([A-Z][a-z]{2})([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})([A-Z][a-z]{2})",
        hgvs
    )

    if m:
        ref1, ref2, pos, alt1, alt2 = m.groups()
        pos = int(pos)

        return (
            f"{AA[ref1]}{pos}{AA[alt1]}"
            f"/"
            f"{AA[ref2]}{pos+1}{AA[alt2]}"
        )

    # If nothing matched, keep original
    return hgvs


def fix_indels(row):
    """
    Convert frameshifts into DEL/INS notation using CODON information.
    """

    aa = row["EFF[*].AA"]
    codon = row["EFF[*].CODON"]

    if pd.isna(aa):
        return ""

    if "fs" not in aa:
        return convert_aa(aa)

    # deletion
    if isinstance(codon, str) and "del" in codon:
        pos = re.search(r"(\d+)", aa).group(1)

        deleted = re.search(r"del([A-Z]+)", codon)

        if deleted:
            n = len(deleted.group(1))
            if n == 1:
                return f"DEL{pos}"
            else:
                end = int(pos) + n - 1
                return f"DEL{pos}-{end}"

        return f"DEL{pos}"

    # insertion
    if isinstance(codon, str) and "ins" in codon:
        pos = re.search(r"(\d+)", aa).group(1)

        inserted = re.search(r"ins([A-Z]+)", codon)

        if inserted:
            #bases = inserted.group(1)

            # If SnpSift has already translated amino acids,
            # replace with appropriate one-letter AA here.
            return f"INS{pos}"

        return f"INS{pos}"

    return convert_aa(aa)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # Remove synonymous variants
    df = df[df["EFF[*].EFFECT"] != "synonymous_variant"].copy()

    # Create simplified AA column
    df["AA_change"] = df.apply(fix_indels, axis=1)

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()