#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Exercițiu: Aliniere globală (Needleman–Wunsch)
"""

from pathlib import Path
import argparse
from Bio import SeqIO


def init_score_matrix_global(m: int, n: int, gap: int):
    score = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        score[i][0] = score[i - 1][0] + gap
    for j in range(1, n + 1):
        score[0][j] = score[0][j - 1] + gap
    return score


def score_cell_global(score, i: int, j: int, a: str, b: str, match: int, mismatch: int, gap: int):
    diag = score[i - 1][j - 1] + (match if a == b else mismatch)
    up = score[i - 1][j] + gap
    left = score[i][j - 1] + gap
    return max(diag, up, left)


def needleman_wunsch(seq1: str, seq2: str, match=1, mismatch=-1, gap=-2):
    m, n = len(seq1), len(seq2)
    score = init_score_matrix_global(m, n, gap)

    for i in range(1, m + 1):
        ai = seq1[i - 1]
        for j in range(1, n + 1):
            bj = seq2[j - 1]
            score[i][j] = score_cell_global(score, i, j, ai, bj, match, mismatch, gap)

    align1, align2 = "", ""
    i, j = m, n
    while i > 0 and j > 0:
        current = score[i][j]
        diag = score[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
        up = score[i - 1][j] + gap
        left = score[i][j - 1] + gap

        if current == diag:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif current == up:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    while i > 0:
        align1 = seq1[i - 1] + align1
        align2 = "-" + align2
        i -= 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j - 1] + align2
        j -= 1

    return align1, align2, score[m][n]


def load_two_sequences(fasta_path: Path, i1: int, i2: int):
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 2:
        raise SystemExit("[eroare] Fișierul trebuie să conțină cel puțin 2 secvențe.")
    if not (0 <= i1 < len(recs) and 0 <= i2 < len(recs)):
        raise SystemExit(f"[eroare] Indici invalizi (0..{len(recs)-1}).")
    return str(recs[i1].seq), str(recs[i2].seq), recs[i1].id, recs[i2].id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--i1", type=int, default=0)
    ap.add_argument("--i2", type=int, default=1)
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise SystemExit(f"[eroare] Nu găsesc fișierul: {fasta_path}")

    s1, s2, id1, id2 = load_two_sequences(fasta_path, args.i1, args.i2)
    a1, a2, sc = needleman_wunsch(s1, s2)

    print("=== Aliniere globală (Needleman–Wunsch) ===")
    print(f"{id1} vs {id2}")
    print(a1)
    print(a2)
    print("Score:", sc)


if __name__ == "__main__":
    main()
