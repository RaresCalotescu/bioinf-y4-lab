from pathlib import Path
import pandas as pd
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator

def main():
    aln_path = Path("data/aligned.fasta")
    out_csv = Path("results/distances.csv")
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    aln = AlignIO.read(str(aln_path), "fasta")
    if len(aln) < 10:
        raise SystemExit(f"[EROARE] Ai {len(aln)} secvențe. Trebuie ≥10.")

    calc = DistanceCalculator("identity")  
    ids = dm.names
    mat = [[dm[i, j] for j in range(len(ids))] for i in range(len(ids))]
    df = pd.DataFrame(mat, index=ids, columns=ids)

    df.to_csv(out_csv, float_format="%.6f")
    print(f"[OK] Salvat: {out_csv}")
    print(df.head())

if __name__ == "__main__":
    main()
