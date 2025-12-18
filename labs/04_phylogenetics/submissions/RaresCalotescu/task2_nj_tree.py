from pathlib import Path
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def main():
    aln_path = Path("data/aligned.fasta")
    out_nwk = Path("results/nj_tree.nwk")
    out_png = Path("results/nj_tree.png")
    out_nwk.parent.mkdir(parents=True, exist_ok=True)

    # 1) Citește alinierea
    aln = AlignIO.read(str(aln_path), "fasta")
    print(f"[INFO] Loaded alignment with {len(aln)} sequences, length={aln.get_alignment_length()}")

    # 2) Matrice distanțe (identity -> dist = 1 - identity)
    calc = DistanceCalculator("identity")
    dm = calc.get_distance(aln)

    # 3) Neighbor-Joining
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # 4) Salvează arborele (Newick)
    Phylo.write(tree, str(out_nwk), "newick")
    print(f"[OK] Saved Newick: {out_nwk}")

    # 5) Print ASCII (pentru raport/notes)
    print("\n=== NJ Tree (ASCII preview) ===")
    Phylo.draw_ascii(tree)

    # 6) Bonus: salvează imagine
    fig = plt.figure(figsize=(9, 9), dpi=160)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    plt.tight_layout()
    fig.savefig(out_png)
    print(f"[OK] Saved plot: {out_png}")


if __name__ == "__main__":
    main()
