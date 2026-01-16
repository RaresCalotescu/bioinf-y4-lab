from __future__ import annotations

from pathlib import Path
from typing import Tuple, Dict, Set, List
import warnings

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# --------------------
# Config
# --------------------
HANDLE = "RaresCalotescu"

BASE_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}")
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

DRUG_GENE_CSV = DATA_DIR / f"drug_gene_{HANDLE}.csv"
DISEASE_GENES_TXT = DATA_DIR / f"disease_genes_{HANDLE}.txt"

OUT_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"
OUT_NET_PNG = OUT_DIR / f"network_drug_gene_{HANDLE}.png"


# --------------------
# Helpers
# --------------------
def _find_col(df: pd.DataFrame, candidates: List[str]) -> str:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    raise ValueError(f"Nu gasesc nicio coloana din: {candidates}. Coloane disponibile: {list(df.columns)}")


def load_drug_gene(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Lipseste input-ul: {path}")

    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("drug_gene CSV este gol.")

    drug_col = _find_col(df, ["Drug", "drug", "compound", "medicament"])
    gene_col = _find_col(df, ["Gene", "gene", "target", "TargetGene"])

    df = df[[drug_col, gene_col]].copy()
    df.columns = ["Drug", "Gene"]
    df["Drug"] = df["Drug"].astype(str)
    df["Gene"] = df["Gene"].astype(str)

    df = df.dropna().drop_duplicates()
    return df


def build_bipartite_graph(df: pd.DataFrame) -> Tuple[nx.Graph, Set[str], Set[str]]:
    drugs = set(df["Drug"].unique().tolist())
    genes = set(df["Gene"].unique().tolist())

    G = nx.Graph()
    for d in drugs:
        G.add_node(d, bipartite="drug", node_type="drug")
    for g in genes:
        G.add_node(g, bipartite="gene", node_type="gene")

    # edges drug-gene
    G.add_edges_from(df[["Drug", "Gene"]].itertuples(index=False, name=None))
    return G, drugs, genes


def task1_summary(G: nx.Graph, drugs: Set[str]) -> pd.DataFrame:
    rows = []
    for d in sorted(drugs):
        # in bipartite, degree(d) = nr gene tinta
        rows.append({"Drug": d, "n_target_genes": int(G.degree(d))})
    out = pd.DataFrame(rows).sort_values("n_target_genes", ascending=False).reset_index(drop=True)
    out.to_csv(OUT_SUMMARY, index=False)
    print(f"[OK] Saved: {OUT_SUMMARY}")
    return out


# --------------------
# Task 2 — drug similarity (Jaccard)
# --------------------
def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def task2_drug_similarity(G: nx.Graph, drugs: Set[str], min_sim: float = 0.0) -> pd.DataFrame:
    # drug -> set of gene neighbors
    drug2genes: Dict[str, Set[str]] = {}
    for d in drugs:
        drug2genes[d] = set(G.neighbors(d))

    dlist = sorted(drugs)
    rows = []
    for i in range(len(dlist)):
        for j in range(i + 1, len(dlist)):
            d1, d2 = dlist[i], dlist[j]
            sim = jaccard(drug2genes[d1], drug2genes[d2])
            if sim > min_sim:
                rows.append({"Drug1": d1, "Drug2": d2, "Jaccard": sim})

    sim_df = pd.DataFrame(rows).sort_values("Jaccard", ascending=False).reset_index(drop=True)
    sim_df.to_csv(OUT_SIMILARITY, index=False)
    print(f"[OK] Saved: {OUT_SIMILARITY}")
    return sim_df


# --------------------
# Task 3 — disease proximity (average shortest path to disease genes)
# --------------------
def load_disease_genes(path: Path) -> Set[str]:
    if not path.exists():
        raise FileNotFoundError(f"Lipseste input-ul: {path}")
    genes = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            genes.add(line)
    if not genes:
        raise ValueError("disease_genes txt e gol.")
    return genes


def task3_disease_proximity(G: nx.Graph, drugs: Set[str], disease_genes: Set[str]) -> pd.DataFrame:
    present_dg = sorted([g for g in disease_genes if g in G.nodes()])
    if not present_dg:
        warnings.warn("Niciuna din disease_genes nu exista in graf. Ranking-ul nu va fi relevant.")

    rows = []
    for d in sorted(drugs):
        dists = []
        for g in present_dg:
            try:
                dist = nx.shortest_path_length(G, source=d, target=g)
                dists.append(dist)
            except nx.NetworkXNoPath:
                continue

        if dists:
            avg_dist = sum(dists) / len(dists)
            rows.append({"Drug": d, "avg_distance_to_disease_genes": avg_dist, "n_reachable_disease_genes": len(dists)})
        else:
            rows.append({"Drug": d, "avg_distance_to_disease_genes": float("inf"), "n_reachable_disease_genes": 0})

    pr = pd.DataFrame(rows).sort_values(
        ["avg_distance_to_disease_genes", "n_reachable_disease_genes"],
        ascending=[True, False],
    ).reset_index(drop=True)

    pr.to_csv(OUT_PRIORITY, index=False)
    print(f"[OK] Saved: {OUT_PRIORITY}")
    return pr


# --------------------
# Task 4 — visualization
# --------------------
def task4_visualize(G: nx.Graph, drugs: Set[str], genes: Set[str], summary: pd.DataFrame) -> None:
    # node sizes proportional to nr gene target (for drugs)
    drug_size = dict(zip(summary["Drug"], summary["n_target_genes"]))
    sizes = []
    colors = []

    for n in G.nodes():
        if n in drugs:
            sizes.append(200 + 200 * drug_size.get(n, 0))
            colors.append("royalblue")
        else:
            sizes.append(80)
            colors.append("crimson")

    pos = nx.spring_layout(G, seed=42)

    plt.figure(figsize=(12, 9))
    nx.draw_networkx_edges(G, pos, alpha=0.15, width=0.8)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes, linewidths=0)

    # etichete doar pentru drugs (mai lizibil)
    drug_labels = {d: d for d in drugs}
    nx.draw_networkx_labels(G, pos, labels=drug_labels, font_size=8)

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_NET_PNG, dpi=300)
    plt.close()
    print(f"[OK] Saved: {OUT_NET_PNG}")


def main():
    df = load_drug_gene(DRUG_GENE_CSV)
    G, drugs, genes = build_bipartite_graph(df)

    print(f"[INFO] Graph nodes: {G.number_of_nodes()}  edges: {G.number_of_edges()}")
    print(f"[INFO] Drugs: {len(drugs)}  Genes: {len(genes)}")

    # Task 1
    summary = task1_summary(G, drugs)

    # Task 2
    task2_drug_similarity(G, drugs, min_sim=0.0)

    # Task 3 (doar daca exista fisierul)
    if DISEASE_GENES_TXT.exists():
        disease_genes = load_disease_genes(DISEASE_GENES_TXT)
        task3_disease_proximity(G, drugs, disease_genes)
    else:
        warnings.warn(f"Nu exista {DISEASE_GENES_TXT}. Sar peste Task 3.")

    # Task 4
    task4_visualize(G, drugs, genes, summary)

    print("[DONE] Assignment 8 outputs generated in results/.")


if __name__ == "__main__":
    main()
