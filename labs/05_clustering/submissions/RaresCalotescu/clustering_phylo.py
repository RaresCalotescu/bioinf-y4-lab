from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO, Phylo

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from sklearn.metrics import silhouette_score

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster


BASE = Path("labs/05_clustering/submissions/RaresCalotescu")
FASTA_PATH = BASE / "data" / "input.fasta"
TREE_PATH = BASE / "data" / "tree_RaresCalotescu.nwk"
OUT_DIR = BASE / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def fasta_to_features(fasta_path: Path) -> pd.DataFrame:
    """
    Feature extraction simplu si robust:
    - lungime secventa
    - GC content
    - frecventele A/C/G/T
    """
    recs = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(recs) < 10:
        raise SystemExit(f"[EROARE] Ai doar {len(recs)} secvente. Trebuie >=10.")

    rows = []
    for r in recs:
        seq = str(r.seq).upper().replace("-", "")
        L = len(seq)

        a = seq.count("A")
        c = seq.count("C")
        g = seq.count("G")
        t = seq.count("T")
        denom = max(L, 1)

        gc = (g + c) / denom
        rows.append([r.id, L, gc, a / denom, c / denom, g / denom, t / denom])

    return pd.DataFrame(
        rows,
        columns=["Sequence_ID", "Length", "GC", "Freq_A", "Freq_C", "Freq_G", "Freq_T"],
    )


def plot_pca(X_scaled: np.ndarray, labels: np.ndarray, title: str, out_png: Path):
    pca = PCA(n_components=2, random_state=42)
    X2 = pca.fit_transform(X_scaled)

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(X2[:, 0], X2[:, 1], c=labels, alpha=0.85)
    plt.title(title)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.colorbar(sc, label="Cluster")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()


def save_dendrogram(Z, labels, out_png: Path):
    plt.figure(figsize=(12, 6))
    dendrogram(Z, labels=labels, leaf_rotation=90)
    plt.title("Hierarchical clustering (average linkage) - dendrogram")
    plt.xlabel("Sequences")
    plt.ylabel("Distance")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()


def annotate_tree_with_clusters(tree_path: Path, cluster_map: dict, out_png: Path):
    tree = Phylo.read(str(tree_path), "newick")

    for leaf in tree.get_terminals():
        if leaf.name in cluster_map:
            leaf.name = f"[C{cluster_map[leaf.name]}] {leaf.name}"

    fig = plt.figure(figsize=(10, 10), dpi=160)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)
    plt.tight_layout()
    fig.savefig(out_png)
    plt.close()


def safe_silhouette(X_scaled: np.ndarray, labels: np.ndarray):
    # silhouette cere >=2 clustere si fara un singur punct in cluster
    uniq = set(labels)
    if len(uniq) < 2:
        return None
    # daca un cluster are 1 punct, silhouette poate esua uneori
    try:
        return float(silhouette_score(X_scaled, labels))
    except Exception:
        return None


def main():
    # 1) Features
    df = fasta_to_features(FASTA_PATH)
    df.to_csv(OUT_DIR / "features.csv", index=False)

    X = df.drop(columns=["Sequence_ID"])
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # 2) Hierarchical (average linkage)
    Z = linkage(X_scaled, method="average")
    save_dendrogram(Z, df["Sequence_ID"].tolist(), OUT_DIR / "hierarchical_dendrogram.png")

    # alegem un k pentru hierarchical (3 e un default bun pentru discutie)
    hier_k = 3
    hier_labels = fcluster(Z, t=hier_k, criterion="maxclust") - 1
    sil_hier = safe_silhouette(X_scaled, hier_labels)

    # 3) KMeans - incercam mai multe K si alegem cel mai bun silhouette
    best_k, best_sil, best_km_labels = None, -1.0, None
    for k in [2, 3, 4, 5, 6]:
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        labs = km.fit_predict(X_scaled)
        sil = safe_silhouette(X_scaled, labs)
        if sil is not None and sil > best_sil:
            best_k, best_sil, best_km_labels = k, sil, labs

    plot_pca(X_scaled, best_km_labels, f"KMeans (best K={best_k}) - PCA", OUT_DIR / "kmeans_pca.png")

    # 4) DBSCAN (bonus) - testam cateva setari
    db_cfgs = [(0.8, 2), (1.0, 2), (1.2, 2), (1.5, 3)]
    best_db, best_db_sil, best_db_labels = None, None, None

    for eps, ms in db_cfgs:
        db = DBSCAN(eps=eps, min_samples=ms)
        labs = db.fit_predict(X_scaled)

        # silhouette doar pe punctele non-noise
        mask = labs != -1
        if mask.sum() > 2 and len(set(labs[mask])) > 1:
            sil = safe_silhouette(X_scaled[mask], labs[mask])
        else:
            sil = None

        if sil is not None and (best_db_sil is None or sil > best_db_sil):
            best_db, best_db_sil, best_db_labels = (eps, ms), sil, labs

    if best_db_labels is None:
        # fallback
        best_db = db_cfgs[0]
        best_db_labels = DBSCAN(eps=best_db[0], min_samples=best_db[1]).fit_predict(X_scaled)

    plot_pca(X_scaled, best_db_labels, f"DBSCAN (eps={best_db[0]}, min_samples={best_db[1]}) - PCA", OUT_DIR / "dbscan_pca.png")

    # 5) Salveaza cluster assignments
    out = pd.DataFrame({
        "Sequence_ID": df["Sequence_ID"],
        "Hierarchical_Cluster": hier_labels,
        "KMeans_Cluster": best_km_labels,
        "DBSCAN_Cluster": best_db_labels,
    })
    out.to_csv(OUT_DIR / "clusters.csv", index=False)

    # 6) Scoruri pentru raport
    scores = pd.DataFrame([{
        "silhouette_hierarchical_k3": sil_hier,
        "kmeans_best_k": best_k,
        "silhouette_kmeans_best": best_sil,
        "dbscan_best_eps": best_db[0],
        "dbscan_best_min_samples": best_db[1],
        "silhouette_dbscan_no_noise": best_db_sil,
    }])
    scores.to_csv(OUT_DIR / "scores.csv", index=False)

    # 7) Integrare cu arborele: adnotam frunzele cu clusterul (KMeans)
    cluster_map = dict(zip(out["Sequence_ID"], out["KMeans_Cluster"]))
    annotate_tree_with_clusters(TREE_PATH, cluster_map, OUT_DIR / "tree_annotated_kmeans.png")

    print("[OK] Generated outputs in:", OUT_DIR)
    print(" - features.csv, clusters.csv, scores.csv")
    print(" - hierarchical_dendrogram.png")
    print(" - kmeans_pca.png, dbscan_pca.png")
    print(" - tree_annotated_kmeans.png")


if __name__ == "__main__":
    main()
