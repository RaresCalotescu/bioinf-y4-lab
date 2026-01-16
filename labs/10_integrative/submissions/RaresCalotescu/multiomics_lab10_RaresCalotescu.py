from __future__ import annotations

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


# --------------------
# Config
# --------------------
HANDLE = "RaresCalotescu"

BASE_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}")
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SNP_FILE = DATA_DIR / f"snp_matrix_{HANDLE}.csv"
EXPR_FILE = DATA_DIR / f"expression_matrix_{HANDLE}.csv"
LABELS_FILE = DATA_DIR / f"labels_{HANDLE}.csv"

OUT_CONCAT = OUT_DIR / f"multiomics_concat_{HANDLE}.csv"
OUT_PCA_SNP = OUT_DIR / f"pca_snp_{HANDLE}.png"
OUT_PCA_EXPR = OUT_DIR / f"pca_expr_{HANDLE}.png"
OUT_PCA_JOINT = OUT_DIR / f"pca_joint_{HANDLE}.png"
OUT_PAIRS = OUT_DIR / f"snp_gene_pairs_{HANDLE}.csv"
OUT_CLUSTER = OUT_DIR / f"cluster_joint_{HANDLE}.csv"  # bonus


# --------------------
# Helpers
# --------------------
def read_matrix_csv(path: Path) -> pd.DataFrame:
    """
    Reads a matrix CSV and returns a DataFrame.
    Accepts either:
      - features x samples (preferred): rows=features, columns=samples
      - samples x features: rows=samples, columns=features
    We'll auto-detect later using intersection of sample names.
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    df = pd.read_csv(path, index_col=0)
    # Clean accidental whitespace
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    return df


def detect_samples_axis(df_a: pd.DataFrame, df_b: pd.DataFrame) -> tuple[str, str]:
    """
    Decide whether samples are in columns or in index by checking overlaps between A and B.
    Returns (a_samples_axis, b_samples_axis) each in {"columns","index"}.
    """
    a_col = set(df_a.columns)
    a_idx = set(df_a.index)
    b_col = set(df_b.columns)
    b_idx = set(df_b.index)

    score_a_cols_b_cols = len(a_col & b_col)
    score_a_cols_b_idx = len(a_col & b_idx)
    score_a_idx_b_cols = len(a_idx & b_col)
    score_a_idx_b_idx = len(a_idx & b_idx)

    # pick best pairing by max overlap
    scores = {
        ("columns", "columns"): score_a_cols_b_cols,
        ("columns", "index"): score_a_cols_b_idx,
        ("index", "columns"): score_a_idx_b_cols,
        ("index", "index"): score_a_idx_b_idx,
    }
    best = max(scores, key=scores.get)

    # if overlap is tiny everywhere, assume both are features x samples (samples in columns)
    if scores[best] < 2:
        warnings.warn("[WARN] Could not confidently detect sample axis; assuming samples are in columns.")
        return ("columns", "columns")

    return best


def to_features_by_samples(df: pd.DataFrame, samples_axis: str) -> pd.DataFrame:
    """Return a DataFrame with features as rows and samples as columns."""
    if samples_axis == "columns":
        return df
    elif samples_axis == "index":
        return df.T
    else:
        raise ValueError("samples_axis must be 'columns' or 'index'")


def zscore_features(df_fxS: pd.DataFrame) -> pd.DataFrame:
    """
    Z-score each feature across samples.
    Input: features x samples
    Output: features x samples (standardized)
    """
    X = df_fxS.values.astype(float)
    scaler = StandardScaler(with_mean=True, with_std=True)
    Xz = scaler.fit_transform(X.T).T  # scale per row(feature) across samples
    return pd.DataFrame(Xz, index=df_fxS.index, columns=df_fxS.columns)


def load_labels(samples: list[str]) -> pd.Series:
    """
    Try to load labels/subtypes. Accepts a CSV with columns:
      - Sample + Subtype  OR
      - Sample + Label
    Returns a Series indexed by sample.
    If missing, generates dummy labels (only for visualization).
    """
    if LABELS_FILE.exists():
        lab = pd.read_csv(LABELS_FILE)
        lab.columns = [c.strip() for c in lab.columns]
        if "Sample" not in lab.columns:
            raise ValueError(f"{LABELS_FILE} must contain a 'Sample' column.")

        if "Subtype" in lab.columns:
            y = lab.set_index("Sample")["Subtype"]
        elif "Label" in lab.columns:
            y = lab.set_index("Sample")["Label"]
        else:
            raise ValueError(f"{LABELS_FILE} must contain 'Subtype' or 'Label' column.")

        y.index = y.index.astype(str).str.strip()
        y = y.reindex(samples)
        return y

    warnings.warn("[WARN] labels_<handle>.csv not found -> generating dummy labels (for plots only).")
    dummy = [f"class_{i%3}" for i in range(len(samples))]
    return pd.Series(dummy, index=samples, name="Subtype")


def plot_pca(X_samples_by_features: pd.DataFrame, y: pd.Series, title: str, out_png: Path) -> None:
    """
    PCA(2) scatter colored by y (subtype).
    """
    # keep only samples with labels
    common = [s for s in X_samples_by_features.index if s in y.index and pd.notna(y.loc[s])]
    X = X_samples_by_features.loc[common].values.astype(float)
    y2 = y.loc[common].astype(str)

    pca = PCA(n_components=2, random_state=42)
    Z = pca.fit_transform(X)

    plt.figure(figsize=(7, 5))
    for cls in sorted(y2.unique()):
        m = (y2 == cls).values
        plt.scatter(Z[m, 0], Z[m, 1], label=cls, edgecolor="k", linewidths=0.3)

    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def compute_snp_gene_correlations(
    snp_fxS_z: pd.DataFrame,
    expr_fxS_z: pd.DataFrame,
    threshold: float = 0.5,
    chunk_genes: int = 2000,
) -> pd.DataFrame:
    """
    Pearson correlations between each SNP (feature) and each gene (feature)
    using standardized matrices:
      r = (X_snp @ X_gene.T) / (n_samples - 1)
    Returns a long table with pairs where |r| > threshold.
    Uses chunking to avoid RAM issues.
    """
    samples = snp_fxS_z.columns
    assert list(samples) == list(expr_fxS_z.columns)

    Xs = snp_fxS_z.values.astype(float)   # snps x samples (standardized)
    Xe = expr_fxS_z.values.astype(float)  # genes x samples (standardized)
    n = len(samples)

    denom = max(n - 1, 1)
    snp_names = snp_fxS_z.index.to_numpy()
    gene_names = expr_fxS_z.index.to_numpy()

    rows = []
    for start in range(0, Xe.shape[0], chunk_genes):
        end = min(start + chunk_genes, Xe.shape[0])
        Xe_chunk = Xe[start:end, :]  # chunk_genes x samples

        # snps x chunk_genes
        R = (Xs @ Xe_chunk.T) / denom

        # find indices where abs(r) > threshold
        ii, jj = np.where(np.abs(R) > threshold)
        for a, b in zip(ii.tolist(), jj.tolist()):
            rows.append((snp_names[a], gene_names[start + b], float(R[a, b])))

    out = pd.DataFrame(rows, columns=["SNP", "Gene", "r"])
    out["abs_r"] = out["r"].abs()
    out = out.sort_values(["abs_r"], ascending=False).reset_index(drop=True)
    return out


# --------------------
# Main
# --------------------
def main() -> None:
    # Load
    snp_raw = read_matrix_csv(SNP_FILE)
    expr_raw = read_matrix_csv(EXPR_FILE)

    # Detect where samples are
    snp_axis, expr_axis = detect_samples_axis(snp_raw, expr_raw)
    snp_fxS = to_features_by_samples(snp_raw, snp_axis)
    expr_fxS = to_features_by_samples(expr_raw, expr_axis)

    # Common samples
    common_samples = sorted(set(snp_fxS.columns) & set(expr_fxS.columns))
    if len(common_samples) < 3:
        raise ValueError(f"Too few common samples: {len(common_samples)}. Check that both files share sample IDs.")

    snp_fxS = snp_fxS[common_samples]
    expr_fxS = expr_fxS[common_samples]

    print(f"[INFO] SNP matrix (features x samples): {snp_fxS.shape}")
    print(f"[INFO] Expr matrix (features x samples): {expr_fxS.shape}")
    print(f"[INFO] Common samples: {len(common_samples)}")

    # Z-score per feature
    snp_z = zscore_features(snp_fxS)
    expr_z = zscore_features(expr_fxS)

    # Joint concat (features stacked)
    joint_fxS = pd.concat([snp_z, expr_z], axis=0)  # (snp+genes) x samples

    # Save integrated matrix as samples x features (easier for ML)
    joint_SxF = joint_fxS.T
    joint_SxF.to_csv(OUT_CONCAT)
    print(f"[OK] Saved integrated matrix: {OUT_CONCAT}")

    # Labels
    y = load_labels(common_samples)

    # PCA plots (samples x features)
    plot_pca(snp_z.T, y, "PCA — SNPs only", OUT_PCA_SNP)
    print(f"[OK] Saved: {OUT_PCA_SNP}")

    plot_pca(expr_z.T, y, "PCA — Expression only", OUT_PCA_EXPR)
    print(f"[OK] Saved: {OUT_PCA_EXPR}")

    plot_pca(joint_SxF, y, "PCA — Joint (SNP + Expression)", OUT_PCA_JOINT)
    print(f"[OK] Saved: {OUT_PCA_JOINT}")

    # Cross-omics correlations
    pairs = compute_snp_gene_correlations(snp_z, expr_z, threshold=0.5, chunk_genes=2000)
    pairs.to_csv(OUT_PAIRS, index=False)
    print(f"[OK] Saved SNP–gene pairs |r|>0.5: {OUT_PAIRS} (n={len(pairs)})")

    # Bonus: clustering on joint
    try:
        X = joint_SxF.values.astype(float)
        kmeans = KMeans(n_clusters=3, random_state=42, n_init="auto")
        cl = kmeans.fit_predict(X)
        cl_df = pd.DataFrame({"Sample": joint_SxF.index, "Cluster": cl})
        cl_df.to_csv(OUT_CLUSTER, index=False)
        print(f"[OK] Saved clustering: {OUT_CLUSTER}")
    except Exception as e:
        warnings.warn(f"[WARN] Clustering skipped due to: {e}")

    print("[DONE] Lab 10 finished.")


if __name__ == "__main__":
    main()
