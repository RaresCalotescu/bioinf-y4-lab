import pandas as pd
import numpy as np
from pathlib import Path

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt

HANDLE = "RaresCalotescu"

BASE_DIR = Path("data/sample")
SNP_FILE = BASE_DIR / "snp_data.csv"              # Samples x SNPs
EXPR_FILE = BASE_DIR / "expression_data.csv"      # Genes x Samples (de obicei)
PROT_FILE = BASE_DIR / "proteomics_data.csv"      # Proteins x Samples
PHENO_FILE = BASE_DIR / "phenotypes.csv"          # Samples x Phenotype

OUT_DIR = Path(f"labs/10_integrative/submissions/{HANDLE}/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PCA_PNG = OUT_DIR / f"pca_clusters_{HANDLE}.png"
OUT_FEATURES_CSV = OUT_DIR / f"integrated_features_{HANDLE}.csv"
OUT_SUMMARY_TXT = OUT_DIR / f"summary_{HANDLE}.txt"

def load_csv(path: Path, index_col=0) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, index_col=index_col)

def main():
    # --------------------------
    # Step 1: Load Data
    # --------------------------
    snp_df = load_csv(SNP_FILE, index_col=0)      # samples x SNPs
    expr_df = load_csv(EXPR_FILE, index_col=0)    # genes x samples
    prot_df = load_csv(PROT_FILE, index_col=0)    # proteins x samples
    pheno_df = load_csv(PHENO_FILE, index_col=0)  # samples x phenotype

    # --------------------------
    # Step 2: Preprocessing
    # --------------------------
    # Log2(x+1) expression (RNA-seq counts)
    expr_df = np.log2(expr_df + 1)

    # Z-score proteomics across samples (pe fiecare proteină)
    prot_scaled = StandardScaler().fit_transform(prot_df.T).T
    prot_df = pd.DataFrame(prot_scaled, index=prot_df.index, columns=prot_df.columns)

    # --------------------------
    # Step 2.5: Align samples
    # --------------------------
    # expr_df: genes x samples => samples are columns
    common_samples = (
        snp_df.index
        .intersection(expr_df.columns)
        .intersection(prot_df.columns)
        .intersection(pheno_df.index)
    )

    print(f"[WARN] Few common samples detected: {len(common_samples)} (toy dataset)")


    snp_df = snp_df.loc[common_samples]
    expr_df = expr_df[common_samples]
    prot_df = prot_df[common_samples]
    y = pheno_df.loc[common_samples].iloc[:, 0]  # prima coloană = phenotype

    # Encode phenotype dacă e text
    if y.dtype == "object":
        le = LabelEncoder()
        y_enc = le.fit_transform(y.astype(str))
        classes = list(le.classes_)
    else:
        y_enc = y.values
        classes = sorted(pd.unique(y))

    # --------------------------
    # Step 3: Feature selection
    # --------------------------
    top_genes = expr_df.var(axis=1).sort_values(ascending=False).head(200).index
    expr_top = expr_df.loc[top_genes]  # genes x samples

    top_proteins = prot_df.var(axis=1).sort_values(ascending=False).head(100).index
    prot_top = prot_df.loc[top_proteins]  # proteins x samples

    # --------------------------
    # Step 4: Integration (concat features)
    # X: samples x features
    # --------------------------
    X = pd.concat(
        [expr_top.T, prot_top.T, snp_df],
        axis=1,
        join="inner"
    )

    # Salvăm o copie a feature-urilor integrate (util în lab)
    X.to_csv(OUT_FEATURES_CSV)

    # --------------------------
    # Step 5: PCA + KMeans
    # --------------------------
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=min(5, X.shape[0], X.shape[1]), random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    k = 2 if len(common_samples) >= 2 else 1  # simplu: 2 clustere
    kmeans = KMeans(n_clusters=k, random_state=42, n_init="auto")
    clusters = kmeans.fit_predict(X_pca)

    # --------------------------
    # Step 6: RandomForest (demo, pe PCA)
    # (atenție: e “training accuracy”, nu evaluare corectă fără split)
    # --------------------------
    rf = RandomForestClassifier(n_estimators=200, random_state=42)
    rf.fit(X_pca, y_enc)
    acc = (rf.predict(X_pca) == y_enc).mean()

    # --------------------------
    # Step 7: Visualization
    # --------------------------
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=clusters, edgecolor="k")
    plt.title("PCA on Integrated Data (SNP + Expression + Proteomics)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    plt.savefig(OUT_PCA_PNG, dpi=300)
    plt.close()

    # Summary text
    with open(OUT_SUMMARY_TXT, "w", encoding="utf-8") as f:
        f.write(f"Common samples: {len(common_samples)}\n")
        f.write(f"Integrated X shape (samples x features): {X.shape}\n")
        f.write(f"PCA explained variance (first 2 PCs): {pca.explained_variance_ratio_[:2].tolist()}\n")
        f.write(f"KMeans k={k}\n")
        f.write(f"RF training accuracy (on PCA): {acc:.3f}\n")
        f.write(f"Phenotype classes: {classes}\n")

    print(f"[OK] Saved: {OUT_PCA_PNG}")
    print(f"[OK] Saved: {OUT_FEATURES_CSV}")
    print(f"[OK] Saved: {OUT_SUMMARY_TXT}")

if __name__ == "__main__":
    main()
