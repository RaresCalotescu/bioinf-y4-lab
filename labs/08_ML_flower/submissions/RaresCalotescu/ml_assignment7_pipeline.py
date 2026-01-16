from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import networkx as nx  # not required, but harmless if installed (ignore if not)
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


# =========================
# Config
# =========================
HANDLE = "RaresCalotescu"

BASE_DIR = Path(f"labs/08_ml/submissions/{HANDLE}")
DATA_DIR = BASE_DIR / "data"
OUT_DIR = BASE_DIR / "results"

EXPR_CSV = DATA_DIR / f"expression_matrix_{HANDLE}.csv"
LABELS_CSV = DATA_DIR / f"labels_{HANDLE}.csv"  # optional: columns Sample,Label

OUT_DIR.mkdir(parents=True, exist_ok=True)

# RF params
RF_N_ESTIMATORS = 300
RF_RANDOM_STATE = 42

# Split
TEST_SIZE = 0.3
USE_STRATIFY_IF_POSSIBLE = True

# PCA/KMeans
PCA_N_COMPONENTS = 2
KMEANS_K = 3

# Semi-supervised
UNKNOWN_FRACTION = 0.4  # 30–50% e ok
MIN_LABELED_PER_CLASS = 2  # if less, semi-supervised becomes meaningless

# If too few samples, we can augment
AUGMENT_IF_FEW_SAMPLES = True
MIN_SAMPLES_TARGET = 30       # try to reach ~30 samples total
AUGMENT_NOISE_STD = 0.02      # small noise
AUGMENT_RANDOM_STATE = 42


# =========================
# Helpers
# =========================
def die(msg: str) -> None:
    raise SystemExit(msg)


def read_expression_matrix(expr_path: Path) -> pd.DataFrame:
    if not expr_path.exists():
        die(f"[ERROR] Nu găsesc expresia: {expr_path}")
    df = pd.read_csv(expr_path, index_col=0)
    if df.empty:
        die("[ERROR] expression_matrix este gol.")
    return df


def load_labels_if_any(X_samples: pd.DataFrame) -> pd.Series | None:
    """
    Returns y as a pandas Series indexed by sample name, or None if not found.
    Accepts:
      - labels_<handle>.csv with columns Sample,Label (recommended)
      - or if X already has a column Label (rare; not your case usually)
    """
    # 1) labels_<handle>.csv
    if LABELS_CSV.exists():
        lab = pd.read_csv(LABELS_CSV)
        if not {"Sample", "Label"}.issubset(lab.columns):
            warnings.warn(
                f"[WARN] {LABELS_CSV} există dar nu are coloanele Sample,Label. Ignor."
            )
        else:
            lab["Sample"] = lab["Sample"].astype(str)
            lab["Label"] = lab["Label"].astype(str)
            lab = lab.set_index("Sample")["Label"]
            # align to X samples
            common = X_samples.index.intersection(lab.index)
            if len(common) == 0:
                warnings.warn("[WARN] labels csv nu se potrivește cu sample-urile din X.")
            else:
                return lab.loc[common]

    # 2) Column Label in X (not typical)
    if "Label" in X_samples.columns:
        y = X_samples["Label"].astype(str)
        return y

    return None


def make_fake_labels(X_samples: pd.DataFrame, n_classes: int = 3) -> pd.Series:
    """
    Generate artificial labels only to make pipeline runnable.
    """
    warnings.warn(
        "[WARN] NU am găsit label-uri reale (nici coloană Label, nici labels_<handle>.csv). "
        "Generez label-uri artificiale doar ca să ruleze pipeline-ul. "
        "Pentru NOTE pe assignment, trebuie label-uri reale!"
    )
    y = pd.Series(
        [f"class_{i % n_classes}" for i in range(len(X_samples))],
        index=X_samples.index,
        name="Label",
    )
    return y


def can_stratify(y_enc: np.ndarray) -> bool:
    vals, counts = np.unique(y_enc, return_counts=True)
    # stratify needs at least 2 samples per class (and enough to split)
    return np.all(counts >= 2) and len(vals) >= 2


def augment_samples(X: pd.DataFrame, y: pd.Series, target_n: int, noise_std: float, random_state: int):
    """
    Duplicate samples with small gaussian noise, keeping labels.
    Works well for toy / teaching purposes.
    """
    if len(X) >= target_n:
        return X, y

    rng = np.random.default_rng(random_state)
    needed = target_n - len(X)

    X_list = [X]
    y_list = [y]

    # we will sample rows to duplicate
    base_idx = X.index.to_list()
    for i in range(needed):
        src = base_idx[i % len(base_idx)]
        row = X.loc[src].to_numpy(dtype=float)
        noisy = row + rng.normal(0.0, noise_std, size=row.shape)
        new_name = f"{src}__aug{i+1}"
        X_list.append(pd.DataFrame([noisy], index=[new_name], columns=X.columns))
        y_list.append(pd.Series([y.loc[src]], index=[new_name], name=y.name))

    X_aug = pd.concat(X_list, axis=0)
    y_aug = pd.concat(y_list, axis=0)

    return X_aug, y_aug


def save_confusion_png(cm: np.ndarray, classes: list[str], out_png: Path, title: str):
    plt.figure(figsize=(6, 5))
    plt.imshow(cm, interpolation="nearest")
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, rotation=45, ha="right")
    plt.yticks(tick_marks, classes)

    # annotate
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, str(cm[i, j]), ha="center", va="center")

    plt.ylabel("True")
    plt.xlabel("Predicted")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


# =========================
# Tasks
# =========================
def run_supervised_rf(X: pd.DataFrame, y: pd.Series, out_dir: Path):
    le = LabelEncoder()
    y_enc = le.fit_transform(y.astype(str))

    # split
    stratify = y_enc if (USE_STRATIFY_IF_POSSIBLE and can_stratify(y_enc)) else None
    if stratify is None and USE_STRATIFY_IF_POSSIBLE:
        warnings.warn("[WARN] Stratified split a eșuat (prea puține exemple într-o clasă). Folosesc split normal.")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y_enc,
        test_size=TEST_SIZE,
        random_state=RF_RANDOM_STATE,
        stratify=stratify
    )

    rf = RandomForestClassifier(
        n_estimators=RF_N_ESTIMATORS,
        random_state=RF_RANDOM_STATE,
        n_jobs=-1
    )
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)

    # classification report (safe labels)
    labels_present = sorted(set(y_test) | set(y_pred))
    target_names = [str(le.inverse_transform([c])[0]) for c in labels_present]

    report_txt = classification_report(
        y_test, y_pred,
        labels=labels_present,
        target_names=target_names,
        zero_division=0
    )

    (out_dir / f"classification_report_{HANDLE}.txt").write_text(report_txt)
    print(f"[OK] Saved classification report: {out_dir}/classification_report_{HANDLE}.txt")

    # confusion matrix png
    cm = confusion_matrix(y_test, y_pred, labels=labels_present)
    out_png = out_dir / f"confusion_rf_{HANDLE}.png"
    save_confusion_png(cm, target_names, out_png, title=f"Confusion Matrix — RF ({HANDLE})")
    print(f"[OK] Saved confusion matrix: {out_png}")

    # feature importances
    fi = pd.DataFrame({
        "gene": X.columns,
        "importance": rf.feature_importances_
    }).sort_values("importance", ascending=False)
    out_fi = out_dir / f"feature_importance_{HANDLE}.csv"
    fi.to_csv(out_fi, index=False)
    print(f"[OK] Saved feature importances: {out_fi}")

    return le, rf


def run_unsupervised_pca_kmeans(X: pd.DataFrame, y: pd.Series, out_dir: Path):
    # PCA
    pca = PCA(n_components=PCA_N_COMPONENTS, random_state=RF_RANDOM_STATE)
    X_pca = pca.fit_transform(X.to_numpy())

    # KMeans
    k = min(KMEANS_K, max(2, len(X)))  # keep sane
    km = KMeans(n_clusters=k, random_state=RF_RANDOM_STATE, n_init=10)
    clusters = km.fit_predict(X_pca)

    # Scatter colored by clusters, shaped by label (optional simple)
    plt.figure(figsize=(7, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=clusters)
    plt.title(f"PCA + KMeans (k={k}) — {HANDLE}")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    out_png = out_dir / f"sup_vs_unsup_scatter_{HANDLE}.png"
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved PCA+KMeans scatter: {out_png}")

    # Crosstab Label x Cluster
    ct = pd.crosstab(pd.Series(y.values, name="Label"), pd.Series(clusters, name="Cluster"))
    out_ct = out_dir / f"cluster_crosstab_{HANDLE}.csv"
    ct.to_csv(out_ct)
    print(f"[OK] Saved crosstab: {out_ct}")

    return clusters


def run_semi_supervised(X: pd.DataFrame, y: pd.Series, out_dir: Path):
    """
    Pseudo-labeling:
      - hide 30-50% labels
      - train RF on labeled only
      - predict unknown
      - retrain on full pseudo-labeled set
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y.astype(str))

    # hide a fraction
    rng = np.random.default_rng(RF_RANDOM_STATE)
    n = len(y_enc)
    idx = np.arange(n)
    rng.shuffle(idx)

    unknown_n = int(np.floor(n * UNKNOWN_FRACTION))
    unknown_idx = idx[:unknown_n]
    labeled_idx = idx[unknown_n:]

    # check labeled per class
    labeled_counts = np.bincount(y_enc[labeled_idx], minlength=len(le.classes_))
    if np.any(labeled_counts < MIN_LABELED_PER_CLASS):
        warnings.warn("[WARN] Prea puține probe în train pentru semi-supervised. Sar peste Task 5.")
        return

    X_labeled = X.iloc[labeled_idx]
    y_labeled = y_enc[labeled_idx]
    X_unknown = X.iloc[unknown_idx]

    rf1 = RandomForestClassifier(n_estimators=RF_N_ESTIMATORS, random_state=RF_RANDOM_STATE, n_jobs=-1)
    rf1.fit(X_labeled, y_labeled)

    pseudo = rf1.predict(X_unknown)

    # retrain on full pseudo-labeled
    y_full = y_enc.copy()
    y_full[unknown_idx] = pseudo

    rf2 = RandomForestClassifier(n_estimators=RF_N_ESTIMATORS, random_state=RF_RANDOM_STATE, n_jobs=-1)
    rf2.fit(X, y_full)

    # quick self-eval on labeled portion (not perfect, but gives a number)
    pred_labeled = rf2.predict(X_labeled)
    labels_present = sorted(set(y_labeled) | set(pred_labeled))
    target_names = [str(le.inverse_transform([c])[0]) for c in labels_present]

    report_txt = classification_report(
        y_labeled, pred_labeled,
        labels=labels_present,
        target_names=target_names,
        zero_division=0
    )
    out_rep = out_dir / f"semi_report_labeled_only_{HANDLE}.txt"
    out_rep.write_text(report_txt)
    print(f"[OK] Saved semi-supervised report (labeled-only eval): {out_rep}")


# =========================
# Main
# =========================
def main():
    expr = read_expression_matrix(EXPR_CSV)          # genes x samples
    X = expr.T.copy()                               # samples x genes

    print(f"[INFO] X shape (samples x genes): {X.shape}")

    # labels
    y = load_labels_if_any(X)
    if y is None:
        y = make_fake_labels(X, n_classes=3)

    # align X,y
    X = X.loc[y.index]
    print("[INFO] y distribution:")
    print(y.value_counts())

    # augment if too few samples
    if AUGMENT_IF_FEW_SAMPLES and len(X) < MIN_SAMPLES_TARGET:
        X, y = augment_samples(X, y, target_n=MIN_SAMPLES_TARGET, noise_std=AUGMENT_NOISE_STD, random_state=AUGMENT_RANDOM_STATE)
        print(f"[WARN] Augmented dataset to {len(X)} samples (toy augmentation).")
        print("[INFO] y distribution after augment:")
        print(y.value_counts())

    # Task 2
    run_supervised_rf(X, y, OUT_DIR)

    # Task 4
    run_unsupervised_pca_kmeans(X, y, OUT_DIR)

    # Task 5
    run_semi_supervised(X, y, OUT_DIR)

    print("[DONE] Pipeline finished. Check results/ folder.")


if __name__ == "__main__":
    main()
