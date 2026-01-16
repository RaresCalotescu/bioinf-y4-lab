import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# --------------------
# Config
# --------------------
HANDLE = "RaresCalotescu"
DATA_PATH = Path(f"labs/08_ML_flower/submissions/{HANDLE}/data/expression_matrix_{HANDLE}.csv")
OUT_DIR = Path(f"labs/08_ML_flower/submissions/{HANDLE}/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------
# Load expression matrix
# --------------------
df = pd.read_csv(DATA_PATH, index_col=0)
print("Original shape (genes x samples):", df.shape)

# --------------------
# Transpose → samples x genes
# --------------------
X = df.T
print("Transposed shape (samples x genes):", X.shape)

# --------------------
# Simulate labels (LAB PURPOSE)
# --------------------
# 3 clase artificiale
y = [i % 3 for i in range(len(X))]

le = LabelEncoder()
y_enc = le.fit_transform(y)

print("Classes:", le.classes_)

# --------------------
# Train / test split
# --------------------
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y_enc,
    test_size=0.3,
    random_state=42
)


# --------------------
# Random Forest
# --------------------
rf = RandomForestClassifier(
    n_estimators=100,
    random_state=42
)
rf.fit(X_train, y_train)

y_pred = rf.predict(X_test)

# --------------------
# Evaluation
# --------------------
print("\nClassification report:")
print(classification_report(y_test, y_pred))

cm = confusion_matrix(y_test, y_pred)

plt.figure(figsize=(5, 4))
sns.heatmap(
    cm,
    annot=True,
    fmt="d",
    cmap="Blues"
)
plt.xlabel("Predicted")
plt.ylabel("True")
plt.title("Confusion Matrix — Lab 8")
plt.tight_layout()
plt.savefig(OUT_DIR / "confusion_matrix_lab8.png", dpi=200)
plt.close()

print("✔ Saved confusion matrix to results/")
