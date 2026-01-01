import pandas as pd
import numpy as np

def main():
    # load data
    df = pd.read_csv(
        "labs/06_wgcna/submissions/RaresCalotescu/data/expression_data.csv",
        index_col=0
    )
    print(f"[INFO] Loaded expression matrix: {df.shape}")

    
    df_log = np.log2(df + 1)

   
    variances = df_log.var(axis=1)
    threshold = 0.0
    df_filt = df_log.loc[variances > threshold]

  
    if df_filt.shape[0] < 5:
        print("[WARN] Too few genes after filtering, keeping all genes")
        df_filt = df_log.copy()

    print(f"[INFO] After variance filtering: {df_filt.shape}")

  
    out = "labs/06_wgcna/submissions/RaresCalotescu/data/expression_preprocessed.csv"
    df_filt.to_csv(out)
    print(f"[OK] Saved preprocessed data to {out}")

if __name__ == "__main__":
    main()
