#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx

def try_louvain(G: nx.Graph):
    """
    Prefer NetworkX built-in Louvain if available, else fall back to python-louvain.
    Returns dict: node -> module_id (int).
    """
    
    try:
        from networkx.algorithms.community import louvain_communities
        communities = louvain_communities(G, seed=42)
        node2mod = {}
        for i, comm in enumerate(communities):
            for n in comm:
                node2mod[n] = i
        return node2mod
    except Exception:
        pass

   
    try:
        import community as community_louvain  
        part = community_louvain.best_partition(G, random_state=42)
        # part: node -> module_id already
        return part
    except Exception as e:
        raise RuntimeError(
            "Nu pot rula Louvain. Incearca sa instalezi python-louvain:\n"
            "  pip install python-louvain\n"
            f"Detalii: {e}"
        )

def main():
    handle = "RaresCalotescu"

    base_dir = f"labs/06_wgcna/submissions/{handle}"
    data_path = os.path.join(base_dir, "data", "expression_preprocessed.csv")
    results_dir = os.path.join(base_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    if not os.path.exists(data_path):
        print(f"[ERROR] Nu gasesc: {data_path}")
        sys.exit(1)

   
    df = pd.read_csv(data_path, index_col=0)

    print(f"[INFO] Loaded preprocessed matrix: {df.shape} (genes x samples)")

    # Basic sanity
    if df.shape[0] < 2:
        print("[ERROR] Prea putine gene pentru corelatii (trebuie >= 2).")
        sys.exit(1)

    if df.shape[1] < 3:
        print("[WARN] Ai <3 sample-uri. Corelatiile pot fi instabile, dar continuam (pentru tema).")

  
    corr = df.T.corr(method="pearson") 
    corr = corr.fillna(0.0)

    
   
    abs_vals = np.abs(corr.values)
    
    abs_vals = abs_vals[~np.eye(abs_vals.shape[0], dtype=bool)]
    auto_thr = np.quantile(abs_vals, 0.90) if abs_vals.size > 0 else 0.7
    thr = max(0.7, float(auto_thr))

    print(f"[INFO] Using correlation threshold |r| >= {thr:.3f}")

    genes = corr.index.tolist()
    G = nx.Graph()
    G.add_nodes_from(genes)

    # Add edges
    added = 0
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            r = float(corr.iat[i, j])
            if abs(r) >= thr:
                G.add_edge(genes[i], genes[j], weight=abs(r), corr=r)
                added += 1

    print(f"[INFO] Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges (added={added})")

    if G.number_of_edges() == 0:
        print("[WARN] Nu s-a creat niciun edge la pragul ales. Scad pragul la 0.5 si reincerc...")
        thr = 0.5
        G = nx.Graph()
        G.add_nodes_from(genes)
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                r = float(corr.iat[i, j])
                if abs(r) >= thr:
                    G.add_edge(genes[i], genes[j], weight=abs(r), corr=r)
        print(f"[INFO] Graph (retry): {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

 
    node2mod = try_louvain(G)

    modules_df = pd.DataFrame({
        "gene": list(node2mod.keys()),
        "module": [node2mod[g] for g in node2mod.keys()]
    }).sort_values(["module", "gene"])

    out_modules = os.path.join(results_dir, f"modules_tp53_{handle}.csv")
    modules_df.to_csv(out_modules, index=False)

  
    out_gexf = os.path.join(results_dir, f"network_tp53_{handle}.gexf")
    nx.write_gexf(G, out_gexf)

    print(f"[OK] Saved modules: {out_modules}")
    print(f"[OK] Saved network:  {out_gexf}")

if __name__ == "__main__":
    main()
