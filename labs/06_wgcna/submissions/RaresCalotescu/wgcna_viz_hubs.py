import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

HANDLE = "RaresCalotescu"

BASE_DIR = f"labs/06_wgcna/submissions/{HANDLE}"
RESULTS_DIR = os.path.join(BASE_DIR, "results")

GEXF_PATH = os.path.join(RESULTS_DIR, f"network_tp53_{HANDLE}.gexf")
MODULES_PATH = os.path.join(RESULTS_DIR, f"modules_tp53_{HANDLE}.csv")

OUT_PNG = os.path.join(RESULTS_DIR, f"network_tp53_{HANDLE}.png")
OUT_HUBS = os.path.join(RESULTS_DIR, f"hubs_tp53_{HANDLE}.csv")


def safe_int(x, default=-1):
    """Converteste valori gen '0', '0.0', 0.0, NaN -> int."""
    try:
        if pd.isna(x):
            return default
        return int(float(x))
    except Exception:
        return default


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)

    if not os.path.exists(GEXF_PATH):
        raise FileNotFoundError(f"Nu gasesc GEXF: {GEXF_PATH}")
    if not os.path.exists(MODULES_PATH):
        raise FileNotFoundError(f"Nu gasesc modules CSV: {MODULES_PATH}")

    G = nx.read_gexf(GEXF_PATH)
  
    G = nx.relabel_nodes(G, {n: str(n) for n in G.nodes()})

  
    mod_df = pd.read_csv(MODULES_PATH)

 
    cols_lower = {c.lower(): c for c in mod_df.columns}
    gene_col = cols_lower.get("gene", None)
    module_col = cols_lower.get("module", None)

    if gene_col is None or module_col is None:
 
        gene_col = mod_df.columns[0]
        module_col = mod_df.columns[1]

    mod_map = {}
    for _, row in mod_df.iterrows():
        gene = str(row[gene_col])
        mod_map[gene] = safe_int(row[module_col], default=-1)

 
    missing = 0
    for n in G.nodes():
        if n in mod_map:
            G.nodes[n]["module"] = mod_map[n]
        else:
            G.nodes[n]["module"] = -1
            missing += 1

    if missing > 0:
        print(f"[WARN] {missing} noduri nu au fost gasite in modules.csv (module=-1)")

 
    deg = dict(G.degree())
  
    module_to_nodes = {}
    for n in G.nodes():
        m = safe_int(G.nodes[n].get("module", -1), default=-1)
        module_to_nodes.setdefault(m, []).append(n)

    hub_rows = []
    hub_set = set()

   
    for m, nodes in module_to_nodes.items():
        nodes_sorted = sorted(nodes, key=lambda x: deg.get(x, 0), reverse=True)
        top_k = max(1, int(round(0.10 * len(nodes_sorted))))
        top_nodes = nodes_sorted[:top_k]

        for rank, n in enumerate(top_nodes, start=1):
            hub_set.add(n)
            hub_rows.append(
                {
                    "gene": n,
                    "module": m,
                    "degree": deg.get(n, 0),
                    "rank_in_module": rank,
                    "module_size": len(nodes_sorted),
                }
            )

    hubs_df = pd.DataFrame(hub_rows).sort_values(["module", "degree"], ascending=[True, False])
    hubs_df.to_csv(OUT_HUBS, index=False)
    print(f"[OK] Saved hubs: {OUT_HUBS}")

  
    pos = nx.spring_layout(G, seed=42, k=None)

    modules = [safe_int(G.nodes[n].get("module", -1), default=-1) for n in G.nodes()]
    min_m = min(modules) if modules else -1
    max_m = max(modules) if modules else 0

 
    sizes = []
    for n in G.nodes():
        base = 150 + 60 * deg.get(n, 0)
        if n in hub_set:
            base *= 1.6
        sizes.append(base)

    plt.figure(figsize=(12, 9))
    ax = plt.gca()
    ax.set_title(f"TP53 co-expression network (modules + hubs) - {HANDLE}")
    ax.axis("off")

   
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.8)

   
    nodes_collection = nx.draw_networkx_nodes(
        G,
        pos,
        node_size=sizes,
        node_color=modules,
        cmap=plt.get_cmap("tab10"),
        vmin=min_m,
        vmax=max_m,
        alpha=0.95,
    )


    if len(hub_set) > 0:
        nx.draw_networkx_nodes(
            G,
            pos,
            nodelist=list(hub_set),
            node_size=[200 + 80 * deg.get(n, 0) for n in hub_set],
            node_color=[safe_int(G.nodes[n].get("module", -1), default=-1) for n in hub_set],
            cmap=plt.get_cmap("tab10"),
            vmin=min_m,
            vmax=max_m,
            edgecolors="red",
            linewidths=2.0,
            alpha=1.0,
        )

    
    cbar = plt.colorbar(nodes_collection, shrink=0.75)
    cbar.set_label("Module (Louvain)")

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=200)
    plt.close()
    print(f"[OK] Saved network plot: {OUT_PNG}")


if __name__ == "__main__":
    main()
