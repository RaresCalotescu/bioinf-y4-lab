import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


HANDLE = "RaresCalotescu"

BASE = Path("labs/06_wgcna/submissions") / HANDLE
RESULTS = BASE / "results"

OUT_DIR = Path("labs/07_networkviz/submissions") / HANDLE
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    # 1. Load network
    G = nx.read_gexf(RESULTS / f"network_tp53_{HANDLE}.gexf")

    # 2. Load module assignments
    modules = pd.read_csv(RESULTS / f"modules_tp53_{HANDLE}.csv")
    gene2module = dict(zip(modules["gene"], modules["module"]))

    # 3. Assign module as node attribute
    for node in G.nodes():
        G.nodes[node]["module"] = gene2module.get(node, -1)

    # 4. Compute hub score (degree)
    degrees = dict(G.degree())
    nx.set_node_attributes(G, degrees, "degree")

    # 5. Layout
    pos = nx.spring_layout(G, seed=42)

    # 6. Node colors by module
    modules_list = [G.nodes[n]["module"] for n in G.nodes()]
    sizes = [300 + 200 * G.nodes[n]["degree"] for n in G.nodes()]

    plt.figure(figsize=(10, 10))
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=modules_list,
        cmap="tab10",
        node_size=sizes,
        alpha=0.85,
    )
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    nx.draw_networkx_labels(G, pos, font_size=9)

    plt.title("TP53 Co-expression Network (Modules + Hub Genes)")
    plt.axis("off")

    out_path = OUT_DIR / f"network_{HANDLE}.png"
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Saved network visualization to {out_path}")


if __name__ == "__main__":
    main()
