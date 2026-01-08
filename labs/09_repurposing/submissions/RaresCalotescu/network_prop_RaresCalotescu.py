import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pathlib import Path

HANDLE = "RaresCalotescu"
DATA_CSV = Path(f"labs/09_repurposing/submissions/{HANDLE}/data/drug_disease_interactions.csv")
OUT_DIR = Path(f"labs/09_repurposing/submissions/{HANDLE}/results")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PNG = OUT_DIR / f"network_{HANDLE}.png"
OUT_TOP_CENTRALITY = OUT_DIR / f"top_centrality_{HANDLE}.csv"
OUT_PATH_TXT = OUT_DIR / f"shortest_path_{HANDLE}.txt"

def node_type(n: str) -> str:
    if n.startswith("Drug_"):
        return "Drug"
    if n.startswith("Disease_"):
        return "Disease"
    if n.startswith("Protein_"):
        return "Protein"
    if n.startswith("Gene_"):
        return "Gene"
    return "Other"

def main():
    # Step 1: Load + Build graph
    df = pd.read_csv(DATA_CSV)
    if df.isna().any().any():
        df = df.dropna()

    G = nx.Graph()  # nedirecționat pentru lab
    for _, row in df.iterrows():
        src = str(row["Source"])
        tgt = str(row["Target"])
        itype = str(row["Interaction Type"])
        G.add_edge(src, tgt, interaction=itype)

    print(f"[INFO] Nodes: {G.number_of_nodes()} | Edges: {G.number_of_edges()}")

    # Step 2: Visualize (color nodes by type)
    colors = []
    for n in G.nodes():
        t = node_type(n)
        if t == "Drug":
            colors.append("red")
        elif t == "Disease":
            colors.append("blue")
        elif t == "Protein":
            colors.append("orange")
        else:
            colors.append("green")  # Gene/Other

    pos = nx.spring_layout(G, seed=42)

    plt.figure(figsize=(11, 8))
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=1.0)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=650, linewidths=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8)
    plt.title("Drug–Target–Gene–Disease Network (Toy) — Lab 9")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=250)
    plt.close()
    print(f"[OK] Saved network figure: {OUT_PNG}")

    # Step 3: Centrality metrics
    degree_c = nx.degree_centrality(G)
    betw_c = nx.betweenness_centrality(G)

    top_degree = sorted(degree_c.items(), key=lambda x: x[1], reverse=True)[:5]
    top_betw = sorted(betw_c.items(), key=lambda x: x[1], reverse=True)[:5]

    top_df = pd.DataFrame({
        "TopDegree_Node": [x[0] for x in top_degree],
        "TopDegree_Value": [x[1] for x in top_degree],
        "TopBetweenness_Node": [x[0] for x in top_betw],
        "TopBetweenness_Value": [x[1] for x in top_betw],
    })
    top_df.to_csv(OUT_TOP_CENTRALITY, index=False)
    print("[INFO] Top 5 degree:", top_degree)
    print("[INFO] Top 5 betweenness:", top_betw)
    print(f"[OK] Saved centrality table: {OUT_TOP_CENTRALITY}")

    # Step 4: Shortest path (select drug + disease)
    drug = "Drug_A"
    disease = "Disease_LungCancer"

    with open(OUT_PATH_TXT, "w", encoding="utf-8") as f:
        f.write(f"Selected drug: {drug}\nSelected disease: {disease}\n\n")
        if drug in G.nodes() and disease in G.nodes():
            try:
                path = nx.shortest_path(G, source=drug, target=disease)
                f.write("Shortest path:\n" + " -> ".join(path) + "\n")
                f.write(f"Path length: {len(path)-1}\n")
                print(f"[OK] Shortest path {drug} -> {disease}: {path}")
            except nx.NetworkXNoPath:
                f.write("No path exists between selected nodes.\n")
                print("[WARN] No path exists between selected nodes.")
        else:
            f.write("Drug or disease node not found in graph.\n")
            print("[WARN] Drug or disease node not found in graph.")

        # Step 5: Very short interpretation (for lab notes)
        f.write("\nInterpretare (scurt):\n")
        f.write("- Nodurile cu degree mare sunt hub-uri (ex. proteine/genes conectate la multe entități).\n")
        f.write("- Nodurile cu betweenness mare sunt „bottleneck”-uri (leagă subrețele → pot fi ținte cheie).\n")
        f.write("- Un drug mai aproape de o boală (drum scurt) poate sugera repurposing (prin target/gene comun).\n")
        f.write("- Limitări: rețeaua depinde de calitatea datelor, tipurile de interacțiuni, și faptul că graful e simplificat.\n")

    print(f"[OK] Saved shortest path + interpretation: {OUT_PATH_TXT}")

if __name__ == "__main__":
    main()

