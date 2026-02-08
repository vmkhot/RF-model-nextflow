# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
df_data = pd.read_csv('RF_results/raw/cluster_0_data_df.tsv',sep='\t')
df_data

# %%

# Example dataframe (replace with your own)
# data = {
#     "index": ["152_G","133_S","187_P","228_V","150_P","623_A","512_T",
#               "571_A","455_E","161_E","293_V","127_N","661_I"],
#     "shap_value": [0]*13,
#     "gini_imp": [0.0067,0.00566,0.0056,0.00534,0.0052,0.00507,0.00488,
#                  2.26e-05,2.26e-05,2.26e-05,2.26e-05,2.26e-05,2.26e-05],
#     "cluster_num": ["cluster_0"]*7 + ["cluster_1"]*6
# }
df = pd.read_csv('RF_results/Feature_importance.tsv',sep='\t')

# Amino acid groups (ordered)
aa_groups = [
    ("Positive charged", ["R", "H", "K"]),
    ("Negative charged", ["D", "E"]),
    ("Polar uncharged", ["S", "T", "N", "Q"]),
    ("Hydrophobic", ["A", "V", "I", "L", "M"]),
    ("Aromatic", ["F", "Y", "W"]),
    ("Special cases", ["C", "U", "G", "P"])
]
aa_order = [aa for _, group in aa_groups for aa in group]

# Function to extract AA from index (e.g., "152_G" → "G")
def get_aa(idx):
    return idx.split("_")[1]

df["AA"] = df["index"].apply(get_aa)

# Radar plot function
def radar_plot_for_cluster(cluster_df, cluster_name, ks=[10,20,50]):
    # Angles for each amino acid position
    N = len(aa_order)
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angles += angles[:1]  # close circle

    plt.figure(figsize=(8,8))
    ax = plt.subplot(111, polar=True)

    for k in ks:
        topk = cluster_df.nlargest(k, "gini_imp") if len(cluster_df) > k else cluster_df
        counts = topk["AA"].value_counts().reindex(aa_order, fill_value=0)
        values = counts.values.tolist()
        values += values[:1]  # close circle

        ax.plot(angles, values, linewidth=2, alpha=0.5)
        ax.fill(angles, values, alpha=0.2)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(aa_order, fontsize=12)
    ax.set_title(f"Radar plot for {cluster_name}", size=14, pad=20)
    plt.show()

# Generate one radar plot per cluster
for cluster_name, cluster_df in df.groupby("cluster_num"):
    radar_plot_for_cluster(cluster_df, cluster_name)


# %%
# ---- Example inputs ----
# Top features dataframe (replace with your real one)
cluster_num = "cluster_1"

feat_df = pd.read_csv('RF_results/Feature_importance.tsv',sep='\t')

# Biome dataframe (your big one)
df_biome = pd.read_csv(f"RF_results/raw/{cluster_num}_data_df.tsv",sep='\t')

# ---- Define AA order ----
# --- Define amino acid order by chemical properties ---
aa_groups = [
    ("Positive charged", ["R", "H", "K"], "crimson"),
    ("Negative charged", ["D", "E"], "royalblue"),
    ("Polar uncharged", ["S", "T", "N", "Q"], "seagreen"),
    ("Hydrophobic", ["A", "V", "I", "L", "M"], "goldenrod"),
    ("Aromatic", ["F", "Y", "W"], "purple"),
    ("Special cases", ["C", "U", "G", "P"], "darkorange")
]
aa_order = [aa for _, group, _ in aa_groups for aa in group]

def get_aa(idx):
    """Extract AA letter from feature index like '152_G' -> 'G'"""
    return idx.split("_")[1]

def radar_biome_plot(cluster_df, biome_df, cluster_name, top_k=50):
    # Top-k features for this cluster
    top_feats = cluster_df.nlargest(top_k, "gini_imp")["index"].tolist()
    feat_cols = [f for f in top_feats if f in biome_df.columns]

    if not feat_cols:
        print(f"No matching features for {cluster_name}")
        return

    sub = biome_df[feat_cols + ["ecosystem_subtype"]]

    biome_counts = {}
    for biome in ["Lake", "Oceanic"]:
        mask = sub["ecosystem_subtype"] == biome
        counts = sub.loc[mask, feat_cols].sum()
        aa_counts = counts.rename(lambda x: get_aa(x)).groupby(level=0).sum()
        aa_counts = aa_counts.reindex(aa_order, fill_value=0)

        total = aa_counts.sum()
        if total > 0:
            aa_counts = aa_counts / total
        biome_counts[biome] = aa_counts

    # ---- Radar plot ----
    N = len(aa_order)
    dtheta = 2*np.pi / N
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    plt.figure(figsize=(9,9))
    ax = plt.subplot(111, polar=True)

    for biome, color in zip(["Lake", "Oceanic"], ["#2B8CBF", "#9dcc7e"]):
        if biome not in biome_counts:
            continue
        values = biome_counts[biome].tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=2, label=biome, color=color, alpha=1)
        ax.fill(angles, values, alpha=0.5, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(aa_order, fontsize=11)

    # ---- Group coloured bands OUTSIDE the radar ----
    outer_radius = ax.get_rmax() * 1.05   # just outside the circle
    bar_height = ax.get_rmax() * 0.03     # thickness of band
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    legend_patches = []
    for group_name, group_aas, group_color in aa_groups:
        start_idx = aa_order.index(group_aas[0])
        end_idx = aa_order.index(group_aas[-1])

        # Start halfway before first AA, end halfway after last AA
        start_angle = start_idx * dtheta - dtheta/2
        end_angle = (end_idx+1) * dtheta - dtheta/2
        width = end_angle - start_angle

        ax.bar(
            x=start_angle, height=bar_height, width=width,
            bottom=outer_radius * 1.05, color=group_color, alpha=0.6, align="edge"
        )
        legend_patches.append(Patch(color=group_color, alpha=0.6, label=group_name))

    # ---- Legend for groups ----
    # Legend for biomes
    biome_handles = [
        Line2D([0], [0], color="#2B8CBF", lw=3, label="Lake"),
        Line2D([0], [0], color="#9dcc7e", lw=3, label="Oceanic"),
    ]
    # Combine group and biome legends
    handles = legend_patches + biome_handles
    ax.legend(handles=handles, title="AA Groups & Biomes",
              loc="upper right", bbox_to_anchor=(1.2, 1.1))
    ax.set_title(f"Top {top_k} features - {cluster_name}", size=14, pad=20)
    plt.show()


radar_biome_plot(feat_df, df_biome, cluster_num, top_k=100)

# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Top features dataframe (replace with your real one)
cluster_num = "cluster_mixed_model"

feat_df = pd.read_csv('RF_results/Feature_importance.tsv',sep='\t')

# Biome dataframe (your big one)
df_biome = pd.read_csv(f"RF_results/raw/{cluster_num}_data_df.tsv",sep='\t')

# --- Define amino acid order by chemical properties ---
aa_groups = [
    ("Positive charged", ["R", "H", "K"], "crimson"),
    ("Negative charged", ["D", "E"], "royalblue"),
    ("Polar uncharged", ["S", "T", "N", "Q"], "seagreen"),
    ("Hydrophobic", ["A", "V", "I", "L", "M"], "goldenrod"),
    ("Aromatic", ["F", "Y", "W"], "purple"),
    ("Special cases", ["C", "U", "G", "P"], "darkorange")
]
aa_order = [aa for _, group, _ in aa_groups for aa in group]

def get_aa(idx):
    """Extract amino acid from feature index like '152_G' -> 'G'"""
    return idx.split("_")[1]

def compare_biome_aas(cluster_df, biome_df, cluster_name, top_k=50, biomes=("Lake","Oceanic")):
    # --- Select top features ---
    top_feats = cluster_df.nlargest(top_k, "gini_imp")["index"].tolist()
    feat_cols = [f for f in top_feats if f in biome_df.columns]

    if not feat_cols:
        print(f"No matching features for {cluster_name}")
        return

    sub = biome_df[feat_cols + ["ecosystem_subtype"]]

    # --- Count amino acid frequencies per biome ---
    biome_freqs = {}
    for biome in biomes:
        mask = sub["ecosystem_subtype"] == biome
        counts = sub.loc[mask, feat_cols].sum()
        aa_counts = counts.rename(lambda x: get_aa(x)).groupby(level=0).sum()
        aa_counts = aa_counts.reindex(aa_order, fill_value=0)

        total = aa_counts.sum()
        if total > 0:
            aa_counts = aa_counts / total  # normalize to relative frequency
        biome_freqs[biome] = aa_counts

    # --- Prepare plotting data ---
    df_plot = pd.DataFrame(biome_freqs)
    df_plot = df_plot.fillna(0)

    # --- Plot ---
    fig, axes = plt.subplots(2, 1, figsize=(12,10), sharex=True)

    # Grouped bar chart
    bar_width = 0.35
    x = np.arange(len(aa_order))

    for i, biome in enumerate(biomes):
        axes[0].bar(x + i*bar_width, df_plot[biome], width=bar_width, label=biome, alpha=0.7)

    axes[0].set_xticks(x + bar_width/2)
    axes[0].set_xticklabels(aa_order, fontsize=12)
    axes[0].set_ylabel("Relative frequency")
    axes[0].set_title(f"Top {top_k} features - {cluster_name} (AA frequencies)")
    axes[0].legend()

    # Difference bar chart (Lake - Oceanic)
    diff = df_plot[biomes[0]] - df_plot[biomes[1]]
    colors = ["blue" if v > 0 else "green" for v in diff]

    axes[1].bar(x, diff, color=colors, alpha=0.7)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(aa_order, fontsize=12)
    axes[1].set_ylabel(f"{biomes[0]} - {biomes[1]}")
    axes[1].set_title("Difference in normalized frequencies")

    plt.tight_layout()
    plt.show()

compare_biome_aas(feat_df, df_biome, cluster_num, top_k=100)


# %%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

# Top features dataframe (replace with your real one)
cluster_num = "cluster_1"

feat_df = pd.read_csv('RF_results/Feature_importance.tsv',sep='\t')

# Biome dataframe (your big one)
df_biome = pd.read_csv(f"RF_results/raw/{cluster_num}_data_df.tsv",sep='\t')

# --- Define AA groups and base colors ---
aa_groups = [
    ("Positive charged", ["R", "H", "K"], "crimson"),
    ("Negative charged", ["D", "E"], "royalblue"),
    ("Polar uncharged", ["S", "T", "N", "Q"], "seagreen"),
    ("Hydrophobic", ["A", "V", "I", "L", "M"], "goldenrod"),
    ("Aromatic", ["F", "Y", "W"], "purple"),
    ("Special cases", ["C", "U", "G", "P"], "darkorange")
]
aa_order = [aa for _, group, _ in aa_groups for aa in group]
base_colors = {aa: color for _, group, color in aa_groups for aa in group}

def lighten_color(color, factor=0.5):
    """Lighten color by blending with white"""
    c = np.array(mcolors.to_rgb(color))
    return tuple(c + (1 - c) * factor)

def darken_color(color, factor=0.5):
    """Darken color by blending with black"""
    c = np.array(mcolors.to_rgb(color))
    return tuple(c * (1 - factor))

def get_aa(idx):
    """Extract amino acid from feature index like '152_G' -> 'G'"""
    return idx.split("_")[1]

def compare_biome_aas(cluster_df, biome_df, cluster_name, top_k=50, biomes=("Lake","Oceanic")):
    # --- Select top features ---
    top_feats = cluster_df.nlargest(top_k, "gini_imp")["index"].tolist()
    feat_cols = [f for f in top_feats if f in biome_df.columns]

    if not feat_cols:
        print(f"No matching features for {cluster_name}")
        return

    sub = biome_df[feat_cols + ["ecosystem_subtype"]]

    # --- Count amino acid frequencies per biome ---
    biome_freqs = {}
    for biome in biomes:
        mask = sub["ecosystem_subtype"] == biome
        counts = sub.loc[mask, feat_cols].sum()
        aa_counts = counts.rename(lambda x: get_aa(x)).groupby(level=0).sum()
        aa_counts = aa_counts.reindex(aa_order, fill_value=0)

        total = aa_counts.sum()
        if total > 0:
            aa_counts = aa_counts / total  # normalize to relative frequency
        biome_freqs[biome] = aa_counts

    df_plot = pd.DataFrame(biome_freqs).fillna(0)

    # --- Plot ---
    fig, axes = plt.subplots(2, 1, figsize=(14,12), sharex=False)

       # 1. Grouped bar chart
    bar_width = 0.4
    x = np.arange(len(aa_order))

    for biome, offset in zip(biomes, [-bar_width/2, bar_width/2]):
        values = df_plot[biome]
        if biome == "Lake":
            colors = [lighten_color(base_colors[aa], 0.6) for aa in aa_order]
        else:  # Oceanic
            colors = [darken_color(base_colors[aa], 0.3) for aa in aa_order]

        axes[0].bar(x + offset, values, width=bar_width, label=biome, alpha=0.9, color=colors)

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(aa_order, fontsize=12)
    axes[0].set_ylabel("Relative frequency")
    axes[0].set_title(f"Top {top_k} features - {cluster_name} (AA frequencies)")

    # --- Build legend manually ---
    biome_patches = [
        mpatches.Patch(color="lightgrey", label="Lake (lighter)"),
        mpatches.Patch(color="dimgray", label="Oceanic (darker)")
    ]
    group_patches = [mpatches.Patch(color=color, label=name) for name, _, color in aa_groups]

    axes[0].legend(handles=group_patches + biome_patches,
                   loc="upper right", bbox_to_anchor=(1.15, 1.05))


    # 2. Difference plot (Lake - Oceanic)
    diff = df_plot[biomes[0]] - df_plot[biomes[1]]
    diff_colors = [
        lighten_color(base_colors[aa], 0.6) if v > 0 else darken_color(base_colors[aa], 0.3)
        for aa, v in zip(aa_order, diff)
    ]

    axes[1].bar(x, diff, color=diff_colors, alpha=0.9)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(aa_order, fontsize=12)
    axes[1].set_ylabel(f"{biomes[0]} - {biomes[1]}")
    axes[1].set_title("Difference in normalized frequencies (enrichment)")

    plt.tight_layout()
    plt.show()


compare_biome_aas(feat_df, df_biome, cluster_num, top_k=50)



