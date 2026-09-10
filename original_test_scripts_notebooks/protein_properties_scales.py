# %%
from Bio import SeqIO
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd

# %%
# how to add a new scale for the protein scales
charge_scale = {
    'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'E': -1, 'Q': 0, 'G': 0,
    'H': 0.5, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0, '-':0
}
# isoelectric points profile (Source: Lehninger Biochemistry 5th ed.; and Expasy amino acid property tables)
# these are more for amino acids in free solution (based on empirical data) not within a peptide
pI_scale = {
    'A': 6.00, 'R': 10.76, 'N': 5.41, 'D': 2.77, 'C': 5.07,
    'E': 3.22, 'Q': 5.65, 'G': 5.97, 'H': 7.59, 'I': 6.02,
    'L': 5.98, 'K': 9.74, 'M': 5.74, 'F': 5.48, 'P': 6.30,
    'S': 5.68, 'T': 5.60, 'W': 5.89, 'Y': 5.66, 'V': 5.96,'-':5.74
}

# scale from the isoelectricpoint module - it is calculated at different pHs therefore dynamic. It is more accurate to use for amino acids within peptides
positive_pKs = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
negative_pKs = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
pKcterminal = {"D": 4.55, "E": 4.75}
pKnterminal = {
    "A": 7.59,
    "M": 7.0,
    "S": 6.93,
    "P": 8.36,
    "T": 6.82,
    "V": 7.44,
    "E": 7.7,
}
charged_aas = ("K", "R", "H", "D", "E", "C", "Y")


pI_scale_dynamic = {
    'A': 5.57, 'R': 9.75, 'N': 5.53, 'D': 4.30, 'C': 5.52,
    'E': 4.6, 'Q': 5.53, 'G': 5.53, 'H': 6.74, 'I': 5.53,
    'L': 5.53, 'K': 8.75, 'M': 5.28, 'F': 5.53, 'P': 5.95,
    'S': 5.24, 'T': 5.18, 'W': 5.53, 'Y': 5.52, 'V': 5.50, '-': 5.28
}

# %%
# all_prot_features = {}
# with open ("RF_results/cluster_fastas/cluster_0_1_mixed_copy.faa", 'r') as fasta_in:
#     records = SeqIO.parse(fasta_in, "fasta")
#     for record in records:
#         # print(record)
#         X = ProteinAnalysis(str(record.seq))
#         all_prot_features[record.id] = [X.count_amino_acids(),X.amino_acids_percent,X.molecular_weight(),
#                                         X.aromaticity(),X.instability_index(),X.isoelectric_point(),X.charge_at_pH(7), X.gravy(),
#                                         X.secondary_structure_fraction(),X.molar_extinction_coefficient(),
#                                         X.protein_scale(ProtParamData.Flex, window=9),X.protein_scale(ProtParamData.pa, window=9), X.protein_scale(charge_scale, window=9),X.protein_scale(pI_scale, window=9)]
#         # polarity_profile = ProteinAnalysis(str(record.seq)).protein_scale(ProtParamData.pa, window=15)
#         # charge_profile = ProteinAnalysis(str(record.seq)).protein_scale(charge_scale, window=15)
#         #X.protein_scale(vihinen_flexibility, window=20)

# %%
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
protein = IP("V")
print(f"IEP of peptide {protein.sequence} is {protein.pi():.2f}")
print(f"Its charge at pH 7 is {protein.charge_at_pH(7.0):.2f}")


# %%
all_prot_features = {}
with open ("RF_results/foldmason/cluster_phyloglm_db_out/cluster_phyloglm_MSTA_aa.fa", 'r') as fasta_in:
    records = SeqIO.parse(fasta_in, "fasta")
    for record in records:
        # print(record)
        X = ProteinAnalysis(str(record.seq))
        all_prot_features[record.id] = [X.protein_scale(ProtParamData.pa, window=9), X.protein_scale(charge_scale, window=9),X.protein_scale(pI_scale_dynamic, window=9),X.protein_scale(ProtParamData.Flex, window=9)]
        # polarity_profile = ProteinAnalysis(str(record.seq)).protein_scale(ProtParamData.pa, window=15)
        # charge_profile = ProteinAnalysis(str(record.seq)).protein_scale(charge_scale, window=15)
        #X.protein_scale(vihinen_flexibility, window=20)

# %%
# all_prot_features['MGYP003333944615'][9]
# for k,v in all_prot_features.items():
#     print(len(v[7]))

# %%
# df_analysis = pd.DataFrame.from_dict(all_prot_features, orient='index',
#                             columns=['count_amino_acids','amino_acids_percent','molecular_weight',
#                                      'aromaticity','instability_index','isoelectric_point','charge_at_pH7','gravy',
#                                      'secondary_structure_fraction','molar_extinction_coefficient','flexibility', 'polarity_profile', 'charge_profile','pI_profile'])
# df_analysis

df_analysis = pd.DataFrame.from_dict(all_prot_features, orient='index',
                            columns=['gravy_profile', 'charge_profile','iso_profile','flexibility_profile'])
df_analysis

# %% [markdown]
# Add in information about the ecosystems for the sequences

# %%
data_df = pd.read_csv("RF_results/raw/cluster_0_1_mixed_data_df.tsv",usecols=['protein','ecosystem_subtype'],sep='\t').set_index('protein')
data_df

# %%
# merge dataframes
df = df_analysis.join(data_df,how='left')
df

# %% [markdown]
# ## Plotting
# ### Plotting the tiny plots with averages

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from itertools import combinations
from scipy.stats import ttest_ind
from scipy.stats import mannwhitneyu


# %%
# Compute the mean flexibility per row
df["median_flexibility"] = df["flexibility"].apply(
    lambda x: np.nanmedian(x) if isinstance(x, (list, np.ndarray)) else np.nan)

# --- 1. Expand nested columns ---
def expand_dict_column(df, col):
    expanded = pd.json_normalize(df[col])
    expanded.columns = [f"{col}_{c}" for c in expanded.columns]
    expanded.index = df.index
    return expanded

def expand_list_column(df, col):
    expanded = pd.DataFrame(df[col].tolist(), index=df.index)
    expanded.columns = [f"{col}_{i}" for i in range(len(expanded.columns))]
    expanded.index = df.index
    return expanded

dict_cols = ["count_amino_acids", "amino_acids_percent"]
list_cols = ["flexibility", "secondary_structure_fraction"]

expanded_parts = []

for col in dict_cols:
    expanded = expand_dict_column(df, col)
    expanded_parts.append(expanded.reset_index(drop=True))

for col in list_cols:
    expanded = expand_list_column(df, col)
    expanded_parts.append(expanded.reset_index(drop=True))

# Reset df index to make sure it's unique
df_reset = df.reset_index(drop=True)

# Safe concatenation
df_expanded = pd.concat(expanded_parts + [df_reset[["ecosystem_subtype"]]], axis=1)

def add_ttest_annotations(ax, df, feature, group_col="ecosystem_subtype"):
    groups = df[group_col].unique()
    pairs = list(combinations(groups, 2))
    
    y_max = df[feature].max()
    y_min = df[feature].min()
    y_range = y_max - y_min
    offset = y_range * 0.05  # vertical spacing for annotations
    current_y = y_max + offset

    for g1, g2 in pairs:
        data1 = df.loc[df[group_col] == g1, feature].dropna()
        data2 = df.loc[df[group_col] == g2, feature].dropna()
        if len(data1) < 2 or len(data2) < 2:
            continue
        stat, pval = ttest_ind(data1, data2, equal_var=False)  # Welch's t-test
        # stat, pval = mannwhitneyu(data1, data2,) # Mann-Whitney U test for non-parametric (e.g. gravy index is not normally idstributed)
        # Add annotation
        if pval < 0.001:
            text = "***"
        elif pval < 0.01:
            text = "**"
        elif pval < 0.05:
            text = "*"
        else:
            text = "ns"
        x1, x2 = groups.tolist().index(g1), groups.tolist().index(g2)
        ax.plot([x1, x1, x2, x2], [current_y, current_y + offset/2, current_y + offset/2, current_y], color='k')
        ax.text((x1 + x2) / 2, current_y + offset/2, text, ha='center', va='bottom', fontsize=8)
        current_y += offset

def plot_pca(data, label, title, ax, colour1,colour2,show_arrows=False, n_arrows=5):
    X = data.select_dtypes(include=[np.number])

    # Drop rows with NaNs and align labels
    mask = ~X.isna().any(axis=1)
    X = X[mask]
    label = label[mask]

    if X.shape[0] < 2 or X.shape[1] < 2:
        ax.set_title(f"{title} (not enough data for PCA)")
        ax.axis("off")
        return

    # Standardize
    X_scaled = StandardScaler().fit_transform(X)

    # PCA
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_scaled)
    pcs_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'])
    pcs_df['ecosystem_subtype'] = label.values

    # Plot samples
    sns.scatterplot(
        data=pcs_df,
        x='PC1', y='PC2',
        hue='ecosystem_subtype',
        palette=[colour1,colour2],
        ax=ax,
        s=40, alpha=0.8
    )

    ax.set_title(
        f"{title} (PCA: {pca.explained_variance_ratio_[0]*100:.1f}% + {pca.explained_variance_ratio_[1]*100:.1f}%)"
    )

    # Optional: add arrows for loadings
    if show_arrows:
        loadings = pca.components_.T  # shape: features x PCs
        feature_names = X.columns

        # Scale arrows for visibility
        arrow_scale = 3 * np.mean(np.abs(pcs)) / np.mean(np.abs(loadings))

        # Pick top n features by length of vector
        magnitudes = np.linalg.norm(loadings[:, :2], axis=1)
        top_idx = np.argsort(magnitudes)[-n_arrows:]

        for i in top_idx:
            x, y = loadings[i, 0] * arrow_scale, loadings[i, 1] * arrow_scale
            ax.arrow(0, 0, x, y, color='black', alpha=0.7, width=0.002, head_width=0.1)
            ax.text(x * 1.1, y * 1.1, feature_names[i][-1], fontsize=8, color='black', ha='center', va='center')

    ax.legend(fontsize=7, loc='best')


# --- 3. Prepare figure grid ---
fig, axes = plt.subplots(3, 4, figsize=(15, 10))
axes = axes.ravel()

# --- 4. PCA plots (4 of them) ---
# lake: "#9dcc7e" ocean:"#2B8CBF"
plot_pca(expand_dict_column(df, "count_amino_acids"), df["ecosystem_subtype"], 
         "Count amino acids", axes[0],"#9dcc7e","#2B8CBF", show_arrows=True)
plot_pca(expand_dict_column(df, "amino_acids_percent"), df["ecosystem_subtype"], 
         "Amino acids %", axes[1], "#9dcc7e","#2B8CBF", show_arrows=True)
# plot_pca(expand_list_column(df, "flexibility"), df["ecosystem_subtype"], "Flexibility", axes[2])
plot_pca(expand_list_column(df, "secondary_structure_fraction"), df["ecosystem_subtype"], "Secondary structure", axes[3], "#9dcc7e","#2B8CBF",show_arrows=True)

# --- 5. Violin plots (6 of them) ---
violin_features = [
    "molecular_weight", "aromaticity", "instability_index",
    "isoelectric_point", "charge_at_pH7", "gravy", "median_flexibility"
]

for i, feature in enumerate(violin_features):
    ax = axes[4 + i]
    sns.violinplot(data=df, x="ecosystem_subtype", y=feature, ax=ax, inner="quartile", scale="width", palette=["#9dcc7e","#2B8CBF"])
    ax.set_title(feature.replace("_", " ").capitalize())
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend_ = None
    add_ttest_annotations(ax, df, feature)


# --- 6. Final touches ---
for j in range(len(violin_features) + 4, 10):
    fig.delaxes(axes[j])  # remove any unused subplots

plt.tight_layout()
plt.show()


# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind

plt.rcParams.update({'font.size': 12})

# define aa groups once
AA_GROUPS = [
    ("Positive charged", ["R", "H", "K"], "crimson"),
    ("Negative charged", ["D", "E"], "royalblue"),
    ("Polar uncharged", ["S", "T", "N", "Q"], "seagreen"),
    ("Hydrophobic", ["A", "V", "I", "L", "M"], "goldenrod"),
    ("Aromatic", ["F", "Y", "W"], "purple"),
    ("Special cases", ["C", "U", "G", "P"], "darkorange")
]
# flat map AA -> color
AA_TO_COLOR = {aa: color for _, aas, color in AA_GROUPS for aa in aas}


def expand_dict_column(df, col):
    expanded = pd.json_normalize(df[col])
    expanded.columns = [f"{col}_{c}" for c in expanded.columns]
    expanded.index = df.index
    return expanded

def expand_list_column(df, col):
    expanded = pd.DataFrame(df[col].tolist(), index=df.index)
    expanded.columns = [f"{col}_{i}" for i in range(len(expanded.columns))]
    expanded.index = df.index
    return expanded

def add_ttest_annotations(ax, df, feature, group_col="ecosystem_subtype"):
    groups = df[group_col].unique()
    pairs = list(combinations(groups, 2))
    y_max = df[feature].max()
    y_min = df[feature].min()
    y_range = (y_max - y_min) if (y_max - y_min) != 0 else 1.0
    offset = y_range * 0.06  # vertical spacing for annotations
    current_y = y_max + offset
    for g1, g2 in pairs:
        data1 = df.loc[df[group_col] == g1, feature].dropna()
        data2 = df.loc[df[group_col] == g2, feature].dropna()
        if len(data1) < 2 or len(data2) < 2:
            continue
        stat, pval = ttest_ind(data1, data2, equal_var=False)
        if pval < 0.001:
            text = "***"
        elif pval < 0.01:
            text = "**"
        elif pval < 0.05:
            text = "*"
        else:
            text = "ns"
        x1, x2 = groups.tolist().index(g1), groups.tolist().index(g2)
        ax.plot([x1, x1, x2, x2], [current_y, current_y + offset/2, current_y + offset/2, current_y], color='k')
        ax.text((x1 + x2) / 2, current_y + offset/2, text, ha='center', va='bottom', fontsize=8)
        current_y += offset

def _extract_aa_from_feature(feature_name):
    """
    Try to obtain an amino acid single-letter code from a feature name.
    Works for names like 'count_amino_acids_A' or 'amino_acids_percent_A' or '..._A'.
    Fallback: return last character if it is an uppercase letter A-Z, otherwise '?'
    """
    parts = str(feature_name).split('_')
    candidate = parts[-1]
    if len(candidate) == 1 and candidate.isalpha() and candidate.isupper():
        return candidate
    # sometimes last char of the string is AA
    if len(feature_name) >= 1 and feature_name[-1].isalpha() and feature_name[-1].isupper():
        return feature_name[-1]
    return None

def plot_pca_with_colored_arrows(data, label, title, ax, colour1, colour2, aa_groups=AA_GROUPS, n_arrows=21, top_n_color=5):
    """
    PCA scatterplot + arrows (loadings). Arrows for top_n_color loadings are colored
    by their AA group, the rest are grey.
    """
    X = data.select_dtypes(include=[np.number]).copy()
    # drop rows with NaNs
    mask = ~X.isna().any(axis=1)
    X = X.loc[mask]
    label = label.loc[mask]

    if X.shape[0] < 2 or X.shape[1] < 2:
        ax.set_title(f"{title} (not enough data for PCA)")
        ax.axis("off")
        return

    # standardize and PCA
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X_scaled)
    pcs_df = pd.DataFrame(pcs, columns=['PC1', 'PC2'])
    pcs_df['ecosystem_subtype'] = label.values

    sns.scatterplot(
        data=pcs_df, x='PC1', y='PC2',
        hue='ecosystem_subtype', palette=[colour1, colour2],
        ax=ax, s=40, alpha=0.9, legend=False
    )

    ax.set_title(f"{title}\n(PCA: {pca.explained_variance_ratio_[0]*100:.1f}% + {pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.axvline(0, color='k', lw=0.4, alpha=0.6); ax.axhline(0, color='k', lw=0.4, alpha=0.6)

    # loadings (features x PCs)
    loadings = pca.components_.T  # shape: n_features x 2
    feature_names = X.columns.to_list()
    # compute magnitudes and pick top n_arrows by magnitude for display arrows (but we'll color only top_n_color)
    magnitudes = np.linalg.norm(loadings[:, :2], axis=1)
    # pick indices to draw arrows for (largest magnitudes). If fewer features than n_arrows, draw all.
    draw_idx = np.argsort(magnitudes)[-min(n_arrows, len(magnitudes)):]

    # identify which of drawn arrows are the top_n_color by magnitude (global for this PCA)
    top_idx = np.argsort(magnitudes)[-min(top_n_color, len(magnitudes)):]

    # arrow scaling — scale arrows to the scatter cloud size
    # choose a scale factor that makes arrows visible but not too large
    # base scale on variance of PCs and average loading magnitude
    avg_pc = np.mean(np.abs(pcs))
    avg_load = np.mean(np.abs(loadings[:, :2]))
    if avg_load == 0:
        arrow_scale = 1.0
    else:
        arrow_scale = 2.5 * avg_pc / (avg_load + 1e-9)

    # map AA letter to color
    aa_to_color = {aa: color for _, aas, color in aa_groups for aa in aas}

    # draw arrows for chosen indices
    for i in draw_idx:
        fx, fy = loadings[i, 0] * arrow_scale, loadings[i, 1] * arrow_scale
        fname = feature_names[i]
        aa = _extract_aa_from_feature(fname)
        # color: colored only if this feature belongs to top_idx AND AA has mapping; else grey
        color = 'lightgray'
        if i in top_idx and aa is not None:
            color = aa_to_color.get(aa, 'black')
        # draw arrow
        ax.arrow(0, 0, fx, fy, color=color, alpha=0.9, width=0.0015, head_width=0.05, length_includes_head=True, zorder=3)
        # annotate arrow tip with AA if available
        if aa is not None:
            ax.text(fx * 1.12, fy * 1.12, aa, color=color if color != 'lightgray' else 'dimgray',
                    fontsize=10, ha='center', va='center', weight='bold')

    return ax


def plot_subplots(df, colour1, colour2):
    """
    Main plotting function:
      - Top row: 3 PCA plots (count_amino_acids, amino_acids_percent, secondary_structure_fraction)
        larger than the violins below.
      - Bottom rows: 6 violin plots arranged 2 x 3.
      - Arrows on PCA: only top 5 loadings colored by AA group, rest grey.
    """
    # Precompute median flexibility
    df["median_flexibility"] = df["flexibility"].apply(
        lambda x: np.nanmedian(x) if isinstance(x, (list, np.ndarray)) else np.nan
    )

    # Expand columns
    count_aa_df = expand_dict_column(df, "count_amino_acids")
    pct_aa_df = expand_dict_column(df, "amino_acids_percent")
    secstr_df = expand_list_column(df, "secondary_structure_fraction")

    # Build figure with GridSpec: top row bigger than bottom (ratio 2:1)
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 3, height_ratios=[2, 1], hspace=0.25, wspace=0.25)

    # top: three PCA axes
    ax_pca1 = fig.add_subplot(gs[0, 0])
    ax_pca2 = fig.add_subplot(gs[0, 1])
    ax_pca3 = fig.add_subplot(gs[0, 2])

    # bottom: 2 rows × 3 cols for violins using subgridspec
    bottom_sub = gs[1, :].subgridspec(2, 3, hspace=0.35, wspace=0.25)
    axs_bottom = [fig.add_subplot(bottom_sub[i, j]) for i in range(2) for j in range(3)]

    # PCA plots (show arrows colored per AA groups; only top 5 colored)
    plot_pca_with_colored_arrows(count_aa_df, df["ecosystem_subtype"], "Count amino acids", ax_pca1, colour1, colour2, aa_groups=AA_GROUPS, n_arrows=40, top_n_color=5)
    plot_pca_with_colored_arrows(pct_aa_df, df["ecosystem_subtype"], "Amino acids %", ax_pca2, colour1, colour2, aa_groups=AA_GROUPS, n_arrows=40, top_n_color=5)
    plot_pca_with_colored_arrows(secstr_df, df["ecosystem_subtype"], "Secondary structure fraction", ax_pca3, colour1, colour2, aa_groups=AA_GROUPS, n_arrows=40, top_n_color=5)

    # Violin features (6)
    violin_features = [
        "molecular_weight", "aromaticity", "instability_index",
        "isoelectric_point", "charge_at_pH_7", "gravy"
    ]

    for ax, feature in zip(axs_bottom, violin_features):
        sns.violinplot(data=df, x="ecosystem_subtype", y=feature, ax=ax, inner="quartile", scale="width", palette=[colour1, colour2])
        ax.set_title(feature.replace("_", " ").capitalize())
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.legend_ = None
        try:
            add_ttest_annotations(ax, df, feature)
        except Exception:
            # ignore annotation errors (e.g., not enough data)
            pass

    # remove unused bottom axes (if any) - but we filled exactly 6
    # final layout tweaks
    plt.tight_layout()
    plt.show()
    return fig


# %%
plot_subplots()

# %% [markdown]
# ### Plotting scales 

# %%
def plot_scales(df, scale_var, colour1, colour2, title):
    # Keep only rows with flexibility lists
    flex_df = df[df[scale_var].apply(lambda x: isinstance(x, (list, np.ndarray)))].copy()

    # Convert each list into a long-format DataFrame
    records = []
    for idx, row in flex_df.iterrows():
        values = row[scale_var]
        eco = row['ecosystem_subtype']
        for i, val in enumerate(values):
            records.append({'residue_index': i, scale_var: val, 'ecosystem_subtype': eco})

    flex_long = pd.DataFrame(records)

    import seaborn as sns
    import matplotlib.pyplot as plt

    plt.figure(figsize=(15,6))

    sns.lineplot(
        data=flex_long,
        x='residue_index',
        y=scale_var,
        hue='ecosystem_subtype',
        estimator=np.median,    # median per residue index
        ci=95,                  # 95% confidence interval
        palette={colour2,colour1}
    )

    plt.xlabel("Residue index")
    plt.ylabel(scale_var.capitalize())
    plt.title(title)
    plt.tight_layout()
    plt.show()
    return plt


# %%
plt = plot_scales(df, "polarity_profile","#9dcc7e","#2B8CBF", title=f"Median GRAVY index defined by Kyte & Doolittle hydropathy with 95% CI,\n sliding window = 9")

# %%
plt = plot_scales(df, "iso_profile","#9dcc7e","#2B8CBF", title=f"Median Isoelectric point with 95% CI,\n sliding window = 9")

# %%
plt = plot_scales(df, "charge_profile","#9dcc7e","#2B8CBF", title=f"Median charge profile with 95% CI,\n sliding window = 9")

# %%
plt = plot_scales(df, "flexibility_profile","#9dcc7e","#2B8CBF", title=f"Median flexibility profile with 95% CI,\n sliding window = 9")

# %%
def plot_scale_difference(df, scale_var, eco1, eco2, color, title):
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    # Filter out rows that actually have scale lists
    flex_df = df[df[scale_var].apply(lambda x: isinstance(x, (list, np.ndarray)))].copy()

    # Expand to long format
    records = []
    for _, row in flex_df.iterrows():
        values = row[scale_var]
        eco = row['ecosystem_subtype']
        for i, val in enumerate(values):
            records.append({'residue_index': i, scale_var: val, 'ecosystem_subtype': eco})

    flex_long = pd.DataFrame(records)

    # Compute median per residue index per ecosystem
    median_by_eco = (
        flex_long
        .groupby(['ecosystem_subtype', 'residue_index'])[scale_var]
        .median()
        .unstack('ecosystem_subtype')
    )

    # Compute the difference: eco2 - eco1
    median_by_eco['difference'] = median_by_eco[eco2] - median_by_eco[eco1]

    # Optional: compute rolling average to smooth the line (e.g., window=5)
    median_by_eco['smoothed_diff'] = median_by_eco['difference'].rolling(window=10, center=True, min_periods=1).mean()

    # Plot
    plt.figure(figsize=(15,6))
    sns.lineplot(
        x=median_by_eco.index,
        y=median_by_eco['smoothed_diff'],
        color=color,
        linewidth=2
    )
    plt.axhline(0, color='gray', linestyle='--', linewidth=1)
    plt.xlabel("Residue index")
    plt.ylabel(f"Δ {scale_var.capitalize()} ({eco2} - {eco1})")
    plt.title(title)
    plt.tight_layout()
    return plt
    # plt.show()


# %%
plot_scale_difference(
    df,
    scale_var='charge_profile',
    eco1='Lake',
    eco2='Oceanic',
    color='teal',
    title='Difference in GRAVY median per residue (Oceanic - Lake)'
)

# plt = plot_scales(df, "pI_profile","#9dcc7e","#2B8CBF", title=f"Median Isoelectric point with 95% CI,\n sliding window = 9")


