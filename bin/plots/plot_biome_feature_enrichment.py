import pandas as pd
import logomaker
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import seaborn as sns

def make_matrices(cluster_num, biome1, biome2, x=50):

    # read in the data_df for the cluster
    data_df_merged_balanced = pd.read_csv(f"RF_results/raw/{cluster_num}_data_df.tsv", sep='\t')
    data_df_merged_balanced['biome'] = data_df_merged_balanced['biome'].str.replace(":",";")
    # print(data_df_merged_balanced.groupby('ecosystem_subtype')['179_G'].sum())

    # Read in feature importance scores
    feat_imp_df = pd.read_csv(f"RF_results/raw/{cluster_num}_Feature_importance.tsv", sep='\t', usecols=[0,2],names=['col_name','feature_importance_vals'],header=0)
    top_50_ft_df = feat_imp_df.head(n=x)

    top_50_ft_df[['position','aa']] = top_50_ft_df['col_name'].str.split("_",expand=True)
    top_50_ft_df['position'] = top_50_ft_df['position'].astype(int)


    data_df_merged_balanced_1 = data_df_merged_balanced.groupby('ecosystem_subtype').apply(lambda x: (x==1).sum()).T
    data_df_merged_balanced_1
    data_df_merged_balanced_1[f"Proportion_{biome1}"] = (
        data_df_merged_balanced_1[biome1] / data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome1].shape[0] * 100)
    data_df_merged_balanced_1[f"Proportion_{biome2}"] = (
        data_df_merged_balanced_1[biome2] / data_df_merged_balanced[data_df_merged_balanced['ecosystem_subtype'] == biome2].shape[0] * 100)

    # print(data_df_merged_balanced_1)
    data_df_merged_balanced_1 = data_df_merged_balanced_1.T

    top_50_ft_df_t = top_50_ft_df.set_index('col_name').T
    top_features = data_df_merged_balanced_1[top_50_ft_df_t.columns]
    # print(top_features.T.head(n=20))
    top_features_T = top_features.T.reset_index(names='pos_AA')
    top_features_T['Proportion_Lake'] = top_features_T['Proportion_Lake']*-1
    top_features_T_m = pd.melt(top_features_T[['pos_AA','Proportion_Lake','Proportion_Oceanic']],id_vars="pos_AA")
    top_features_T_m

    return top_features_T_m

def plot_bars(top_features_T_m, cluster_num,colour1, colour2, x=50):
    sns.set_theme(style="whitegrid")
    my_palette = sns.diverging_palette(h_neg=100, h_pos=200,s=100,l=60, n=2)
    min,max = 2, 50
    plt.figure(figsize=(10,x/2))

    b = sns.barplot(
        data=top_features_T_m,
        x="value",
        y="pos_AA",
        hue="ecosystem_subtype",
        palette=[colour1,colour2],
        dodge=False,
        width=-.5,  
    )

    # show the graph
    b.set_title(f"{cluster_num} Top{x} features",fontsize=14)
    b.set_xlabel("Proportion of feature in sequences from either biome (%)",fontsize=14)
    b.set_ylabel(f"Top {x} Important features",fontsize=14)

    # plt.show()
    return plt

def plot(cluster_num, biome1, biome2, colour1, colour2,k):
    top_features_T_m = make_matrices(cluster_num, biome1, biome2, k)
    plt = plot_bars(top_features_T_m, cluster_num,colour1, colour2,k)
    return plt