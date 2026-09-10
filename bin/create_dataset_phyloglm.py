import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Create dataset for PhyloGLM")
parser.add_argument("feature_df", help="Path to the feature.tsv file")
parser.add_argument("biome_file", help="Path to the Biome file")
parser.add_argument("-o", "--output", default="output_dataset.tsv", help="Output file path")
args = parser.parse_args()

feat_imp_df = pd.read_csv(args.feature_df, sep='\t', usecols=[0,1],names=['col_name','feature_importance_vals'], dtype={'col_name':'str','feature_importance_vals':np.float64},header=0) 
#print(feat_imp_df.shape)
x = len(feat_imp_df)
top_50_ft_df = feat_imp_df.head(n=x)
#top_50_ft_df

# feat_imp_df.loc[feat_imp_df['feature_importance_vals'] >= 0.00001]

df_data = pd.read_csv(args.biome_file, sep='\t')
#print(df_data)
df_data.drop_duplicates('protein',inplace=True)

df_biome = df_data[["protein","ecosystem_subtype"]]
# df_biome.to_csv("comb_cluster_biome.tsv", sep='\t', index=False)
df_biome_num = pd.get_dummies(df_biome.set_index("protein"), dtype=int)
#print(df_biome_num)

# get all values of ecosystem_subtype, choose the one which is not Lake
unique_biomes = df_biome['ecosystem_subtype'].dropna().unique()
non_lake_biomes = [b for b in unique_biomes if b != "Lake"]
target_biome = non_lake_biomes[0]
target_col = f"ecosystem_subtype_{target_biome}"

df_data_filtered = df_data[["protein"] + [col for col in top_50_ft_df['col_name'] if col in df_data.columns]]
df_data_filtered_merged = df_biome_num[[target_col]].join(df_data_filtered.set_index("protein"),how="left")
df_data_filtered_merged.reset_index(inplace=True)
#print(df_data_filtered_merged)

# save look up table with feature and gini
gini_lookup = top_50_ft_df.rename(columns={'col_name': 'feature', 'feature_importance_vals': 'gini'})
gini_output = args.output.replace('_phyloglm_dataset.tsv', '_gini.tsv')
gini_lookup.to_csv(gini_output, sep='\t', index=False)

df_data_filtered_merged.to_csv(args.output, sep='\t', index=False)

# df_phylo_res = pd.read_csv("/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/Microviridae_analysis/multiple_cluster_analysis/combined_model_phyloglm/phyloglm_results_combined_model_pruned_full_new.csv",index_col=0).set_index("feature")

# df_phylo_res

# df_phylo_res_filt = df_phylo_res.reset_index()
# # df_phylo_res_filt['abs_estimate'] = abs(df_phylo_res_filt['Estimate'])
# # df_phylo_res_filt = df_phylo_res_filt.loc[df_phylo_res_filt['z.value'] > 0]     # change to get positive or negative zvalues for oceanic/lake respectively

# df_phylo_res_filt['abs_z.value'] = abs(df_phylo_res_filt['z.value'])
# # df_phylo_res_filt = df_phylo_res_filt.loc[df_phylo_res_filt['padj'] <= 0.01]
# df_phylo_res_filt
# df_phylo_res_filt = df_phylo_res_filt[['feature', 'abs_z.value']]
# df_phylo_res_filt.sort_values('abs_z.value',ascending=False,inplace=True)
# df_phylo_res_filt
# df_phylo_res_filt.fillna(0,inplace=True)
# df_phylo_res_filt['fake_shap'] = 0
# df_phylo_res_filt['cluster_num'] = 'cluster_0_1_mixed_MSTA_aa_plus_ec6'

# df_phylo_res_filt.rename(columns={'feature':'index',
#                                   'abs_z.value' : 'gini_imp'},inplace=True)
# df_phylo_res_filt = df_phylo_res_filt[['index','fake_shap','gini_imp','cluster_num']]
# #print(df_phylo_res_filt)

# feat_imp_df_2 = feat_imp_df.reset_index(names=['imp_order'])

# df_merged = df_phylo_res.join(feat_imp_df_2.set_index("col_name"),how='outer')
# df_merged['abs_z.value'] = abs(df_merged['z.value'])
# df_merged['abs_Estimate'] = abs(df_merged['Estimate'])

# df_merged = df_merged.loc[df_merged['padj'] <= 0.01]
# df_merged['biome_predictor'] = np.where(df_merged['Estimate'] <0, "lake","ocean")
# df_merged.sort_values(by='feature_importance_vals',ascending=False,inplace=True)
# df_merged.drop('imp_order',axis=1,inplace=True)
# df_merged.reset_index(names="feature",inplace=True)
# df_merged.reset_index(names="imp_order",inplace=True)
# df_merged.to_csv("phyloglm_intersect_giniscores_new.tsv", sep='\t')        # intersection of phyloglm and gini scores
# df_merged_top_200 = df_merged.sort_values('feature_importance_vals',ascending=False).head(n=200)
# #print(df_merged_top_200)

# plt.figure(figsize=(15,15))
# ax = sns.scatterplot(df_merged, x='feature_importance_vals', y='abs_z.value', hue='biome_predictor')

# def label_point(x, y, val, ax):
#     # Convert val (Index) to Series so pd.concat works
#     a = pd.concat({'x': x, 'y': y, 'val': pd.Series(val, index=x.index)}, axis=1)
#     for i, point in a.iterrows():
#         ax.text(point['x'], point['y'], str(point['val']))

# label_point(df_merged['feature_importance_vals'], df_merged['abs_z.value'], df_merged.index, ax)
# #plt.show()