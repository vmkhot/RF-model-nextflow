# %% [markdown]
# # Random Forest Classifier for Predicting Biomes (AA SEQUENCES)
# 
# Notebook for exploring amino acid sites as "important" features in predicting biomes using random forest
# 
# You need a cluster of proteins --> MSA of this cluster --> trimmed using TrimAl so sequences with a lot of gaps or positions with a lot of gaps are removed

# %%
import pandas as pd
from Bio import SeqIO
import numpy as np

# %% [markdown]
# ## Data Inputs

# %%
fasta = '/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/RF_models_automation/foldmason/cluster_mixed_model_db_out/cluster_mixed_model_MSTA_aa.fa'

seq_dict = {} 

for record in SeqIO.parse(fasta, "fasta"):
    seq_dict[record.id]=record.seq

seq_dict

# %%
# convert to a dataframe
data_df = pd.DataFrame.from_dict(seq_dict,orient='index').reset_index(names='protein')
data_df.loc[data_df['protein'].str.contains('IMG'), 'img_split'] = (data_df.loc[data_df['protein'].str.contains('IMG'), 'protein']
      .str.split('|').str[0])
data_df['img_split'] = data_df['img_split'].fillna(data_df['protein'])
data_df.set_index('img_split', inplace=True)
data_df.drop('protein', axis=1, inplace=True)
data_df

# %%
# check the cardinality of the data
# i.e. how much variation is there at each position?
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy


def shannon_entropy(column):
    counts = column.value_counts()  # Count occurrences of each unique character
    probs = counts / counts.sum()   # Convert counts to probabilities
    return entropy(probs, base=2)   # Shannon entropy formula

# Compute entropy for each column
entropy_per_column = data_df.apply(shannon_entropy, axis=0)

# Print results
print(entropy_per_column)


sns.lineplot(x=entropy_per_column.index, y=entropy_per_column)

# Add labels to your graph
plt.xlabel('MSA Position')
plt.ylabel('Shannon Entropy')
plt.title("Shannon entropy across sequence position")
plt.legend()
plt.show()



# %%
data_df_numerical = pd.get_dummies(data_df, dtype=int)
data_df_numerical = data_df_numerical.loc[:,~data_df_numerical.columns.str.contains('-')]       # remove the position-gap columns

data_df_numerical

# %%
# check variance
from sklearn.feature_selection import VarianceThreshold

selector = VarianceThreshold(0.9 * (1 - 0.9))
selector.fit_transform(data_df_numerical)

variance_array = selector.variances_
variance_array.sort()
print(variance_array)

# %% [markdown]
# ## Biomes input and filtering

# %%
biomes_metadata_df = pd.read_csv('/Users/varadakhot/Library/CloudStorage/OneDrive-Friedrich-Schiller-UniversitätJena/MCP_struct/amino_acid_analyis/imgvr_mgnify_combined/best_biome_per_mcp_protein.csv', header=0, sep=',')

biomes_metadata_df

# %%
### Check the overall enrichment of AA in each of biomes chosen. 
data_df_melted = pd.melt(data_df.reset_index(),id_vars='img_split')
data_df_melted

data_df_grouped = data_df_melted.groupby(['img_split','value'])['value'].count().reset_index(level=1,name='Count_of_AA')
data_df_grouped

data_df_grouped_w_biomes = data_df_grouped.join(biomes_metadata_df.set_index('protein'), how='left')
data_df_grouped_w_biomes

# Reset index to ensure img_split is a column


# %%
# del min

# df_reset = data_df_grouped_w_biomes.reset_index()

# # Sample 1000 unique img_split IDs per ecosystem_subtype
# sampled_img_split_ids = (
#     df_reset.groupby("ecosystem_subtype")["img_split"]
#     .apply(lambda x: x.drop_duplicates().sample(n=min(76, len(x)), random_state=42))
#     .values
# )

# # Filter the original dataframe to keep all rows for those img_split IDs
# df_sampled = df_reset[df_reset["img_split"].isin(sampled_img_split_ids)]

# # Optional: reset index for cleaner display or further use
# df_sampled = df_sampled.reset_index()
# df_sampled
# df_merged_grouped_by_biome = df_sampled.groupby(['ecosystem_subtype','value'])['Count_of_AA'].sum().reset_index()
# df_merged_grouped_by_biome.sort_values(['ecosystem_subtype','Count_of_AA'], ascending=[True, False], inplace=True)
# df_merged_grouped_by_biome

# # Assuming your DataFrame is called df
# pivot = df_merged_grouped_by_biome.pivot(index="value", columns="ecosystem_subtype", values="Count_of_AA")

# # Calculate difference
# pivot["Lake_minus_Hypersaline"] = pivot["Lake"] - pivot["Hot (42-90C)"]

# # Reset index so 'value' becomes a column again
# pivot_reset = pivot.reset_index()

# # Melt back to long format
# df_long = pivot_reset.melt(id_vars=["value", "Lake_minus_Hypersaline"], 
#                            value_vars=["Lake", "Hot (42-90C)"],
#                            var_name="ecosystem_subtype", 
#                            value_name="Count_of_AA")

# # Sort by Lake_minus_Hypersaline
# df_long_sorted = df_long.sort_values("Lake_minus_Hypersaline", ascending=True)
# df_long_sorted

# sns.set_theme(style="whitegrid")
# my_palette = sns.diverging_palette(h_neg=100, h_pos=200,s=100,l=60, n=2)
# min,max = 2, 50
# plt.figure(figsize=(6,10))

# b = sns.barplot(
#     data=df_long_sorted,
#     x="Count_of_AA",
#     y="value",
#     hue="ecosystem_subtype",
#     palette=my_palette,
#     dodge=True,
#     width=-.5,  
# )

# # show the graph
# # b.set_yticklabels(top_200_features_melted.variable)
# b.set_xlabel("Amino acid",fontsize=15)
# b.set_ylabel("count_of_amino acids in by biome",fontsize=15)
# plt.show()


# %%
data_df_merged = data_df_numerical.join(biomes_metadata_df.set_index('protein'), how='left') #.dropna(subset=['ecosystem_subtype'])
# print(data_df_merged)
temp = data_df_merged.groupby('ecosystem_subtype').count() # only to count the number of samples in each class before subsampling
temp


# %%
data_df_merged_balanced = data_df_merged[~data_df_merged['ecosystem_type'].isin(['Landfill','Non-marine Saline and Alkaline','Activated Sludge','Defined media','Digestive system','Industrial wastewater','Terephthalate','Water and sludge','Nutrient removal','Sediment','Estuary','Lentic','Aquaculture'])]
data_df_merged_balanced

temp = data_df_merged_balanced.groupby('ecosystem_type').count() # only to count the number of samples in each class before subsampling
print(temp)


# %%
data_df_merged_balanced.loc[data_df_merged_balanced['ecosystem_type'] == 'Thermal springs']

# %%
# MGYP000612766211
# data_df_merged_balanced.loc['MGYP000612766211']

# %%
# data_df_merged_balanced.to_csv('data_df_merged_balanced_lake_hot.csv')

# %%
# sample randomly 1000 from each class
max_size = data_df_merged_balanced['ecosystem_type'].value_counts().max()

data_df_merged_balanced = data_df_merged_balanced.groupby('ecosystem_type').sample(n=max_size, replace=True)
# data_df_merged_balanced.to_csv('data_df_merged_balanced_lake_hot_0_1_2_nonloop.csv',sep='\t')


# %%
print(data_df_numerical.shape)
print(data_df_merged_balanced.shape)

# %% [markdown]
# ## RF Classifier and Accuracy Scores

# %%
# Import train_test_split function
from sklearn.model_selection import train_test_split

X=data_df_merged_balanced.iloc[:,0:data_df_numerical.columns.shape[0]] # Features, features are the number of columns in the data_numerical
print(X)
y=data_df_merged_balanced['ecosystem_type']  # Labels
y = pd.get_dummies(y, dtype=int)
print(y)
# Split dataset into training set and test set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3,stratify=y) # 70% training and 30% test

# %%
# #Import Random Forest Model
# from sklearn import metrics
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.metrics import confusion_matrix
# #Create a Gaussian Classifier
# # n_estimators is the number of the decision trees in the forest. Increasing number generally improves the performance of the model but also increases the computational cost
# clf=RandomForestClassifier(n_estimators=100,max_depth=2, class_weight='balanced')

# #Train the model using the training sets y_pred=clf.predict(X_test)
# clf.fit(X_train,y_train)

# # test the model using the 30% test dataset
# y_pred=clf.predict(X_test)
# print("Accuracy:",metrics.accuracy_score(y_test, y_pred))


# %% [markdown]
# ### Accuracy scores
# 
# #### With training dataset
# Check for overfitting (when accuracy is 100%)

# %% [markdown]
# ### Hyperparameter fine tuning
# {'bootstrap': True,
#  'ccp_alpha': 0.0,
#  'class_weight': 'balanced_subsample',
#  'criterion': 'gini',
#  'max_depth': 5,
#  'max_features': 'sqrt',
#  'max_leaf_nodes': None,
#  'max_samples': None,
#  'min_impurity_decrease': 0.0,
#  'min_samples_leaf': 1,
#  'min_samples_split': 2,
#  'min_weight_fraction_leaf': 0.0,
#  'monotonic_cst': None,
#  'n_estimators': 10000,
#  'n_jobs': None,
#  'oob_score': False,
#  'random_state': None,
#  'verbose': 0,
#  'warm_start': False}

# %%
from sklearn.model_selection import RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier

param_grid = {
 'max_depth': [5, 10, None],
 'max_features': ['auto', 'sqrt', 'log2'],
 'min_samples_leaf': [1, 2, 4],
 'min_samples_split': [2, 5, 10],
 'n_estimators': [500,1000, 2000, 5000],
 'class_weight': ['balanced',None]}

random_search = RandomizedSearchCV(RandomForestClassifier(),
                                   param_grid)
random_search.fit(X_train, y_train)
print(random_search.best_estimator_)


# %%
#Import scikit-learn metrics module for accuracy calculation
from sklearn import metrics
from sklearn.metrics import confusion_matrix

# updated model check
clf_opt=random_search.best_estimator_

#Train the model using the training sets y_pred=clf.predict(X_test)
clf_opt.fit(X_train,y_train)

# test the model using the 30% test dataset
y_pred_opt=clf_opt.predict(X_test)

# Model Accuracy, how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_test, y_pred_opt))
print("Classification report:",metrics.classification_report(y_test, y_pred_opt))
print("confusion matrix: \n", confusion_matrix(y_test.values.argmax(axis=1), y_pred_opt.argmax(axis=1)))

# %%
report = metrics.classification_report(y_test, y_pred_opt,output_dict=True)
report['0']['f1-score']


# %%
# check for overfitting

y_pred_tr_opt=clf_opt.predict(X_train)

# Model Accuracy, how often is the classifier correct?
print("Accuracy:",metrics.accuracy_score(y_train, y_pred_tr_opt))
print("Classification report:",metrics.classification_report(y_train, y_pred_tr_opt))
print("confusion matrix: \n", confusion_matrix(y_train.values.argmax(axis=1), y_pred_tr_opt.argmax(axis=1)))

# %% [markdown]
# ## Feature importance
# 
# ### GINI importance

# %%

import pandas as pd
feature_imp = pd.Series(clf_opt.feature_importances_,index=data_df_numerical.columns).sort_values(ascending=False)
print(feature_imp.head(n=20))
print(feature_imp.head(n=20).sum())


# %%
import matplotlib.pyplot as plt
import seaborn as sns
#matplotlib inline

plt.figure(figsize=(10,6))
# Creating a bar plot
sns.histplot(feature_imp, bins = 100, stat='count', log_scale=(False,True))

# Add labels to your graph
plt.xlabel('Feature Importance Score')
plt.ylabel('Count of features')
plt.title("Histo of Features using GINI scores")
plt.legend()
plt.show()

# %%
# temp = data_df_merged_balanced.loc[data_df_merged_balanced['333_K'] == 1]
# temp
print(data_df_merged_balanced.shape)

# %% [markdown]
# ### SHAP Values

# %%
# SHAP values
import shap

# Create the SHAP explainer for the Random Forest model
explainer = shap.TreeExplainer(clf_opt)

# Calculate SHAP values for the test set
shap_values = explainer.shap_values(X_test)

# print(shap_values)

def shap_values_to_list(shap_values, clf_opt):
    shap_as_list=[]
    for i in range(len(clf_opt.classes_)):
        shap_as_list.append(shap_values[:,:,i])
    return shap_as_list


# Generate a SHAP summary plot
print(shap_values.shape)
print(X_test.shape)
# shap_values_transposed = shap_values.transpose(1, 0, 2)
# shap.summary_plot(shap_values[0], X_test)
# shap.summary_plot(shap_values[:,:,0], X_test, class_names=feature_imp.index,plot_type="bar")

shap_as_list = shap_values_to_list(shap_values, clf_opt)
# print(shap_as_list)

shap.summary_plot(shap_as_list[0], X_test, class_names=['Freshwater','marine'])
# shap.waterfall_plot(explainer, X_test)
# shap.summary_plot(shap_as_list[1], X_test, class_names=['freshwater','marine'])



# %%
feature_names = X_train.columns
rf_resultX = pd.DataFrame(shap_as_list[0], columns = feature_names)

vals = np.abs(rf_resultX.values).mean(0)

shap_importance = pd.DataFrame(list(zip(feature_names, vals)),
                                  columns=['col_name','feature_importance_vals'])
shap_importance.sort_values(by=['feature_importance_vals'],
                               ascending=False, inplace=True)
shap_importance

# %%
# shap_values = explainer.shap_values(X)

# # Interaction values (expensive)
# interaction_values = explainer.shap_interaction_values(X)

# # You can visualize interactions
# # shap.summary_plot(interaction_values[1], X)  # for class 1

# # Feature labels from your input DataFrame
# feature_names = X.columns.tolist()

# sns.clustermap(interaction_matrix, 
#                xticklabels=feature_names, 
#                yticklabels=feature_names, 
#                cmap='viridis',
#                figsize=(14, 12))
# plt.title("Clustered SHAP Feature Interaction Heatmap")
# plt.show()

# %%


# %%
plt.figure(figsize=(10,6))
# Creating a bar plot
sns.histplot(shap_importance, bins = 100, stat='count', log_scale=(False,True))

# Add labels to your graph
plt.xlabel('Feature Importance Score')
plt.ylabel('Count of features')
plt.title("Histo of Features using shap values")
plt.legend()
plt.show()

# %% [markdown]
# ## Feature and biome correlation
# For top x features (based on the feature importance scores from the model) - How often do they appear in a certain biome?

# %%
# get feature importances based on GINI
x = 50
feature_imp = pd.Series(clf_opt.feature_importances_,index=data_df_numerical.columns).sort_values(ascending=False)
top_50_ft_df = feature_imp.head(n=x).to_frame().T
print(top_50_ft_df)

# %%
# get feature importances based on SHAP
x = 50
top_50_ft_df = shap_importance.head(n=x).set_index('col_name').T
print(top_50_ft_df)

# %%
# data_df_merged_balanced[top_50_ft_df.columns].to_csv(f"top{x}_features.csv")

# %%
import numpy as np

# # merge with the bigger dataframe with biomes
# data_df_merged_balanced_1 = data_df_merged_balanced.groupby('ecosystem_type').apply(lambda x: (x==1).sum()).T
# data_df_merged_balanced_1['diff'] = abs(data_df_merged_balanced_1["Thermal springs"] - data_df_merged_balanced_1['Freshwater'])

# data_df_merged_balanced_1

# data_df_merged_balanced_1 = data_df_merged_balanced_1.T
# data_df_merged_balanced_1.loc['Freshwater']*= -1
# data_df_merged_balanced_1
# # top_200_features = data_df_merged_balanced_1.columns.intersection(top_50_ft_df.columns)
# top_200_features = data_df_merged_balanced_1[top_50_ft_df.columns]
# top_200_features


import numpy as np

# merge with the bigger dataframe with biomes
data_df_merged_balanced['ecosystem_type'] = data_df_merged_balanced['ecosystem_type'].str.replace(' ', '_')     # replace a space in Thermal springs to '_'
data_df_merged_balanced_1 = data_df_merged_balanced.groupby('ecosystem_type').apply(lambda x: (x==1).sum()).T

data_df_merged_balanced_1['proportion_freshwater'] = data_df_merged_balanced_1['Freshwater']/data_df_merged_balanced.groupby('ecosystem_type').count().T.Freshwater*100*-1
data_df_merged_balanced_1['proportion_thermal'] = data_df_merged_balanced_1['Thermal_springs']/data_df_merged_balanced.groupby('ecosystem_type').count().T.Thermal_springs*100
data_df_merged_balanced_1

# print(data_df_merged_balanced_1)
data_df_merged_balanced_1 = data_df_merged_balanced_1.T
data_df_merged_balanced_1
top_200_features = data_df_merged_balanced_1[top_50_ft_df.columns]
top_200_features
# only proportions
top_200_features.drop(labels=['Freshwater','Thermal_springs'], inplace=True)


# %% [markdown]
# #### Sorting values for the plot

# %%
top_200_features.T.head(n=20)

# %%
# ORIGINAL BY FEATURE IMPORTANCE (most important features shown at the top)
top_200_features.reset_index(inplace=True)
top_200_features
# melt
top_200_features_melted = pd.melt(top_200_features, id_vars=['ecosystem_type'],ignore_index=False)
# top_200_features_melted


# LARGEST DIFFERENCE IN THE LOG2 FOLD CHANGE BETWEEN FRESHWATER AND MARINE
# top_200_features = top_200_features.sort_values(axis=1, by='diff', ascending=False)     #.drop(index=('diff'),inplace=True)
# top_200_features=top_200_features.loc[['Freshwater','Marine']]
# top_200_features.reset_index(inplace=True)
# print(top_200_features)
# # melt
# top_200_features_melted = pd.melt(top_200_features, id_vars=['ecosystem_type'],ignore_index=False)
# top_200_features_melted


# FEATURE POSITION
# top_200_features_1 = top_200_features.drop(index=('diff'))
# top_200_features_1.reset_index(inplace=True)
# top_200_features_melted = pd.melt(top_200_features_1, id_vars=['ecosystem_type'])     # melt data
# top_200_features_melted[['0','position','AA','1']]= top_200_features_melted['variable'].str.split('(\d+)_([A-Za-z]+|-)', expand=True)   #split the feature into position and the AA
# top_200_features_melted[['position']]=top_200_features_melted[['position']].apply(pd.to_numeric)        # covert the position to a numeric
# top_200_features_melted = top_200_features_melted.drop(['0','1'],axis=1).sort_values("position")     # sort by position


# BY BIOME (KINDA)
# top_200_features.drop(index=('diff'),inplace=True)
# top_200_features.reset_index(inplace=True)
# top_200_features_melted = pd.melt(top_200_features, id_vars=['ecosystem_type'])     # melt data
# top_200_features_melted = top_200_features_melted.sort_values("value")

# %%
top_200_features_melted.reset_index(drop=True, inplace=True)

# %%
print(top_200_features_melted.head(20))

# %%
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid")
my_palette = sns.diverging_palette(h_neg=200, h_pos=11,s=100,l=60, n=2)
min,max = 2, 50
plt.figure(figsize=(6,x/5))

b = sns.barplot(
    data=top_200_features_melted,
    x="value",
    y="variable",
    hue="ecosystem_type",
    palette=my_palette,
    dodge=False,
    width=-.5,  
)

# show the graph
# b.set_yticklabels(top_200_features_melted.variable)
b.set_xlabel("Ecosystem",fontsize=15)
b.set_ylabel(f"Top {x} Features",fontsize=15)
plt.show()

# %% [markdown]
# ## Create a sequence logo of the top (x) positions

# %%
top_200_features_melted

# %%
shap_importance

# %%
# sequence logo creation
import logomaker

# top_200_features_melted = pd.melt(top_200_features.T, id_vars=['ecosystem_type'])     # melt data
top_200_features_melted[['0','position','AA','1']]= top_200_features_melted['variable'].str.split('(\d+)_([A-Za-z]+|-)', expand=True)   #split the feature into position and the AA
top_200_features_melted[['position']]=top_200_features_melted[['position']].apply(pd.to_numeric)        # covert the position to a numeric
top_200_features_melted = top_200_features_melted.drop(['0','1'],axis=1).sort_values("position")     # sort by position

top_200_features_melted
# adding shap or gini importance to dataframe
top_200_features_melted.reset_index(inplace=True, names='importance_order')
xx = top_200_features_melted.merge(shap_importance, left_on='variable', right_on='col_name')
top_200_features_melted = xx.drop('col_name', axis=1)
print(top_200_features_melted)

position_df = top_200_features_melted[['value', 'position', 'AA','ecosystem_type']]
print(position_df)

position_df_marine = position_df.loc[position_df['ecosystem_type'] == 'proportion_thermal']
print(position_df_marine)
position_df_new_marine = position_df_marine.pivot_table(index=["position"], columns='AA', values='value', fill_value=0.00001)
# position_df_new_marine.drop('-',inplace=True,axis=1)


position_df_freshwater = position_df.loc[position_df['ecosystem_type'] == 'proportion_freshwater']
position_df_new_freshwater = position_df_freshwater.pivot_table(index=["position"], columns='AA', values='value', fill_value=0.00001)
print(position_df_new_freshwater)
# position_df_new_freshwater.drop('-',inplace=True, axis=1)

position_df_new_marine.index = pd.to_numeric(position_df_new_marine.index)
position_df_new_freshwater.index = pd.to_numeric(position_df_new_freshwater.index)

# %%
# top_200_features_melted
top_200_features_melted_for_bars = top_200_features_melted.sort_values(
    by=['variable', 'value'], 
    ascending=[True, False], 
    key=lambda x: x if x.name == 'variable' else abs(x)  # Only apply abs to 'value' column
).groupby('variable').head(n=1)

top_200_features_melted_for_bars.loc[top_200_features_melted_for_bars['ecosystem_type'] == 'proportion_freshwater']
top_200_features_melted_for_bars

# %%
# Create the new figure and subplots
fig, axes = plt.subplots(4, 1, sharex=True, figsize=(300, 40))

# Create the marine logo directly on the first subplot
AA_logo_marine = logomaker.Logo(
    position_df_new_marine,
    ax=axes[1],  # Use the first subplot's axes
    font_name='Arial Rounded MT Bold',
    color_scheme='dmslogo_charge'
)

# Create the freshwater logo directly on the second subplot
AA_logo_fresh = logomaker.Logo(
    position_df_new_freshwater,
    ax=axes[2],  # Use the second subplot's axes
    font_name='Arial Rounded MT Bold',
    color_scheme='dmslogo_charge',
    flip_below=False
)

# create marine bars above first seqlogo (axes[0])
data_m=top_200_features_melted_for_bars.loc[top_200_features_melted_for_bars['ecosystem_type'] == 'proportion_thermal']
data_m["position"] = data_m["position"].astype(float)

axes[0].bar(
    data_m["position"], 
    data_m["feature_importance_vals"], 
    color="#ff2f00",
    width=0.75  # Adjust width to match sequence logo spacing
)

# create freshwater bars below last seqlogo (axes[3])
data_f = top_200_features_melted_for_bars.loc[top_200_features_melted_for_bars['ecosystem_type'] == 'proportion_freshwater']
data_f["position"] = data_f["position"].astype(float)

axes[3].bar(
    data_f["position"], 
    data_f["feature_importance_vals"], 
    color="#599900",
    width=0.75
)
axes[3].invert_yaxis()

for ax in axes:
    ax.set_xticks([0, 300])

# Ensure Logomaker uses the same ticks
AA_logo_marine.ax.set_xticks(position_df_new_freshwater.index)
AA_logo_fresh.ax.set_xticks(position_df_new_freshwater.index)

axes[-1].set_xticklabels(position_df_new_freshwater.index, rotation=90, fontsize=100)


axes[0].set_ylabel("Feature Importance Thermal springs",fontsize=50)
axes[1].set_ylabel("SeqLogo Thermal springs",fontsize=50)
axes[2].set_ylabel("SeqLogo Freshwater",fontsize=50)
axes[3].set_ylabel("Feature Importance Freshwater",fontsize=50)

axes[0].tick_params(labelsize=80)
axes[1].tick_params(labelsize=80)
axes[2].tick_params(labelsize=80)
axes[3].tick_params(labelsize=80)

# plt.yticks(fontsize=80)

# Adjust layout
plt.tight_layout()

# Display the combined figure
plt.show()



# %%
top_200_features_melted['value'] = abs(top_200_features_melted['value'])

# top_200_features_melted.loc[top_200_features_melted['value']]
temp = top_200_features_melted_for_bars.sort_values(['variable','value'], ascending=[True,False]).drop_duplicates(['variable'])
temp

# grouped = temp.groupby(['ecosystem_type','AA']).count().reset_index().drop(['value','position'],axis=1)
# grouped.rename(columns={'variable' : 'count_of_AA'}, inplace=True)
# grouped.sort_values(['ecosystem_type','count_of_AA'],ascending=[False, False], inplace=True)

# grouped

# %%
sns.set_theme(style="whitegrid")
my_palette = sns.diverging_palette(h_neg=11, h_pos=200,s=100,l=60, n=2)
min,max = 2, 50
plt.figure(figsize=(6,5))

b = sns.barplot(
    data=grouped,
    x="count_of_AA",
    y="AA",
    hue="ecosystem_type",
    palette=my_palette,
    dodge=True,
    width=-.5,  
)

# show the graph
# b.set_yticklabels(top_200_features_melted.variable)
b.set_xlabel("Amino acid",fontsize=15)
b.set_ylabel("count_of_amino acids in by biome",fontsize=15)
plt.show()


