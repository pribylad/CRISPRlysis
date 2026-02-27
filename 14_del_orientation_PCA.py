# pylint: disable=line-too-long
# pylint: disable=invalid-name

# %%
# Load table /home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/05_pivot_DEL_orientation_normalised_avg260225.tsv
import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
##########################################################################################
output_format = 'svg'
show_plots = True

# # 1 no filter, analysis on all samples after merging 02_mastertable*.tsv tables from CRISPRlysis outputs and generating 05_tables by 06_DEL.py script
# filter_criterion = 'None'
# filter_value = 'None'
# data_source = 'all_samples_no_filters'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/05_pivot_DEL_orientation_normalised_avgNone_None_260227.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples')

# 2 with data filtered for overall efficacy >=5%, input table is result from 02_filter_low_efficacy.py script
filter_criterion = 'efficacy'
filter_value = '95'
data_source = 'all_' + filter_criterion + '_' + filter_value + '%'
df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/05_pivot_DEL_orientation_normalised_avgefficacy_95_260227.tsv'
charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_{filter_criterion}_{filter_value}')
table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_{filter_criterion}_{filter_value}')

# # 3 with data filtered for overall efficacy >=5% and 1ntIN >=80%, input table is result from 13_IN_in_use.py script
# filter_criterion = 'efficacy_1ntIN'
# filter_value = '95_80'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_efficacy_1ntIN_95_80/05_pivot_DEL_orientation_normalised_avgefficacy_1ntIN_95_80_260227.tsv'
# data_source = 'all_' + filter_criterion + '_' + filter_value + '%'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/03_filtered_{filter_criterion}_{filter_value}')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_{filter_criterion}_{filter_value}')

##########################################################################################
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')
##########################################################################################
# Load the dataframe
df = pd.read_csv(df_dir, sep='\t')
# Remove the columne Unnamed: 0 if it exists
if 'Unnamed: 0' in df.columns:
    df = df.drop(columns=['Unnamed: 0'])    

# Make column combining organism and cell_line and hours_post_inoculation for better visualisation in the plot
df['organism_cell_line_hpi'] = df['organism'] + ' - ' + df['cell_line'] + ' - ' + df['hours_post_inoculation']
df
# %%
df['organism_cell_line_hpi'].unique()
# %%
# Make PCA analysis where vlues are from column values and categories are from column variable, hue is from column organism_cell_line_hpi and style is from column variable
# remove columns source, mutant, gene, gRNA_seq, ID_no_rep
df_pca = df.drop(columns=['source', 'mutant', 'gene', 'gRNA_seq', 'ID_no_rep', 'kingdom', 'organism', 'cell_line', 'hours_post_inoculation'])
# set organism_cell_line_hpi as index
df_pca = df_pca.set_index('organism_cell_line_hpi')
# make pivot table with columns as variable and values as values
df_pca = df_pca.pivot_table(index=df_pca.index, columns='variable', values='value')
df_pca
# %%
# Make dictionary for colors and markesrs for each index value. Plant - green, animal - pink, markes based on the time which is last part of the index 8 - circle, 24 - square, 48 - triangle
index_values = df_pca.index.unique()
properties_dict = {
    'Homo sapiens - HCT116 - 4': {'color': "#FF00F2", 'marker': '*'},
    'Homo sapiens - HCT116 - 8': {'color': "#FF00F2", 'marker': 'o'},
    'Homo sapiens - HCT116 - 16': {'color': "#FF00F2", 'marker': '^'},
    'Homo sapiens - HCT116 - 24': {'color': "#FF00F2", 'marker': '<'},
    'Homo sapiens - HCT116 - 48': {'color': "#FF00F2", 'marker': 's'},
    'Homo sapiens - HEK293 - 4': {'color': "#910089", 'marker': '*'},
    'Homo sapiens - HEK293 - 8': {'color': "#910089", 'marker': 'o'},
    'Homo sapiens - HEK293 - 16': {'color': "#910089", 'marker': '^'},
    'Homo sapiens - HEK293 - 24': {'color': "#910089", 'marker': '<'},
    'Homo sapiens - HEK293 - 48': {'color': "#910089", 'marker': 's'},
    'Homo sapiens - K562 - 4': {'color': '#520047', 'marker':'*'},
    'Homo sapiens - K562 - 8': {'color': '#520047', 'marker':'o'},
    'Homo sapiens - K562 - 16': {'color': '#520047', 'marker':'^'}, 
    'Homo sapiens - K562 - 24': {'color': "#520047", 'marker':'<'},
    'Homo sapiens - K562 - 48': {'color': '#520047', 'marker':'s'},
    'Homo sapiens - HAP1 - 48to72': {'color': 'pink', 'marker': 's'},
    'Homo sapiens - Hela - 72': {'color': 'brown', 'marker': 's'},
    'Homo sapiens - Jurkat - 48to72': {'color': 'purple', 'marker': 's'},
    'Arabidopsis thaliana - xxx - 504': {'color': 'green', 'marker': 'o'},
    'Arabidopsis thaliana - xxx - unknown': {'color': 'green', 'marker': '*'},
    'Nicotiana benthamiana - xxx - 120': {'color': 'green', 'marker': 'v'},
    'Petunia hybrida - xxx - 24': {'color': 'green', 'marker': '>'},
    'Silene latifoliaa - xxx - xx': {'color': 'green', 'marker': '<'},
    'Solanum lycopersicum - xxx - 24': {'color': 'green', 'marker': 's'}}

# Check if all index values in df_pca are in properties_dict
missing_indices = [index for index in df_pca.index if index not in properties_dict]
if missing_indices:
    raise ValueError(f"Missing indices in properties_dict: {missing_indices}")
# Remove from poperties_dict any indices that are not in df_pca
properties_dict = {key: value for key, value in properties_dict.items() if key in df_pca.index}
# %%
# Data are normalised, so we can directly apply PCA without standardisation
pca = PCA(n_components=2)
pca_result = pca.fit_transform(df_pca)
# Create a new dataframe with PCA results and index from df_pca
df_pca_result = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=df_pca.index)
df_pca_result
# %%
# Print the explained variance ratio for each principal component
explained_variance_ratio = pca.explained_variance_ratio_
print(f"Explained variance ratio: PC1: {explained_variance_ratio[0]:.2f}, PC2: {explained_variance_ratio[1]:.2f}")
# %%
# Plot PCA results with seaborn, hue is from index and style is from index
sns.set_theme(style='whitegrid')
plt.figure(figsize=(5, 5))
# plot wit properties_dict
for index in df_pca_result.index:
    color = properties_dict[index]['color']
    marker = properties_dict[index]['marker']
    plt.scatter(df_pca_result.loc[index, 'PC1'], df_pca_result.loc[index, 'PC2'], c=color, marker=marker, s=100)
plt.title('PCA of DEL orientation normalised average')
# plot x and y labels as PC1 and PC2 with explained variance ratio*100 in brackets
plt.xlabel(f"PC1 ({explained_variance_ratio[0]*100:.2f}%)")
plt.ylabel(f"PC2 ({explained_variance_ratio[1]*100:.2f}%)")
# plot legend with properties_dict order as key in the dictionary
legend_elements = []
for key in properties_dict:
    color = properties_dict[key]['color']
    marker = properties_dict[key]['marker']
    legend_elements.append(plt.Line2D([0], [0], marker=marker, color='w', label=key, markerfacecolor=color, markersize=10))
plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
# plot grid in the background
plt.grid(axis='both', alpha=0.5)
# Save PCA plot
output_file = os.path.join(charts_output_dir, f'06_PCA_DEL_orientation_normalised_avg_{date}.{output_format}')
plt.savefig(output_file, bbox_inches='tight', dpi=300)
# Save PCA results to TSV
output_file = os.path.join(table_output_dir, f'06_PCA_DEL_orientation_normalised_avg_{date}.tsv')
df_pca_result.to_csv(output_file, sep='\t')
# %%
