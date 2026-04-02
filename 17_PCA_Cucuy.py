# pylint: disable=line-too-long
# pylint: disable=invalid-name

# %%
# Load table /home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/04_split_index_Solanum lycopersicum - M82 - 48_260402.tsv


import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os
##########################################################################################
output_format = 'png'
# output_format = 'svg'
show_plots = True

# TODO, the dataset must be merged with the Cucuy data to have the columns "within_atac_peak", "leaf_exp", "t0_exp", and "description" by 16_merge_Cucuy_tables.py

# with data filtered for overall efficacy >=5%, input table is result from 02_filter_low_efficacy.py script
filter_criterion = 'Cucuy'
filter_value = 'all'
data_source = 'Cucuy_filtered_95'
df_dir = '/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/02_filtered_efficacy_95/04_Solanum_lycopersicum_M82_48_260402/05_pivot_DEL_orientation_normalised_avg_merged_260402.tsv'
charts_output_dir = Path(f'/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/04_charts/02_filtered_efficacy_95/04_Solanum_lycopersicum_M82_48_260402')
table_output_dir = Path(f'/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/02_filtered_efficacy_95/04_Solanum_lycopersicum_M82_48_260402')


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
# Convert df['hours_post_inoculation'] to string
df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)

# Make column combining organism and cell_line and hours_post_inoculation for better visualisation in the plot
# df['organism_cell_line_hpi'] = df['ProExoInt'] + ' - ' + df['cell_line'] + ' - ' + df['hours_post_inoculation']
df.columns
# %%
df.head()
# df['organism_cell_line_hpi'].unique()
# %%
# Make PCA analysis where vlues are from column values and categories are from column variable, hue is from column organism_cell_line_hpi and style is from column variable
# remove columns source, mutant, gene, ID_no_rep
df_pca = df.drop(columns=['source', 'organism', 'cell_line', 'mutant', 'gene','hours_post_inoculation', 'description'])
# Change values in within_atac_peak column to 1 if True and 0 if False
df_pca['within_atac_peak'] = df_pca['within_atac_peak'].apply(lambda x: 1 if x == True else 0)
# change values in leaf_exp and t0_exp to high=3, med=2, low=1
def convert_exp_to_numeric(x):
    if x == 'high':
        return 3
    elif x == 'med':
        return 2
    elif x == 'low':
        return 1
    else:
        return 0  # for any other value, including NaN
df_pca['leaf_exp'] = df_pca['leaf_exp'].apply(convert_exp_to_numeric)
df_pca['t0_exp'] = df_pca['t0_exp'].apply(convert_exp_to_numeric)
df_pca = df_pca.set_index('gRNA_seq')
df_pca
# %%
# Make dictionary with key as index of df_pca and values as within_atac_peak, leaf_exp, and t0_exp for later use in plotting
properties_dict = {}
for index, row in df_pca.iterrows():
    properties_dict[index] = {
        'within_atac_peak': row['within_atac_peak'],
        'leaf_exp': row['leaf_exp'],
        't0_exp': row['t0_exp']
    }
properties_dict
# %%
# make pivot table with columns as variable and values as values
df_pca = df_pca.pivot_table(index=df_pca.index, columns='variable', values='value')
df_pca
# %%
# Plot with clusters based on PCA results and properties_dict styling
# Data are normalised, so we can directly apply PCA without standardisation
pca = PCA(n_components=2)
pca_result = pca.fit_transform(df_pca)

# Create a new dataframe with PCA results and index from df_pca
df_pca_result = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=df_pca.index)

# Add K-means clustering on PCA results (e.g., 3 clusters)
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=3, random_state=42)
df_pca_result['Cluster'] = kmeans.fit_predict(pca_result)
# Add columns from properties_dict to df_pca_result
df_pca_result['within_atac_peak'] = df_pca_result.index.map(lambda x: properties_dict[x]['within_atac_peak'])
df_pca_result['leaf_exp'] = df_pca_result.index.map(lambda x: properties_dict[x]['leaf_exp'])
df_pca_result['t0_exp'] = df_pca_result.index.map(lambda x: properties_dict[x]['t0_exp'])
df_pca_result  # Now includes 'Cluster' column
# %%
# Print the explained variance ratio for each principal component
explained_variance_ratio = pca.explained_variance_ratio_
print(f"Explained variance ratio: PC1: {explained_variance_ratio[0]:.2f}, PC2: {explained_variance_ratio[1]:.2f}")

# %%
def plot_pca_clusters(df_pca, cluster_col, charts_output_dir, date=None):
    """
    Plot PCA scores with clusters colored and custom markers.

    Parameters
    ----------
    df_pca : pd.DataFrame
        Must have at least columns: 'PC1', 'PC2', and the cluster column.
    cluster_col : str
        Name of the column in df_pca containing cluster labels (e.g., 'Cluster').
    charts_output_dir : str
        Directory to save the plot.
    date : str, optional
        Suffix string for the filename (e.g., '2026_04_02').
    """
    sns.set_theme(style='whitegrid')
    plt.figure(figsize=(8, 6))

    # Auto‑generate cluster colors (up to 10 clusters)
    unique_clusters = sorted(df_pca[cluster_col].unique())
    color_list = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink',
                  'gray', 'olive', 'cyan']
    cluster_colors = dict(zip(unique_clusters, color_list[:len(unique_clusters)]))

    # Scatter points one by one (cluster color + marker from properties_dict)
    for index in df_pca.index:
        cluster = df_pca.loc[index, cluster_col]
        color = cluster_colors[cluster]
        plt.scatter(
            df_pca.loc[index, 'PC1'],
            df_pca.loc[index, 'PC2'],
            c=color,
            s=100,
            edgecolors='black',
            linewidth=0.5
        )

    # Labels and title
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA')

    plt.legend(*plt.gca().get_legend_handles_labels(), title={cluster_col})

    # Build filename
    date_str = date if date else "unknown"
    output_file = os.path.join(
        charts_output_dir,
        f'004_1_PCA_DEL_orientation_normalised_avg_with_{cluster_col}_{date_str}.png'
    )
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    return plt.show()

plot_pca_clusters(df_pca_result, 'Cluster', charts_output_dir, date)
plot_pca_clusters(df_pca_result, 'within_atac_peak', charts_output_dir, date)
plot_pca_clusters(df_pca_result, 'leaf_exp', charts_output_dir, date)
plot_pca_clusters(df_pca_result, 't0_exp', charts_output_dir, date)

# %%
# Plot PCA results with clusters colored + properties_dict styling
sns.set_theme(style='whitegrid')
plt.figure(figsize=(8, 6))

cluster_colors = {0: 'red', 1: 'blue', 2: 'green'}  # Define cluster colors
for index in df_pca_result.index:
    cluster = df_pca_result.loc[index, 'Cluster']
    color = cluster_colors[cluster]  # Cluster-based color
    plt.scatter(df_pca_result.loc[index, 'PC1'], df_pca_result.loc[index, 'PC2'], 
                c=color, s=100, edgecolors='black', linewidth=0.5)

plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA with K-means Clusters')
plt.legend(*plt.gca().get_legend_handles_labels(), title='Clusters')
# Save PCA plot with clusters
output_file = os.path.join(charts_output_dir, f'004_1_PCA_DEL_orientation_normalised_avg_with_clusters_{date}.png')
plt.savefig(output_file, bbox_inches='tight', dpi=300)
plt.show()
# %%
# Plot based on the samples
# Data are normalised, so we can directly apply PCA without standardisation
# pca = PCA(n_components=2)
# pca_result = pca.fit_transform(df_pca)
# # Create a new dataframe with PCA results and index from df_pca
# df_pca_result = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=df_pca.index)
# df_pca_result
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
output_file = os.path.join(charts_output_dir, f'004_2_PCA_DEL_orientation_normalised_avg_{date}.{output_format}')
plt.savefig(output_file, bbox_inches='tight', dpi=300)
# Save PCA results to TSV
output_file = os.path.join(table_output_dir, f'06_PCA_DEL_orientation_normalised_avg_{date}.tsv')
df_pca_result.to_csv(output_file, sep='\t')
# %%
