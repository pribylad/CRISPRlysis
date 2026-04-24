# pylint: disable=line-too-long
# pylint: disable=invalid-name

# %%
import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import umap
# %%
output_format = 'png'
show_plots = True

# #  1 no filter, analysis on all samples after merging 02_mastertable*.tsv tables from CRISPRlysis outputs and generating 05_tables by 06_DEL.py script
# filter_criterion = 'None'
# filter_value = 'None'
# data_source = 'all_samples_no_filters'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_all_samples/05_pivot_DEL_orientation_normalised_avg.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_all_samples')

# %%
# filter_criterion = 'efficacy'
# filter_value = '95'
# data_source = 'all_samples_filtered_95'
# # read the DEL_orientation data where the biological replicates are already averaged and the data is filtered for efficacy > 95%
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/01_all_samples/05_pivot_DEL_orientation_normalised_avg.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_{filter_criterion}_{filter_value}/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_{filter_criterion}_{filter_value}/01_all_samples')


# 4 with data from Crisitna Ponce Lilly
filter_criterion = 'efficacy'
filter_value = '95'
data_source = 'all_samples_filtered_95_with_Ponce'
df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/01_all_samples/05_pivot_DEL_orientation_normalised_avg.tsv'
charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_{filter_criterion}_{filter_value}/with_Ponce_samples/01_all_samples')
table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_{filter_criterion}_{filter_value}/with_Ponce_samples/01_all_samples')
# %%
df = pd.read_csv(df_dir, sep='\t', index_col=0)
df.head()
df['variable'].unique()
df.columns
# %%
# Make column '4th_nt' from the 'gRNA_seq' column by taking the 4th character from the right end of the gRNA sequence (which is the 17th character from the left)
df['4th_nt'] = df['gRNA_seq'].str[17] 
df['4th_3rd_nt'] = df['gRNA_seq'].str[17:19]
# %%
metadata_cols = [
    "kingdom", "source", "organism", "cell_line", "mutant",
    "gene", "gRNA_seq", "4th_nt", "4th_3rd_nt", "hours_post_inoculation", "Org_line_hpi"
]

feature_col = "variable"
value_col = "value"
# Replace spaces and . and , in column source with underscores
df['source'] = df['source'].str.replace('.,', '_').str.replace(',', '_').str.replace(' ', '_').str.replace('__', '_')
print(df['source'].unique())
df['organism'] = df['organism'].str.replace('A. thaliana', 'Arabidopsis thaliana').str.replace(' ', '_')
print(df['organism'].unique())
df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str).str.replace('xxx', 'unknown').str.replace('xx', 'unknown').str.replace('unknown', 'xxx')
print(df['hours_post_inoculation'].unique())
df['Org_line_hpi'] = df['organism'] + '_' + df['cell_line'] + '_' + df['hours_post_inoculation']
print(df['Org_line_hpi'].unique())
# %%
df['4th_nt'].unique()
df['4th_3rd_nt'].unique()

# %%
# 1) Create a sample ID from metadata columns
df = df.copy()
df["sample_id"] = df[metadata_cols].astype(str).agg("_".join, axis=1)
df['sample_id'].unique()
# %%
# 2) Pivot long -> wide
wide = df.pivot_table(
    index="sample_id",
    columns=feature_col,
    values=value_col,
    aggfunc="mean"   # use mean in case duplicates exist
)
wide
# %%
# Find duplicate sample_ids
duplicate_sample_ids = wide.index[wide.index.duplicated()].unique()
print(f"Duplicate sample IDs: {duplicate_sample_ids}")
# %%
# 3) Keep one metadata row per sample
meta = df.drop_duplicates("sample_id")[["sample_id"] + metadata_cols].set_index("sample_id")
meta
# %%
# 4) Align metadata with the wide matrix
wide = wide.loc[meta.index]
wide
# %%
# 5) Fill missing values if needed
wide = wide.fillna(0)
# Save the wide matrix and metadata for future use
wide.to_csv(table_output_dir / f'07_UMAP_wide_matrix_{data_source}.tsv', sep='\t')
# %%
# 6) Run UMAP
reducer = umap.UMAP(
    n_neighbors=15,
    min_dist=0.1,
    metric="euclidean",
    random_state=42
)

embedding = reducer.fit_transform(wide)
# %%
# 7) Store results
out = meta.copy()
out["UMAP1"] = embedding[:, 0]
out["UMAP2"] = embedding[:, 1]

# %%
properties_dic = {'stagg2nt':{'color':'#E890B8','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'stagg2nt_likely':{'color':'#E890B8','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'stagg1nt':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'stagg1nt_likely':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'////', 'marker':'o'}
            ,'blunt':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'blunt_likely':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'cannot_be_determined':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'can_not_be_determined':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'can not be determined':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'cannot_be_determined_ambiguous':{'color':'white','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'cannot_be_determined_NO_ambi':{'color':'white','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'distal_stagg2nt_likely':{'color':'#E890B8','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'distal_stagg2nt':{'color':'#E890B8','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'distal_stagg1nt_likely':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'distal_stagg1nt':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'distal_blunt_likely':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'\\\\', 'marker':'o'}
            ,'distal_blunt':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_blunt':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_blunt_likely':{'color':'#6ea7ff','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'proximal_stagg1nt':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_stagg1nt_likely':{'color':'#ffaaaa','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'proximal_stagg2nt':{'color':'#E890B8','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'proximal_stagg2nt_likely':{'color':'#E890B8','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'gap_out_of_cut':{'color':'gray','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'IN':{'color':'#091f92','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'INs':{'color':'#091f92','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'DEL':{'color':'#ba8759','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'DELs':{'color':'#ba8759','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'SUB':{'color':'#ffcc33','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'OTHERS':{'color':'#808080','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'no_mutation':{'color':'#111111','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'U':{'color':'#f0d271','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'B':{'color':'#808080','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'D':{'color':'#f5ad62','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'X':{'color':'#d3d3d3','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'X_gap_out_of_cut':{'color':'#d3d3d3','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'upstream':{'color':'#f0d271','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'upstream_ambiguous':{'color':'#f0d271','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'upstream_NO_ambi':{'color':'#f0d271','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'upstream_both':{'color':'#C4C4C4','edge_color':'#f0d271', 'hatch':'\\\\\\', 'marker':'o'}
            ,'upstream_both_downstream':{'color':'#979797','edge_color':'black', 'hatch':'xx', 'marker':'o'}
            ,'both':{'color':'#808080','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'both_ambiguous':{'color':'#808080','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'both_NO_ambi':{'color':'#808080','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'both_downstream':{'color':'#C4C4C4','edge_color':'#f5ad62', 'hatch':'////', 'marker':'o'}
            ,'downstream':{'color':'#f5ad62','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'downstream_ambiguous':{'color':'#f5ad62','edge_color':'gray', 'hatch':'//', 'marker':'o'}
            ,'downstream_NO_ambi':{'color':'#f5ad62','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'deletion_spans_analysed_window':{'color':'#a9a9a9','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'deletion spans analysed window':{'color':'#a9a9a9','edge_color':'black', 'hatch':'', 'marker':'o'}
            ,'PAM_distal_nt_in_MH':{'color':'#52D8FF','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'PAM_proximal_nt_in_MH':{'color':'#FADA52','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'no_MH_likely':{'color':'#A9A9A9','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'no_MH':{'color':'#D3D3D3','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'A':{'color':'#70c936','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'C':{'color':'#ff0000','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'G':{'color':'#ffd42a','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'T':{'color':'#9955ff','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,True:{'color':'#D3D3D3','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,False:{'color':'#B8BEDB','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'True':{'color':'#D3D3D3','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'False':{'color':'#B8BEDB','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'unmethylated':{'color':'#000000','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'CG_only':{'color':'#4D4D4D','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'heterochromatin':{'color':'#999999','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'RdDM':{'color':'#fafafa','edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'Homo_sapiens_HCT116_4': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Homo_sapiens_HCT116_8': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Homo_sapiens_HCT116_16': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': '^'}
            ,'Homo_sapiens_HCT116_24': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': '<'}
            ,'Homo_sapiens_HCT116_48': {'color': "#E8B3D1", 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo_sapiens_HEK293_4': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Homo_sapiens_HEK293_8': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Homo_sapiens_HEK293_16': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': '^'}
            ,'Homo_sapiens_HEK293_24': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': '<'}
            ,'Homo_sapiens_HEK293_48': {'color': "#CC79A7", 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo_sapiens_K562_4': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'*'}
            ,'Homo_sapiens_K562_8': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'o'}
            ,'Homo_sapiens_K562_16': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'^'}
            ,'Homo_sapiens_K562_24': {'color': "#A55D84", 'edge_color':'gray', 'hatch':'', 'marker':'<'}
            ,'Homo_sapiens_K562_48': {'color': '#A55D84', 'edge_color':'gray', 'hatch':'', 'marker':'s'}
            ,'Homo_sapiens_HAP1_48to72': {'color': '#7D3E60', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo_sapiens_Jurkat_48to72': {'color': '#810F72', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Homo_sapiens_Hela_72': {'color': '#520047', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Arabidopsis_thaliana_xxx_504': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Arabidopsis_thaliana_root culture_504': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Arabidopsis_thaliana_xxx_xxx': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Arabidopsis_thaliana_embryo_xxx': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '*'}
            ,'Nicotiana_benthamiana_xxx_120': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'v'}
            ,'Nicotiana_benthamiana_xxx_336': {'color': "#11775c", 'edge_color':'gray', 'hatch':'', 'marker': 'v'}
            ,'Petunia_hybrida_xxx_24': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '>'}
            ,'Silene_latifolia_xxx_xxx': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '<'}
            ,'Solanum_lycopersicum_xxx_24': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            ,'Solanum_lycopersicum_M82_48': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': '^'}
            ,'plant': {'color': '#009d73', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'animal': {'color': '#cc79a6', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'UniCorn_random_dataset_0;': {'color': '#FFD700', 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Xenopus_tropicalis_xxx_xxx': {'color': "#181766", 'edge_color':'gray', 'hatch':'', 'marker': 'o'}
            ,'Pleurodeles_waltl_xxx_xxx': {'color': "#272592", 'edge_color':'gray', 'hatch':'', 'marker': 's'}
            }
# %%

def plot_umap(out, groupby, color_map, show_plots=False):
    plt.figure(figsize=(8, 8))
    sns.set_style("whitegrid")
    for group, g in out.reset_index().groupby(groupby):
        try:
            plt.scatter(
                g["UMAP1"],
                g["UMAP2"],
                s=35,
                alpha=0.5,
                color=color_map[group]['color'],
                marker=color_map[group]['marker'],
                label=group
            )
        except KeyError:
            print(f"Group {group} not found in color map")
            plt.scatter(
                g["UMAP1"],
                g["UMAP2"],
                s=35,
                alpha=0.3,
                label=group
            )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.xlim(-4,15)
    plt.ylim(1,15)
    plt.title("UMAP of deletion profiles")
    plt.legend(title=groupby, bbox_to_anchor=(1.05, 1), loc='upper left')
    # Save the plot
    plt.savefig(charts_output_dir / f'005_umap_{groupby}_{data_source}.{output_format}', dpi=300, bbox_inches='tight')
    if show_plots:
        plt.show()

plot_umap(out, "kingdom", properties_dic, show_plots=show_plots)
# %%
plot_umap(out, "4th_nt", properties_dic, show_plots=show_plots)
# %%
plot_umap(out, "3rd_4th_nt", properties_dic, show_plots=show_plots)
# %%
plot_umap(out, "Org_line_hpi", properties_dic, show_plots=show_plots)
# %%
# Identify clusters in the UMAP embedding using HDBSCAN
import hdbscan
cluster_size = 12
clusterer = hdbscan.HDBSCAN(min_cluster_size=cluster_size)
out['cluster'] = clusterer.fit_predict(embedding)
# Plot UMAP colored by cluster
plt.figure(figsize=(8, 7))
sns.scatterplot(
    x="UMAP1",
    y="UMAP2",
    hue="cluster",
    palette="tab10",
    data=out,
    legend="full",
    alpha=0.5
)
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.xlim(-4,15)
plt.ylim(1,15)
plt.title("UMAP of deletion profiles colored by HDBSCAN clusters")
plt.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc='upper left')
# Save the plot
plt.savefig(charts_output_dir / f'005_umap_clusters_{data_source}_{cluster_size}.{output_format}', dpi=300, bbox_inches='tight')
# Save the out dataframe with UMAP coordinates and cluster labels
out.to_csv(table_output_dir / f'07_umap_coordinates_clusters_{data_source}_{cluster_size}.tsv', sep='\t', index=True)
if show_plots:
    plt.show()
# %%
# Make UMAP plot highlighting the samples base on Org_line_hpi, rest samples in gray
for org_line_hpi, group in out.groupby("Org_line_hpi"):
    plt.figure(figsize=(8, 6))
    # Plot all samples in gray
    plt.scatter(
        out["UMAP1"],
        out["UMAP2"],
        s=35,
        alpha=0.3,
        color='lightgray',
        label='Other samples'
    )
    # Highlight the current group
    plt.scatter(
        group["UMAP1"],
        group["UMAP2"],
        s=35,
        alpha=0.5,
        color=properties_dic.get(org_line_hpi, {'color': 'red'})['color'],
        marker=properties_dic.get(org_line_hpi, {'marker': 'o'})['marker'],
        label=org_line_hpi
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.xlim(-4,15)
    plt.ylim(1,15)
    plt.title(f"UMAP of deletion profiles highlighting {org_line_hpi}")
    plt.legend(title="Org_line_hpi", bbox_to_anchor=(1.05, 1), loc='upper left')
    # Save the plot
    plt.savefig(charts_output_dir / f'005_umap_highlight_Org_line_hpi_{org_line_hpi}_{data_source}.{output_format}', dpi=300, bbox_inches='tight')
    if show_plots:
        plt.show()
# %%
# Make UMAP plot highlighting the samples base on 4th_3rd_nt, rest samples in gray
for nt, group in out.groupby("4th_3rd_nt"):
    plt.figure(figsize=(8, 6))
    # Plot all samples in gray
    plt.scatter(
        out["UMAP1"],
        out["UMAP2"],
        s=35,
        alpha=0.3,
        color='lightgray',
        label='Other samples'
    )
    # Highlight the current group
    plt.scatter(
        group["UMAP1"],
        group["UMAP2"],
        s=35,
        alpha=0.5,
        color=properties_dic.get(nt, {'color': 'red'})['color'],
        marker=properties_dic.get(nt, {'marker': 'o'})['marker'],
        label=nt
    )
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.xlim(-4,15)
    plt.ylim(1,15)
    plt.title(f"UMAP of deletion profiles highlighting 4th_3rd_nt={nt}")
    plt.legend(title="4th_3rd_nt", bbox_to_anchor=(1.05, 1), loc='upper left')
    # Save the plot
    plt.savefig(charts_output_dir / f'005_umap_highlight_4th_3rd_nt_{nt}_{data_source}.{output_format}', dpi=300, bbox_inches='tight')
    if show_plots:
        plt.show()

# %%
# split dataframe out into separate dataframes for each cluster
clusters = out['cluster'].unique()
cluster_dfs = {cluster: out[out['cluster'] == cluster] for cluster in clusters}
# Get the samples from the wide dataframe that belong to each cluster
cluster_samples = {cluster: cluster_df.index.tolist() for cluster, cluster_df in cluster_dfs.items()}
# Save the sample IDs for each cluster to separate files
for cluster, samples in cluster_samples.items():
    with open(table_output_dir / f'08_cluster_{cluster}_samples_{data_source}.txt', 'w') as f:
        for sample in samples:
            f.write(f"{sample}\n")

# %%
from scripts.plot_functions import plot_sns_violin_and_print, save_stat
from scripts.dict_functions import get_color_dic_35_violin

# filter the 'wide' dataframe to keep only the samples that belong to individual clusters and make individual plots for each cluster
for cluster in clusters:
    print(f"Processing cluster {cluster}...")
    cluster_samples = cluster_dfs[cluster].index.tolist()
    cluster_wide = wide.loc[cluster_samples]
    # chnge cluster_wide to long format for seaborn
    cluster_wide = cluster_wide.reset_index().melt(id_vars='sample_id', var_name='variable', value_name='value')
    # plot violinplots for cluster_wide

    # plot violinplots for cluster_0_wide
    colour_dic = {'upstream_ambiguous':'#f0d271', 'upstream_NO_ambi':'#f0d271', 'both_ambiguous':'#808080', 'both_NO_ambi':'#808080', 'downstream_ambiguous':'#cccccc', 'downstream_NO_ambi':'#cccccc', 'cannot_be_determined_ambiguous':'#ffffff', 'cannot_be_determined_NO_ambi':'#ffffff', 'X_gap_out_of_cut':'#000000', 'plant':'#009d73', 'animal':'#cc79a6'}
    list_of_events = ['upstream_ambiguous', 'both_ambiguous', 'downstream_ambiguous', 'upstream_NO_ambi', 'both_NO_ambi', 'downstream_NO_ambi', 'cannot_be_determined_ambiguous', 'cannot_be_determined_NO_ambi', 'X_gap_out_of_cut']
    ylabel = 'Normalised percentage \nof DEL around the cut'
    title = 'Normalised incidence of DEL according to the orientation,\n mean of replicates'
    name_to_save = f'006_1_DEL_orientation_cluster_{cluster}_violinplot_small'
    data_source = f'Cluster_{cluster}'
    plot = plot_sns_violin_and_print(data=cluster_wide, x='variable', y='value'\
        , number_of_columns=len(cluster_wide['variable'].unique())/2\
        , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
        , show_hue_medians=False\
        , hue=None, hue_only_for_data_points=False \
        , hatches=False\
        , clip_right=True\
        , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
        , order=list_of_events\
        , plot_legend=False\
        , integrated_properties_dic=False, color_dic=colour_dic\
        , text_bellow_chart=None\
        , y_legend_coordinate=-0.4, yaxis_formater=100\
        , x_max=None\
        , change_edge_color=False\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir\
        , data_source=data_source, output_format='png'\
        , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# %%
