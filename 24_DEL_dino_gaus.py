# pylint: disable=line-too-long
# pylint: disable=invalid-name

# This script is for analysis of DELs around the cut site, their orientation and ends participating on the repair.

# %%
# import necessary libraries
from datetime import datetime
from pathlib import Path
import matplotlib
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import seaborn as sns
import statannotations
from statannotations.Annotator import Annotator
import statsmodels

from scripts.plot_functions import plot_sns_violin_and_print, save_stat
from scripts.dict_functions import get_color_dic_35_violin
##########################################################################################
output_format = 'png'
# output_format = 'svg'
show_plots = False
individual_dataset = True    # If False (for settings 1, 2, and 3), Kingdom and ID_no_reps is used in the default_list and additional plots for them are made.
use_3_5_possible = False

# # 1 no filter, analysis on all samples after merging 02_mastertable*.tsv tables from CRISPRlysis outputs
# filter_criterion = 'None'
# filter_value = 'None'
# data_source = 'all_samples_no_filters'
# individual_dataset = False
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_All_260409.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_all_samples')

# 2 with data filtered for overall efficacy >=5%, input table is result from 02_filter_low_efficacy.py script
filter_criterion = 'efficacy'
filter_value = '95'
individual_dataset = False
data_source = 'all_samples_filtered_95'
df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/01_all_samples')
table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/01_all_samples')
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
# data_source = 'all_samples_filtered_95'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/01_all_samples')


##########################################################################################
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')

if charts_output_dir.exists() is False:
    os.makedirs(charts_output_dir)

if table_output_dir.exists() is False:
    os.makedirs(table_output_dir)
##########################################################################################
def create_pivot_table(df, list_of_indexes, columns, values, aggfunc):
    pivot_table = pd.pivot_table(
        df,
        index=list_of_indexes,
        columns=columns,
        values=values,
        aggfunc=aggfunc,
        fill_value=0
    )
    return pivot_table
##########################################################################################
# Load datafreame
df = pd.read_csv(df_dir, sep='\t')
df.columns
# %%
df.shape
# %%
df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)
df['blunt_cut_with_MH_possible'] = df['blunt_cut_with_MH_possible'].astype(str)
# change TRUE to True and FALSE to False in blunt_cut_with_MH_possible column
df['blunt_cut_with_MH_possible'] = df['blunt_cut_with_MH_possible'].replace({'TRUE':'True', 'FALSE':'False'})
print(f"Loaded dataframe with shape {df.shape} from {df_dir}")
# %%
df['blunt_cut_with_MH_possible'].unique()
# %%
df['n_deleted_around_cut'].unique()
# %%
##########################################################################################
#               DDDD    EEEEEE  L
#               D   D   E       L
#               D   D   EEEE    L
#               D   D   E       L
#               DDDD    EEEEEE  LLLLLL
##########################################################################################

# Define list of indexes for pivot table
# default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism', 'replicate']
if individual_dataset == False:
    default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep','n_deleted_around_cut','gRNA_num', 'replicate']
else:
    default_list = ['source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'n_deleted_around_cut','gRNA_num', 'replicate']
list_of_events = ['upstream_ambiguous', 'both_ambiguous', 'downstream_ambiguous', 'upstream_NO_ambi', 'both_NO_ambi', 'downstream_NO_ambi', 'cannot_be_determined_ambiguous', 'cannot_be_determined_NO_ambi', 'X_gap_out_of_cut']
df_DEL_orientation = df[default_list + ['DEL_orientation', 'Gap_around_cut', 'Maximal_MH_len', '%Reads']]
# Filter out rows where Gap_around_cut is False
df_DEL_orientation = df_DEL_orientation[df_DEL_orientation['Gap_around_cut'] == True]
# Replace values in the DEL_orientation column
df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'].replace({
    'U': 'upstream',
    'B': 'both',
    'D': 'downstream',
    'C': 'cannot_be_determined'
})
order = ['upstream_ambiguous', 'both_ambiguous', 'downstream_ambiguous', 'upstream_NO_ambi', 'both_NO_ambi', 'downstream_NO_ambi', 'cannot_be_determined_ambiguous', 'cannot_be_determined_NO_ambi', 'X_gap_out_of_cut']
print(df_DEL_orientation['Maximal_MH_len'].unique())
# Replace in Maximal_MH_len '0' for 'NO_ambi', 'gap_out_of_cut' keep and anything else as 'ambiguous'
df_DEL_orientation['Maximal_MH_len'] = df_DEL_orientation['Maximal_MH_len'].replace({'0': 'NO_ambi', 'gap_out_of_cut': 'gap_out_of_cut'})
df_DEL_orientation['Maximal_MH_len'] = df_DEL_orientation['Maximal_MH_len'].apply(lambda x: 'ambiguous' if x not in ['NO_ambi', 'gap_out_of_cut'] else x)
# Replace column 'DEL_orientation' with merged column of DEL_orientation and Maximal_MH_len
df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'] + '_' + df_DEL_orientation['Maximal_MH_len']
print(df_DEL_orientation['DEL_orientation'].unique())
print(df['n_deleted_around_cut'].unique())

df_DEL_orientation
# %%
# Remove columns Gap_around_cut and Maximal_MH_len as they are not needed anymore
df_DEL_orientation = df_DEL_orientation.drop(columns=['Gap_around_cut', 'Maximal_MH_len'], inplace=False)
df_DEL_orientation
# %%
##########################################################################################
#        Deletion orientation around the cut site
##########################################################################################
# Get data
pivot = pd.pivot_table(df_DEL_orientation, index=default_list, columns='DEL_orientation', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
df_pivot
# %%
# remove index name
df_pivot.columns.name = None
# Save pivot table as tsv
df_pivot.to_csv(Path(str(table_output_dir) + '/' + '09' + '_DEL_orientation_NOTnormalised_pivot_DElsize' + '.tsv'), sep='\t')
# Normalise rows to 100%
df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
df_pivot
# %%
# Melt dataframe
# Check if all list_of_events items are in columns, if not create them with 0 values
for item in list_of_events:
    if item not in df_pivot.columns.tolist():
        df_pivot[item] = 0
df_pivot = df_pivot.melt(id_vars=default_list, value_vars=list_of_events)
# Remove 'replicate' from default_list
default_list.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# Melt pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(table_output_dir) + '/' + '09_pivot_DEL_orientation_normalised_avg' + '.tsv'), sep='\t')
# Make column combining organism and cell_line for better visualisation in the plot
pivot_table_avg['organism_cell_line_hpi'] = pivot_table_avg['organism'] + ' - ' + pivot_table_avg['cell_line'] + ' - ' + pivot_table_avg['hours_post_inoculation']
pivot_table_avg
# %%
# make n_deleted_around_cut as integer 
pivot_table_avg['n_deleted_around_cut'] = pivot_table_avg['n_deleted_around_cut'].astype(int)
pivot_table_avg['n_deleted_around_cut'].unique()
# %%
# Normalise value column to 100%
pivot_table_avg['value'] = pivot_table_avg['value'].div(pivot_table_avg['value'].sum()) * 100
# %%
# Keep only rows with variable in ['upstream_ambiguous', 'upstream_NO_ambi', 'downstream_ambiguous', 'downstream_NO_ambi']
pivot_table_avg = pivot_table_avg[pivot_table_avg['variable'].isin(['upstream_ambiguous', 'upstream_NO_ambi', 'downstream_ambiguous', 'downstream_NO_ambi'])]

# Set n_deleted_around_cut negative for upstream and positive for downstream
pivot_table_avg.loc[pivot_table_avg['variable'].str.contains('upstream'), 'n_deleted_around_cut'] = -pivot_table_avg['n_deleted_around_cut']
pivot_table_avg
# %%

# %%
pivot_table_avg['n_deleted_around_cut'].unique()
# %%

# %%
# Make table keeping only n_deleted_around_cut, variable and value columns and make average of value
table_for_plot = create_pivot_table(pivot_table_avg, ['n_deleted_around_cut'], columns='variable', values='value', aggfunc='sum')
table_for_plot
# %%
# save table_for_plot as tsv
table_for_plot.to_csv(Path(str(table_output_dir) + '/' + '09_table_for_plot_gaussian_DEL' + '.tsv'), sep='\t')
# %%
# Make stacked ber plot with n_deleted_around_cut on x axis and values NO_ambi as they are and on top of them values from ambiguous
plt.figure(figsize=(10, 6))
sns.barplot(data=pivot_table_avg, x='n_deleted_around_cut', y='value', hue='variable', estimator='sum', errorbar=None)
plt.xlabel('Number of bases deleted around cut site')
plt.ylabel('Percentage of reads')
plt.xlim(50, 150)
plt.title('Distribution of deletion orientation around cut site')
plt.legend(title='DEL orientation')
plt.tight_layout()