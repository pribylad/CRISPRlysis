# pylint: disable=line-too-long
# pylint: disable=invalid-name
# %%
# import necessary libraries
from datetime import datetime
from pathlib import Path
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from statannot import add_stat_annotation

from scripts.plot_functions import plot_sns_violin_and_print, save_stat
##########################################################################################
output_format = 'png'
show_plots = True

# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit5efficacy_All_260223.tsv'
df_dir = Path('/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/' + '04_filtered_efficacy_and_IN_samples' + '_' + date + '.tsv')
data_source = 'All samples with efficacy > 5%, average of replicates.'
project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/00_prefinal_cahrts')
##########################################################################################
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')

if project_output_dir.exists() is False:
    os.makedirs(project_output_dir)
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
print(f"Loaded dataframe with shape {df.shape} from {df_dir}")
##########################################################################################
#            I   N    N
#            I   N N  N
#            I   N  N N
#            I   N   NN
#            I   N    N
##########################################################################################
file_number = '002_Fig'
maximum_to_view = 2

# Define list
default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'name', 'replicate']
# Make a list of events with numbers from 1 to maximum_to_view and >maximum_to_view
list_of_events = [str(i) for i in range(1, maximum_to_view+1)]
list_of_events.append(f'>{maximum_to_view}')

df_IN = df[default_list + ['n_inserted', '4th_nt', '4th_nt_inserted', '%Reads']]
# Replace NaN values in cell line with 'None'
df_IN['cell_line'] = df_IN['cell_line'].fillna('None')
#########################################################################################
#        Number of inserted nucleotides around the cut site
#########################################################################################
# Get data
pivot = pd.pivot_table(df_IN, index=default_list, columns='n_inserted', values='%Reads', aggfunc='sum', fill_value=0)
pivot
# %%
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
# remove column 0
df_pivot = df_pivot.drop(columns=[0], inplace=False)
# Make sum from the column len(default_list)+len(list_of_events) to the end as >{maximum_to_view}
df_pivot[f'>{maximum_to_view}'] = df_pivot.iloc[:, len(default_list)+len(list_of_events)-1 : ].sum(axis=1)
# Convert columns to string
df_pivot.columns = df_pivot.columns.map(str)
# If all list_of_events items are not in columns, create them with 0 values
for item in list_of_events:
    if item not in df_pivot.columns.tolist():
        df_pivot[item] = 0
# Keep only default_list + list_of_events columns
df_pivot = df_pivot[default_list + list_of_events]
# Normalise rows to 100%
df_pivot[list_of_events] = df_pivot[list_of_events].div(df_pivot[list_of_events].sum(axis=1), axis=0) * 100
# Melt dataframe
df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list[:-1], columns='variable', values='value', aggfunc='mean')
# Melt pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list[:-1], value_vars=list_of_events)
print(f'Shape of the pivot table with average values: {pivot_table_avg.shape}')
#########################################################################################
# %%
panel = 'A'
color_dic = {'0':'#ffffff', '1':'#A8B7FF', '2':'#536CE8', '3':'#466270', '4':'#514C7A', '5':'#4E8387', f'>{maximum_to_view}':'#AAAAAA'}
title='Normalised numbers of inserted nucleotides around the cut site, mean of replicates'
ylabel='Normalised percentage\n of inserted nucleotides'
no_of_samples = len(pivot_table_avg) / len(pivot_table_avg['variable'].unique())
text_bellow_chart = f'Only IN events around cut.\
                    \nSum of all IN in one sample = 100%\
                    \nNo of samples: {no_of_samples}'
name_to_save=f'{file_number}_{panel}_IN_no_of_IN_nucleotides_summary_violinplot'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value', number_of_columns=len(pivot_table_avg['variable'].unique()), show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, order=None, project_output_dir=project_output_dir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=title, xtickt_rotation=0, ylabel=ylabel, xlabel='No. of inserted nucleotides', yaxis_formater=100, name_to_save=name_to_save, text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=102, output_format=output_format, show_plots=show_plots)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
stat_info = pivot_table_avg.groupby('variable').value.describe().unstack()
save_stat(stat_info, name_to_save, project_output_dir)
# %%
pivot_table_avg
# %%
# Make list of 'name' where Normalised percentage of inserted nucleotides for 1 is higher than 80% and filter the original dataframe fot these names and save it as 04_filtered_efficacy_and_IN_samples.tsv
names_with_1IN_high = pivot_table_avg[(pivot_table_avg['variable'] == '1') & (pivot_table_avg['value'] > 80)]['name'].tolist()
names_with_1IN_high
# %%
df_filtered = df[df['name'].isin(names_with_1IN_high)]
df_filtered.to_csv(Path('/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/' + '04_filtered_efficacy_and_IN_samples' + '_' + date + '.tsv'), sep='\t', index=False)
# %%
df_filtered
# %%
##########################################################################################
#        Is the 4th nucleotide the inserted nucleotide?
##########################################################################################
# list_4th_nt = default_list + ['4th_nt']
# # Filter df_IN for 1nd insertion only
# df_1IN = df_IN[df_IN['n_inserted'] == 1]
# df_1IN
# # %%
# # Get data
# pivot = pd.pivot_table(df_1IN, index=list_4th_nt, columns='4th_nt_inserted', values='%Reads', aggfunc='sum', fill_value=0)
# # reset pivot and remove index
# df_pivot = pivot.reset_index()
# # remove index name
# df_pivot.columns.name = None
# # Normalise rows to 100%
# df_pivot.iloc[:, len(list_4th_nt):] = df_pivot.iloc[:, len(list_4th_nt):].div(df_pivot.iloc[:, len(list_4th_nt):].sum(axis=1), axis=0) * 100
# # If all items are not in columns, create them with 0 values
# for item in [True, False]:
#     if item not in df_pivot.columns.tolist():
#         df_pivot[item] = int()
# # Melt dataframe
# df_pivot = pd.melt(df_pivot, id_vars=list_4th_nt, value_vars=[True, False])
# # remove from list_4th_nt the 'replicate' element
# list_4th_nt.remove('replicate')
# # Make pivot table with all columns without column replicate and make average of values
# pivot_table_avg = create_pivot_table(df_pivot, list_4th_nt, columns='variable', values='value', aggfunc='mean')
# # Met pivot table_avg reset index
# pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=list_4th_nt, value_vars=[True, False])
# # Filter out False values of variable
# pivot_table_avg = pivot_table_avg[pivot_table_avg['variable'] == True]
# # Rename column variable to '4th_nt_inserted'
# pivot_table_avg = pivot_table_avg.rename(columns={'variable': '4th_nt_inserted'})

# panel = 3
# # make hue_order according to present 4th nucleotides
# nucleotide_order_list = {'A':1, 'C':2, 'G':3, 'T':4}
# nucleotide_list = pivot_table_avg['4th_nt'].unique()
# sort_list = []
# for nucleotide in nucleotide_list:
#     sort_list.append(nucleotide_order_list[nucleotide])
# hue_order = [x for _, x in sorted(zip(sort_list, nucleotide_list))]     

# title='Is the 4th nucleotide the inserted nucleotide? Mean of replicates'
# no_of_samples = len(pivot_table_avg) / len(pivot_table_avg['4th_nt_inserted'].unique())
# ylabel='Normalised percentage\n of inserted 4th nucleotides'
# text_bellow_chart = f'Only 1nt IN around the cut are plotted.\
#                     \nSequences without DEL/SUB around the cut.\
#                     \nIn individual sample True + False = 100%.\
#                     \nHorizontal lines show medians for individual 4th nucleotide.\
#                     \nNo of samples: {no_of_samples}'
# name_to_save = f'{file_number}_{panel}_IN_4th_nt_inserted_normalised_violinplot'

# plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
#     , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
#     , show_hue_medians=False\
#     , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
#     , clip_right=True\
#     , title=title, xtickt_rotation=0, ylabel=ylabel\
#     , order=None, plot_legend=True \
#     , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
#     , y_legend_coordinate=-0.05, yaxis_formater=100\
#     , boxplot_distances_with=[0.8, 0.3]\
#     , x_max=1\
#     , name_to_save=name_to_save\
#     , project_output_dir=project_output_dir \
#     , data_source=data_source, output_format=output_format\
#     , dpi=200\
#     , y_min=0\
#     , y_max=102, show_plots=show_plots
#         )
# # Save pivot table as tsv
# pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# # Save statistics
# stat_info = pivot_table_avg.groupby('4th_nt_inserted').describe().unstack()
# save_stat(stat_info, name_to_save, project_output_dir)
# %%