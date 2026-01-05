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
filter_criterion = 'cell_cycle'  # parameter used to filter the df_dir table by
# filter_value = 'dividing'  # value of the parameter used to filter the df_dir table by
filter_value = 'nondividing'  # value of the parameter used to filter the df_dir table by

df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/04_split_{filter_criterion}_{filter_value}_260105.tsv'
data_source = filter_criterion + '_' + filter_value  # 'dividing' or 'nondividing'
project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/{filter_criterion}_{filter_value}')
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
file_number = '02'

# Define list
default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'replicate']
list_of_events = ['1', '2', '3', '4', '5', '>6']
df_IN = df[default_list + ['n_inserted', '4th_nt', '4th_nt_inserted', '%Reads']]

#########################################################################################
#        Number of inserted nucleotides around the cut site
#########################################################################################
# Get data
pivot = pd.pivot_table(df_IN, index=default_list, columns='n_inserted', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
# remove column 0
df_pivot = df_pivot.drop(columns=[0], inplace=False)
# Make sum from the column len(default_list)+len(list_of_events) to the end as >6
df_pivot['>6'] = df_pivot.iloc[:, len(default_list)+len(list_of_events)-1 : ].sum(axis=1)
# Convert columns to string
df_pivot.columns = df_pivot.columns.map(str)
# Keep only default_list + list_of_events columns
df_pivot = df_pivot[default_list + list_of_events]
# Normalise rows to 100%
df_pivot[list_of_events] = df_pivot[list_of_events].div(df_pivot[list_of_events].sum(axis=1), axis=0) * 100
# Melt dataframe
df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
df_pivot.head()
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list[:-1], columns='variable', values='value', aggfunc='mean')
# Met pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list[:-1], value_vars=list_of_events)

counter = 2
color_dic = {'0':'#ffffff', '1':'#A8B7FF', '2':'#536CE8', '3':'#466270', '4':'#514C7A', '5':'#4E8387'\
                    , '>6':'#6F79A8'}
title='Normalised numbers of inserted nucleotides around the cut site, mean of replicates'
ylabel='Normalised percentage\n of inserted nucleotides'
no_of_samples = len(pivot_table_avg['gRNA_seq'].unique())
text_bellow_chart = f'Only IN events around cut.\
                    \nSum of all IN in one sample = 100%\
                    \nNo of samples: {no_of_samples}'
name_to_save=f'{file_number}_{counter}_IN_no_of_IN_nucleotides_summary_violinplot'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value', number_of_columns=6, show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, order=None, project_output_dir=project_output_dir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=title, xtickt_rotation=0, ylabel=ylabel, xlabel='No. of inserted nucleotides', yaxis_formater=100, name_to_save=name_to_save, text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=102, output_format=output_format, show_plots=show_plots)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
stat_info = pivot_table_avg.groupby('variable').value.describe().unstack()
save_stat(stat_info, name_to_save, project_output_dir)

##########################################################################################
#        Is the 4th nucleotide the inserted nucleotide?
##########################################################################################
list_4th_nt = default_list + ['4th_nt']
# Filter df_IN for 1nd insertion only
df_1IN = df_IN[df_IN['n_inserted'] == 1]
# Get data
pivot = pd.pivot_table(df_1IN, index=list_4th_nt, columns='4th_nt_inserted', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
# Normalise rows to 100%
df_pivot.iloc[:, len(list_4th_nt):] = df_pivot.iloc[:, len(list_4th_nt):].div(df_pivot.iloc[:, len(list_4th_nt):].sum(axis=1), axis=0) * 100
# Melt dataframe
df_pivot = pd.melt(df_pivot, id_vars=list_4th_nt, value_vars=[True, False])
# remove from list_4th_nt the 'replicate' element
list_4th_nt.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, list_4th_nt, columns='variable', values='value', aggfunc='mean')
# Met pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=list_4th_nt, value_vars=[True, False])
# Filter out False values of variable
pivot_table_avg = pivot_table_avg[pivot_table_avg['variable'] == True]
# Rename column variable to '4th_nt_inserted'
pivot_table_avg = pivot_table_avg.rename(columns={'variable': '4th_nt_inserted'})

counter = 6
# make hue_order according to present 4th nucleotides
nucleotide_order_list = {'A':1, 'C':2, 'G':3, 'T':4}
nucleotide_list = pivot_table_avg['4th_nt'].unique()
sort_list = []
for nucleotide in nucleotide_list:
    sort_list.append(nucleotide_order_list[nucleotide])
hue_order = [x for _, x in sorted(zip(sort_list, nucleotide_list))]     

title='Is the 4th nucleotide the inserted nucleotide? Mean of replicates'
no_of_samples = len(pivot_table_avg['gRNA_seq'].unique())
ylabel='Normalised percentage\n of inserted 4th nucleotides'
text_bellow_chart = f'Only 1nt IN around the cut are plotted.\
                    \nSequences without DEL/SUB around the cut.\
                    \nIn individual sample True + False = 100%.\
                    \nHorizontal lines show medians for individual 4th nucleotide.\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_IN_4th_nt_inserted_normalised_violinplot'

plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
    , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
    , show_hue_medians=False\
    , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
    , clip_right=True\
    , title=title, xtickt_rotation=0, ylabel=ylabel\
    , order=None, plot_legend=True \
    , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.05, yaxis_formater=100\
    , x_max=1\
    , name_to_save=name_to_save\
    , project_output_dir=project_output_dir \
    , data_source=data_source, output_format=output_format\
    , dpi=200\
    , y_min=0\
    , y_max=102, show_plots=show_plots
        )
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# Save statistics
stat_info = pivot_table_avg.groupby('4th_nt_inserted').describe().unstack()
save_stat(stat_info, name_to_save, project_output_dir)
# %%