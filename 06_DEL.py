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
from scripts.dict_functions import get_color_dic_35_violin
##########################################################################################
output_format = 'png'
show_plots = True

filter_criterion = 'All'
filter_value = 'All'
df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_efficacy_samples_limit5efficacy_260109.tsv'
data_source = 'all_samples'
project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/all_main_samples/{filter_criterion}_{filter_value}')
#####
# filter_criterion = 'cell_cycle'  # parameter used to filter the df_dir table by
# filter_value = 'dividing'  # value of the parameter used to filter the df_dir table by
# filter_value = 'nondividing'  # value of the parameter used to filter the df_dir table by
# df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/04_split_{filter_criterion}_{filter_value}_260109.tsv'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/nondividing_animal/04_split_cell_cycle_nondividing_260110.tsv'
# data_source = filter_criterion + '_' + filter_value  # 'dividing' or 'nondividing'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/{filter_criterion}_{filter_value}')

#####
# filter_criterion = 'Event_around_cut'
# filter_value = 'IN'
# df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/04_split_Event_around_cut_IN_above_treshold_260106.tsv'
# data_source = filter_criterion + '_' + filter_value + '_above_treshold'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/{filter_criterion}_{filter_value}_above_treshold')
# df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/04_split_Event_around_cut_IN_below_treshold_260106.tsv'
# data_source = filter_criterion + '_' + filter_value + '_below_treshold'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/{filter_criterion}_{filter_value}_below_treshold')
#####
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
#               DDDD    EEEEEE  L
#               D   D   E       L
#               D   D   EEEE    L
#               D   D   E       L
#               DDDD    EEEEEE  LLLLLL
##########################################################################################
file_number = '03'

# Define list of indexes for pivot table
default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism','replicate']
list_of_events = ['upstream', 'both', 'downstream', 'cannot_be_determined', 'X']
df_DEL_orientation = df[default_list + ['DEL_orientation', 'Gap_around_cut', '%Reads']]
# Filter out rows where Gap_around_cut is False
df_DEL_orientation = df_DEL_orientation[df_DEL_orientation['Gap_around_cut'] == True]
# Replace values in the DEL_orientation column
df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'].replace({
    'U': 'upstream',
    'B': 'both',
    'D': 'downstream',
    'C': 'cannot_be_determined'
})
order = ['upstream', 'both', 'downstream', 'cannot_be_determined']
df_DEL_orientation

##########################################################################################
#        Deletion orientation around the cut site
##########################################################################################
# Get data
pivot = pd.pivot_table(df_DEL_orientation, index=default_list, columns='DEL_orientation', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
print(df_pivot.head())
# Save pivot table as tsv
df_pivot.to_csv(Path(str(project_output_dir) + '/' + file_number + '_DEL_orientation_NOTnormalised_pivot_' + date + '.tsv'), sep='\t')
# Normalise rows to 100%
df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
# Melt dataframe
df_pivot = df_pivot.melt(id_vars=default_list, value_vars=list_of_events)
# Remove 'replicate' from default_list
default_list.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# Met pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
pivot_table_avg


counter = 2
# Set charts properties
title = 'Normalised incidence of DEL according to the orientation,\n mean of replicates'
ylabel = 'Normalised percentage \nof DEL around the cut'
no_of_samples =  len(pivot_table_avg[['source', 'gRNA_seq', 'organism']].drop_duplicates())
text_bellow_chart = f'Only DELs around the cut are plotted. Without gaps spanning analysed window\
                    \n - see table with source data.\
                    \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_2_DEL_orientation_summary_violinplot'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
    , number_of_columns=len(pivot_table_avg['variable'].unique())/2\
    , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
    , show_hue_medians=False\
    , hue=None, hue_only_for_data_points=False \
    , hatches=True\
    , clip_right=True\
    , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
    , order=order\
    , plot_legend=False\
    , integrated_properties_dic=True, color_dic=None\
    , text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.4, yaxis_formater=100\
    , x_max=None\
    , change_edge_color=True\
    , name_to_save=name_to_save\
    , project_output_dir=project_output_dir\
    , data_source=data_source, output_format=output_format\
    , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
stat_info = pivot_table_avg.groupby('variable').value.describe()
save_stat(stat_info, name_to_save, project_output_dir)

##################################################################################
#       Ends participating on the repair
##################################################################################
file_number = '05'
counter = 3

# Define list of indexes for pivot table
default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'replicate']
list_of_events = ['3\' NOT blunt', '5\' possible blunt', '5\' NOT blunt', 'eider 3\' or 5\' possible blunt', 'eider 3\' or 5\' NOT blunt']
df_end = df[default_list + ['Gap_around_cut', 'Maximal_MH_len', 'End_3_5', 'blunt_cut_with_MH_possible', '%Reads']]
# Remove rows where Maximal_MH_len is 'no_gap', '0' or 'gap_out_of_cut'
df_end = df_end[~df_end['Maximal_MH_len'].isin(['no_gap', '0', 'gap_out_of_cut'])]
# filter out deletions which are not around the cut
df_end = df_end[df_end['Gap_around_cut'] == True]
# Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'eider 3\' or 5\''})
# Change True and False in  blunt_cut_with_MH_possible to 'possible blunt' and 'NOT possible blunt'
df_end['blunt_cut_with_MH_possible'] = df_end['blunt_cut_with_MH_possible'].replace({'TRUE':'possible blunt', 'FALSE':'NOT blunt'})

# Make column where both End_3_5 and blunt_cut_with_MH_possible are combined
df_end['End_3_5_blunt'] = df_end['End_3_5'] + ' ' + df_end['blunt_cut_with_MH_possible']
# remove columns End_3_5, blunt_cut_with_MH_possible
df_end = df_end.drop(columns=['End_3_5', 'blunt_cut_with_MH_possible'], inplace=False)

# Get data
pivot = pd.pivot_table(df_end, index=default_list, columns='End_3_5_blunt', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
# Normalise rows to 100%
df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
df_pivot.head()

# Melt dataframe
df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
# remove 'replicate' from default_list
default_list.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# Melt pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
pivot_table_avg

#########################################################################################
# plot VIOLINPLOT
##########################################################################################
# NORMALISED
# Set charts properties for violinplot
title = 'Normalised ends participating on the repair, mean of replicates'
ylabel = 'Normalised percentage of ends \n around cut events'
text_bellow_chart = f'Only DELs around the cut are plotted.\
                    \nNormalised to: SUM of all ends in individual sample = 100%\
                    \n3\' possible blunt can not be in -2,0,0 setup -> it is eider 3\' or 5\' NOT blunt\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_violinplot'
color_dic = get_color_dic_35_violin()
# TODO split it to True/False for can be blunt cut with MH
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
    , number_of_columns=3\
    , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
    , show_hue_medians=False\
    , hue=None, hue_only_for_data_points=False \
    , hatches=False\
    , clip_right=True\
    , title=title, xtickt_rotation=-55, ylabel=ylabel\
    , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']\
    , plot_legend=False\
    , integrated_properties_dic=False, color_dic=color_dic\
    , text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.45, yaxis_formater=100\
    , x_max=None\
    , change_edge_color=False\
    , name_to_save=name_to_save\
    , project_output_dir=project_output_dir\
    , data_source=data_source, output_format=output_format\
    , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# Save statistics
stat_info = pivot_table_avg.groupby('variable').value.describe()
save_stat(stat_info, name_to_save, project_output_dir)

# %%
