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
from statannotations.Annotator import Annotator

from scripts.plot_functions import plot_sns_violin_and_print, save_stat
from scripts.dict_functions import get_color_dic_35_violin
##########################################################################################
output_format = 'png'
show_plots = True
x = 95  # efficacy threshold percentage

# # 1 no filter, analysis on all samples after merging 02_mastertable*.tsv tables from CRISPRlysis outputs
# filter_criterion = 'None'
# filter_value = 'None'
# data_source = 'all_samples_no_filters'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_merged_02tables_NoFilter_All_260225.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples')

# # 2 with data filtered for overall efficacy >=5%, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy'
# filter_value = '95'
# data_source = 'all_' + filter_criterion + '_' + filter_value + '%'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All_260226.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_{filter_criterion}_{filter_value}')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_{filter_criterion}_{filter_value}')

# 3 with data filtered for overall efficacy >=5% and 1ntIN >=80%, input table is result from 13_IN_in_use.py script
filter_criterion = 'efficacy_1ntIN'
filter_value = '95_80'
df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/04_filtered_efficacy_95_samples_260226.tsv'
data_source = 'all_' + filter_criterion + '_' + filter_value + '%'
charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/03_filtered_{filter_criterion}_{filter_value}')
table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_{filter_criterion}_{filter_value}')

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
df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)
print(f"Loaded dataframe with shape {df.shape} from {df_dir}")
##########################################################################################
#               DDDD    EEEEEE  L
#               D   D   E       L
#               D   D   EEEE    L
#               D   D   E       L
#               DDDD    EEEEEE  LLLLLL
##########################################################################################
# Plot full figure for supplementary material
file_number = '001_Supp'

# Define list of indexes for pivot table
# default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism', 'replicate']
default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep', 'replicate']
list_of_events = ['upstream', 'both', 'downstream', 'cannot be determined', 'deletion spans analysed window']
df_DEL_orientation = df[default_list + ['DEL_orientation', 'Gap_around_cut', '%Reads']]
# Filter out rows where Gap_around_cut is False
df_DEL_orientation = df_DEL_orientation[df_DEL_orientation['Gap_around_cut'] == True]
# Replace values in the DEL_orientation column
df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'].replace({
    'U': 'upstream',
    'B': 'both',
    'D': 'downstream',
    'C': 'cannot be determined',
    'X': 'deletion spans analysed window'
})
order = ['upstream', 'both', 'downstream', 'cannot be determined', 'deletion spans analysed window']

properties_dic = {'upstream':{'color':'#f0d271','edge_color':'gray', 'hatch':''}\
            ,'both':{'color':'#808080','edge_color':'black', 'hatch':''}\
            ,'downstream':{'color':'#f5ad62','edge_color':'gray', 'hatch':''}\
            ,'cannot be determined':{'color':'#d3d3d3','edge_color':'black', 'hatch':'//'}\
            ,'deletion spans analysed window':{'color':'#a9a9a9','edge_color':'black', 'hatch':'xx'}}
# Make color dic from properties_dic
color_dic = {key: properties_dic[key]['color'] for key in properties_dic.keys()}
df_DEL_orientation
# %%
##########################################################################################
#        Deletion orientation around the cut site
##########################################################################################
# Get data
pivot = pd.pivot_table(df_DEL_orientation, index=default_list, columns='DEL_orientation', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
# Save pivot table as tsv
df_pivot.to_csv(Path(str(table_output_dir) + '/' + file_number + '_DEL_orientation_NOTnormalised_pivot_' + filter_criterion + '_' + filter_value + '_' + date + '.tsv'), sep='\t')
# Normalise rows to 100%
df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
# If all list_of_events items are not in columns, create them with 0 values
for item in list_of_events:
    if item not in df_pivot.columns.tolist():
        df_pivot[item] = int()
# Melt dataframe
df_pivot = df_pivot.melt(id_vars=default_list, value_vars=list_of_events)
# Remove 'replicate' from default_list
default_list.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# Met pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
pivot_table_avg.to_csv(Path(str(table_output_dir) + '/' + file_number + '_DEL_orientation_normalised_pivot_' + filter_criterion + '_' + filter_value + '_' + date + '.tsv'), sep='\t')

# %%
counter = 1
# Set charts properties
title = 'Normalised incidence of DEL according to the orientation'
ylabel = 'Normalised percentage \nof DEL around the cut'
no_of_samples =  len(pivot_table_avg[['source', 'gRNA_seq', 'organism']].drop_duplicates())
text_bellow_chart = f'Only DELs around the cut are plotted.\
                    \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                    \nAny replicates are averaged prior plotting.\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_2_DEL_orientation_summary_violinplot_{filter_criterion}_{filter_value}'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
    , number_of_columns=len(pivot_table_avg['variable'].unique())/2\
    , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
    , show_hue_medians=False\
    , hue=None, hue_only_for_data_points=False \
    , hatches=False\
    , clip_right=True\
    , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
    , order=order\
    , plot_legend=False\
    , integrated_properties_dic=properties_dic, color_dic=None\
    , text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.55, yaxis_formater=100\
    , x_max=None\
    , change_edge_color=False\
    , name_to_save=name_to_save\
    , project_output_dir=charts_output_dir\
    , data_source=data_source, output_format=output_format\
    , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
stat_info = pivot_table_avg.groupby('variable').value.describe()
save_stat(stat_info, name_to_save, charts_output_dir)
# %%
# Plot main figure without 'deletion spans analysed window'
file_number = '001_Fig'

# Define list of indexes for pivot table
# default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism','replicate']
default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep', 'replicate']
list_of_events = ['upstream', 'both', 'downstream', 'cannot be determined']
df_DEL_orientation = df[default_list + ['DEL_orientation', 'Gap_around_cut', '%Reads']]
# Filter out rows where Gap_around_cut is False
df_DEL_orientation = df_DEL_orientation[df_DEL_orientation['Gap_around_cut'] == True]
# Replace values in the DEL_orientation column
df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'].replace({
    'U': 'upstream',
    'B': 'both',
    'D': 'downstream',
    'C': 'cannot be determined',
    'X': 'deletion spans analysed window'
})
order = ['upstream', 'both', 'downstream', 'cannot be determined']

properties_dic = {'upstream':{'color':'#f0d271','edge_color':'gray', 'hatch':''}\
            ,'both':{'color':'#808080','edge_color':'black', 'hatch':''}\
            ,'downstream':{'color':'#f5ad62','edge_color':'gray', 'hatch':''}\
            ,'cannot be determined':{'color':'#d3d3d3','edge_color':'black', 'hatch':'//'}\
            ,'deletion spans analysed window':{'color':'#a9a9a9','edge_color':'black', 'hatch':'xx'}}
# Make color dic from properties_dic
color_dic = {key: properties_dic[key]['color'] for key in properties_dic.keys()}
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
# Save pivot table as tsv
df_pivot.to_csv(Path(str(table_output_dir) + '/' + file_number + '_DEL_orientation_NOTnormalised_pivot_' + filter_criterion + '_' + filter_value + '_' + date + '.tsv'), sep='\t')
# Normalise rows to 100%
df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
# If all list_of_events items are not in columns, create them with 0 values
for item in list_of_events:
    if item not in df_pivot.columns.tolist():
        df_pivot[item] = int()
# Melt dataframe
df_pivot = df_pivot.melt(id_vars=default_list, value_vars=list_of_events)
# Remove 'replicate' from default_list
default_list.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# Met pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)

counter = 1
# Set charts properties
title = 'Normalised incidence of DEL according to the orientation'
ylabel = 'Normalised percentage \nof DEL around the cut'
no_of_samples =  len(pivot_table_avg[['source', 'gRNA_seq', 'organism']].drop_duplicates())
text_bellow_chart = f'Only DELs around the cut are plotted. Without "deletion spans analysed window" category - is in supplementary material.\
                    \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                    \nAny replicates are averaged prior plotting.\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_2_DEL_orientation_summary_violinplot_{filter_criterion}_{filter_value}'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
    , number_of_columns=len(pivot_table_avg['variable'].unique())/2\
    , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
    , show_hue_medians=False\
    , hue=None, hue_only_for_data_points=False \
    , hatches=False\
    , clip_right=True\
    , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
    , order=order\
    , plot_legend=False\
    , integrated_properties_dic=properties_dic, color_dic=None\
    , text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.55, yaxis_formater=100\
    , x_max=None\
    , change_edge_color=False\
    , name_to_save=name_to_save\
    , project_output_dir=charts_output_dir\
    , data_source=data_source, output_format=output_format\
    , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
stat_info = pivot_table_avg.groupby('variable').value.describe()
save_stat(stat_info, name_to_save, charts_output_dir)
# %%
# ##################################################################################
# #       Ends participating on the repair
# ##################################################################################
# # %%
# file_number = '05'
# counter = 3

# # Define list of indexes for pivot table
# default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'replicate']
# list_of_events = ['3\' NOT blunt', '5\' possible blunt', '5\' NOT blunt', 'eider 3\' or 5\' possible blunt', 'eider 3\' or 5\' NOT blunt']
# df_end = df[default_list + ['Gap_around_cut', 'Maximal_MH_len', 'End_3_5', 'blunt_cut_with_MH_possible', '%Reads']]
# # Remove rows where Maximal_MH_len is 'no_gap', '0' or 'gap_out_of_cut'
# df_end = df_end[~df_end['Maximal_MH_len'].isin(['no_gap', '0', 'gap_out_of_cut'])]
# # filter out deletions which are not around the cut
# df_end = df_end[df_end['Gap_around_cut'] == True]
# # Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
# df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'eider 3\' or 5\''})
# # Change True and False in  blunt_cut_with_MH_possible to 'possible blunt' and 'NOT possible blunt'
# df_end['blunt_cut_with_MH_possible'] = df_end['blunt_cut_with_MH_possible'].replace({'TRUE':'possible blunt', 'FALSE':'NOT blunt'})

# # Make column where both End_3_5 and blunt_cut_with_MH_possible are combined
# df_end['End_3_5_blunt'] = df_end['End_3_5'] + ' ' + df_end['blunt_cut_with_MH_possible']
# # remove columns End_3_5, blunt_cut_with_MH_possible
# df_end = df_end.drop(columns=['End_3_5', 'blunt_cut_with_MH_possible'], inplace=False)

# # Get data
# pivot = pd.pivot_table(df_end, index=default_list, columns='End_3_5_blunt', values='%Reads', aggfunc='sum', fill_value=0)
# # reset pivot and remove index
# df_pivot = pivot.reset_index()
# # remove index name
# df_pivot.columns.name = None
# # Normalise rows to 100%
# df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
# df_pivot.head()

# # Melt dataframe
# df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
# # remove 'replicate' from default_list
# default_list.remove('replicate')
# # Make pivot table with all columns without column replicate and make average of values
# pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# # Melt pivot table_avg reset index
# pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
# pivot_table_avg

# #########################################################################################
# # plot VIOLINPLOT
# ##########################################################################################
# # NORMALISED
# # Set charts properties for violinplot
# title = 'Normalised ends participating on the repair, mean of replicates'
# ylabel = 'Normalised percentage of ends \n around cut events'
# text_bellow_chart = f'Only DELs around the cut are plotted.\
#                     \nNormalised to: SUM of all ends in individual sample = 100%\
#                     \n3\' possible blunt can not be in -2,0,0 setup -> it is eider 3\' or 5\' NOT blunt\
#                     \nNo of samples: {no_of_samples}'
# name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_violinplot'
# color_dic = get_color_dic_35_violin()
# # TODO split it to True/False for can be blunt cut with MH
# plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
#     , number_of_columns=3\
#     , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
#     , show_hue_medians=False\
#     , hue=None, hue_only_for_data_points=False \
#     , hatches=False\
#     , clip_right=True\
#     , title=title, xtickt_rotation=-55, ylabel=ylabel\
#     , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']\
#     , plot_legend=False\
#     , integrated_properties_dic=False, color_dic=color_dic\
#     , text_bellow_chart=text_bellow_chart\
#     , y_legend_coordinate=-0.45, yaxis_formater=100\
#     , x_max=None\
#     , change_edge_color=False\
#     , name_to_save=name_to_save\
#     , project_output_dir=project_output_dir\
#     , data_source=data_source, output_format=output_format\
#     , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# # Save pivot table as tsv
# pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# # Save statistics
# stat_info = pivot_table_avg.groupby('variable').value.describe()
# save_stat(stat_info, name_to_save, project_output_dir)
# # %%
# # Plot only 5' NOT blunt and 3' NOT blunt from the pivot_table_avg as x variables
# name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_violinplot_filtered'
# color_dic = get_color_dic_35_violin()
# import statannotations
# print(statannotations.__version__)
# plot = plot_sns_violin_and_print(data=pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])], x='variable', y='value'\
#     , number_of_columns=2\
#     , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
#     , show_hue_medians=False\
#     , hue=None, hue_only_for_data_points=False \
#     , hatches=False\
#     , clip_right=True\
#     , title=title, xtickt_rotation=-55, ylabel=ylabel\
#     , order=['5\' NOT blunt', '3\' NOT blunt']\
#     , plot_legend=False\
#     , integrated_properties_dic=False, color_dic=color_dic\
#     , text_bellow_chart=text_bellow_chart\
#     , y_legend_coordinate=-0.25, yaxis_formater=100\
#     , x_max=None\
#     , change_edge_color=False\
#     , add_stat=True, stat_pairs=[("5' NOT blunt", "3' NOT blunt")], stat_test='t-test_ind', stat_format='simple', stat_correction=None\
#     , name_to_save=name_to_save\
#     , project_output_dir=project_output_dir\
#     , data_source=data_source, output_format=output_format\
#     , dpi=200, y_min=0, y_max=119, show_plots=show_plots)


# # %%
# # Plot End participating on the repair with respect of the 4th nucleotid
# counter = 4
# # Define list of indexes for pivot table
# default_list = ['source', 'gene', 'gRNA_seq', '4th_nt', 'mutant', 'replicate']
# list_of_events = ['3\' NOT blunt', '5\' possible blunt', '5\' NOT blunt', 'eider 3\' or 5\' possible blunt', 'eider 3\' or 5\' NOT blunt']
# df_end = df[default_list + ['Gap_around_cut', 'Maximal_MH_len', 'End_3_5', 'blunt_cut_with_MH_possible', '%Reads']]
# # Remove rows where Maximal_MH_len is 'no_gap', '0' or 'gap_out_of_cut'
# df_end = df_end[~df_end['Maximal_MH_len'].isin(['no_gap', '0', 'gap_out_of_cut'])]
# # filter out deletions which are not around the cut
# df_end = df_end[df_end['Gap_around_cut'] == True]
# # Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
# df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'eider 3\' or 5\''})
# # Change True and False in  blunt_cut_with_MH_possible to 'possible blunt' and 'NOT possible blunt'
# df_end['blunt_cut_with_MH_possible'] = df_end['blunt_cut_with_MH_possible'].replace({'TRUE':'possible blunt', 'FALSE':'NOT blunt'})

# # Make column where both End_3_5 and blunt_cut_with_MH_possible are combined
# df_end['End_3_5_blunt'] = df_end['End_3_5'] + ' ' + df_end['blunt_cut_with_MH_possible']
# # remove columns End_3_5, blunt_cut_with_MH_possible
# df_end = df_end.drop(columns=['End_3_5', 'blunt_cut_with_MH_possible'], inplace=False)

# # Get data
# pivot = pd.pivot_table(df_end, index=default_list, columns='End_3_5_blunt', values='%Reads', aggfunc='sum', fill_value=0)
# # reset pivot and remove index
# df_pivot = pivot.reset_index()
# # remove index name
# df_pivot.columns.name = None
# # Normalise rows to 100%
# df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
# df_pivot.head()

# # Melt dataframe
# df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
# # remove 'replicate' from default_list
# default_list.remove('replicate')
# # Make pivot table with all columns without column replicate and make average of values
# pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# # Melt pivot table_avg reset index
# pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
# pivot_table_avg
# # %%
# #########################################################################################
# # plot VIOLINPLOT
# ##########################################################################################
# # NORMALISED 4th as hue
# # Set charts properties for violinplot
# title = 'Normalised ends participating on the repair, mean of replicates'
# ylabel = 'Normalised percentage of ends \n around cut events'
# text_bellow_chart = f'Only DELs around the cut are plotted.\
#                     \nNormalised to: SUM of all ends in individual sample = 100%\
#                     \n3\' possible blunt can not be in -2,0,0 setup -> it is eider 3\' or 5\' NOT blunt\
#                     \nNo of samples: {no_of_samples}'
# name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_4thnt_violinplot'
# color_dic = get_color_dic_35_violin()
# # TODO split it to True/False for can be blunt cut with MH
# plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
#     , number_of_columns=15\
#     , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=True\
#     , show_hue_medians=False\
#     , boxplot_distances_with = [0.8, 0.18]\
#     , hue='4th_nt', hue_only_for_data_points=False \
#     , hatches=False\
#     , clip_right=True\
#     , title=title, xtickt_rotation=-55, ylabel=ylabel\
#     , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']\
#     , hue_order=['A', 'C','G','T']\
#     , plot_legend=True\
#     , integrated_properties_dic=True, color_dic=None\
#     , text_bellow_chart=text_bellow_chart\
#     , y_legend_coordinate=-0.45, yaxis_formater=100\
#     , x_max=None\
#     , change_edge_color=False\
#     , name_to_save=name_to_save\
#     , project_output_dir=project_output_dir\
#     , data_source=data_source, output_format=output_format\
#     , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# # Save pivot table as tsv
# pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# # Save statistics
# stat_info = pivot_table_avg.groupby('variable').value.describe()
# save_stat(stat_info, name_to_save, project_output_dir)

# # %%

# %%
