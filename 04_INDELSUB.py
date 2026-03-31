# pylint: disable=line-too-long
# pylint: disable=invalid-name
# Description: This script takes the  dataframe from 01_merge_filter_02_tables.py and creates a pivot table with the percentage of each event (IN, DEL, SUB, OTHERS) for each sample. Then it creates a violin plot to show the distribution of each event across all samples.
# The script also checks if the sum of all events and no_mutation is 100% for each sample, if not it raises a ValueError and prints the rows where it is not the case.
# The script also checks if all unique values in the 'organism' column of the merged dataframe are present in the organism_dict, if not it raises a ValueError and prints the missing organism.
# The script also adds a column 'OTHERS' to the pivot table which is the sum of all events that are not IN, DEL or SUB (INDEL, INSUB, INDELSUB, DELSUB).
# The script also creates a pivot table with the average of each event for each sample, without the replicate column, and saves it as a tsv file.
# The script also creates a violin plot for the INDEL category only, by removing the OTHERS and SUB rows from the pivot table, and saves it as a png file.
# The script also saves the statistics of the pivot table as a tsv file.

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


#####
# filter_criterion = 'Kingdom'
# filter_value = 'Animals'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit5efficacy_Animals_260126.tsv'
# data_source = 'Animals'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/all_main_samples/{filter_criterion}_{filter_value}')
# #####
# filter_criterion = 'Kingdom'
# filter_value = 'Plants'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit5efficacy_Plants_260126.tsv'
# data_source = 'Plants'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/all_main_samples/{filter_criterion}_{filter_value}')


#####
# Analysis for data without filtering for 5% efficacy
# filter_criterion = 'All'
# filter_value = 'All'
# # df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All_260330.tsv'
# df_dir = '/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_All_260330.tsv'
# data_source = 'all_samples'
# project_output_dir = Path(f'/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/04_charts/01_no_filter_all_samples')

#####
# Analysis for data with filtering for 5% efficacy
filter_criterion = 'Efficacy'
filter_value = '95'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All_260330.tsv'
df_dir = '/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All_260330.tsv'
data_source = 'all_samples_filtered_95'
project_output_dir = Path(f'/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/04_charts/02_filtered_efficacy_95/01_all_samples')

# %%
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
df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)
print(f"Loaded dataframe with shape {df.shape} from {df_dir}")
df['hours_post_inoculation'].unique()
# %%
##########################################################################################
#            I   N    N   DDDD    EEEEEE  L        SSSS   U   U   BBBB
#            I   N N  N   D   D   E       L       S       U   U   B   B
#            I   N  N N   D   D   EEEE    L        SSS    U   U   BBBB
#            I   N   NN   D   D   E       L           S   U   U   B   B
#            I   N    N   DDDD    EEEEEE  LLLLLL  SSSS     UUU    BBBB
##########################################################################################
# Define lists
default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep', 'name', 'replicate']
list_of_events = ['IN', 'DEL', 'SUB', 'OTHERS']

default_list_NOreplicate = default_list.copy()
default_list_NOreplicate.remove('replicate')
default_list_NOreplicate.append('Event_around_cut')

# Get data
df_INDELSUB = df[default_list + ['Event_around_cut', '%Reads']]
list_of_others = ['INDEL', 'INSUB', 'INDELSUB', 'DELSUB']

pivot = pd.pivot_table(df_INDELSUB, index=default_list, columns='Event_around_cut', values='%Reads', aggfunc='sum', fill_value=0)
# Save pivot table as tsv into same folder as file from df_dir came from when df_dir is sting to a file
pivot.to_csv(Path(df_dir.rsplit('/', 1)[0] + '/' + f'02_INDELSUB_pivot_table_{date}.tsv'), sep='\t')
pivot
# %%
# Check if the sum of all events and no_mutation is 100% for each sample, if not print exit with VALUEERROR and print the rows where it is not the case
pivot['total'] = pivot.sum(axis=1)
if (pivot['total'] >= 100).any():
    print("Error: The sum of all events and no_mutation is not 100% for the following rows:")
    print(pivot['total'].max())
    raise ValueError("The sum of all events and no_mutation is not 100% for some samples. Please check the data.")

else:
    print("The sum of all events and no_mutation is 100% for all samples.")
# %%
df_INDELSUB_summ = pivot.reset_index()

# check if all columns exists, if not make an empty one
events_list = ['IN', 'DEL', 'SUB']
for item in events_list:
    if item not in df_INDELSUB_summ.columns.tolist():
        df_INDELSUB_summ[item] = int()   

## make column 'OTHERS' depending if the individual event is present or not
df_INDELSUB_summ['OTHERS'] = int()
for event in list_of_others:
    if event in df_INDELSUB_summ.columns.tolist():
        df_INDELSUB_summ['OTHERS'] = df_INDELSUB_summ['OTHERS'] + df_INDELSUB_summ[event]
# %%

df_INDELSUB_summ = df_INDELSUB_summ[default_list + list_of_events]
df_INDELSUB_summ
# %%
df_INDELSUB_summ = pd.melt(df_INDELSUB_summ, id_vars=default_list, value_vars=list_of_events)
# MERGE REPLICATES BY TAKING THE MEAN OF THE VALUES
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_INDELSUB_summ, default_list_NOreplicate, columns=None, values='value', aggfunc='mean')
pivot_table_avg.reset_index(inplace=True)
print(pivot_table_avg.head())
print(pivot_table_avg['source'].unique())
# %%
# get no of samples by taking no of rows in pivot_table_avg and dividing by unique values in column 'Event_around_cut'
no_of_samples = len(pivot_table_avg) / pivot_table_avg['Event_around_cut'].nunique()
print(f"No of samples: {no_of_samples}")
# %%
##########################################################################################
# plot VIOLINPLOT
##########################################################################################
file_number = '001'
counter = 1
# Set charts properties
title = 'INDELSUB representation in all samples, mean of replicates'
ylabel = 'Percentage'
text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                    \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                    \nIN+DEL+SUB+OTHERS+no_mutation = 100%.\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_INDELSUB_summary_violinplot'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='Event_around_cut', y='value'\
        , number_of_columns=3\
        , project_output_dir=project_output_dir\
        , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False, show_hue_medians=False\
        , hue='Event_around_cut', hue_only_for_data_points=True, transparency=0.7 \
        , clip_right=True, change_edge_color=False\
        , title=title, xtickt_rotation=0, ylabel=ylabel\
        , order=['IN', 'DEL', 'SUB', 'OTHERS'] \
        , integrated_properties_dic=True\
        , text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05\
        , yaxis_formater=100, name_to_save=name_to_save\
        , data_source=data_source, output_format=output_format\
        , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
stat_info = pivot_table_avg.groupby('Event_around_cut').value.describe().unstack()
save_stat(stat_info, name_to_save, project_output_dir)
# %%
# remove from the pivot_table_avg the OTHERS and SUB rows
pivot_table_avg = pivot_table_avg[~pivot_table_avg['Event_around_cut'].isin(['OTHERS', 'SUB'])]
pivot_table_avg
# %%
file_number = '001'
counter = 2
text_bellow_chart = f'INDEL category clasified only for events around the cut - full data in the Supplement.\
                    \nIN+DEL+SUB+OTHERS+no_mutation = 100%.\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_INDEL_summary_violinplot'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='Event_around_cut', y='value'\
        , number_of_columns=2\
        , project_output_dir=project_output_dir\
        , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False, show_hue_medians=False\
        , hue='Event_around_cut', hue_only_for_data_points=True, transparency=0.7 \
        , clip_right=True, change_edge_color=False\
        , title=title, xtickt_rotation=0, ylabel=ylabel\
        , order=['IN', 'DEL'] \
        , integrated_properties_dic=True\
        , text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05\
        , yaxis_formater=100, name_to_save=name_to_save\
        , data_source=data_source, output_format=output_format\
        , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# %%
