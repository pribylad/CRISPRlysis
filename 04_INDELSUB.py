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


#####
# filter_criterion = 'Kingdom'
# filter_value = 'Animals'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit5efficacy_Animals_260126.tsv'
# data_source = 'Animals'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/all_main_samples/{filter_criterion}_{filter_value}')
# #####
filter_criterion = 'Kingdom'
filter_value = 'Plants'
df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit5efficacy_Plants_260126.tsv'
data_source = 'Plants'
project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/all_main_samples/{filter_criterion}_{filter_value}')


#####

# filter_criterion = 'All'
# filter_value = 'All'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_efficacy_samples_limit5efficacy_260109.tsv'
# data_source = 'all_samples'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/all_main_samples/{filter_criterion}_{filter_value}')
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
#            I   N    N   DDDD    EEEEEE  L        SSSS   U   U   BBBB
#            I   N N  N   D   D   E       L       S       U   U   B   B
#            I   N  N N   D   D   EEEE    L        SSS    U   U   BBBB
#            I   N   NN   D   D   E       L           S   U   U   B   B
#            I   N    N   DDDD    EEEEEE  LLLLLL  SSSS     UUU    BBBB
##########################################################################################
file_number = '002_Supp'
counter = 1
# Define lists
default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism', 'replicate']
list_of_events = ['IN', 'DEL', 'SUB', 'OTHERS']

default_list_NOreplicate = default_list.copy()
default_list_NOreplicate.remove('replicate')
default_list_NOreplicate.append('Event_around_cut')

# Get data
df_INDELSUB = df[default_list + ['Event_around_cut', '%Reads']]
list_of_others = ['INDEL', 'INSUB', 'INDELSUB', 'DELSUB']

pivot = pd.pivot_table(df_INDELSUB, index=default_list, columns='Event_around_cut', values='%Reads', aggfunc='sum', fill_value=0)
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

df_INDELSUB_summ = df_INDELSUB_summ[default_list + list_of_events]
df_INDELSUB_summ = pd.melt(df_INDELSUB_summ, id_vars=default_list, value_vars=list_of_events)
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_INDELSUB_summ, default_list_NOreplicate, columns=None, values='value', aggfunc='mean')
pivot_table_avg.reset_index(inplace=True)
print(pivot_table_avg)
print(pivot_table_avg['source'].unique())

##########################################################################################
# plot VIOLINPLOT
##########################################################################################
# Set charts properties
title = 'INDELSUB representation in all samples, mean of replicates'
ylabel = 'Percentage'
no_of_samples = len(pivot_table_avg[['source', 'gRNA_seq', 'organism']].drop_duplicates())
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
file_number = '002_Fig'
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