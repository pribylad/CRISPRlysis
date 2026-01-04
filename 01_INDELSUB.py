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
# data_source = 'dividing'
data_source = 'nondividing'
##########################################################################################
if data_source == 'dividing':
    data_dir = '/home/lefnui/CRISPRlysis_outputs/02_tables/dividing'
    project_output_subdir = '/home/lefnui/CRISPRlysis_outputs/03_charts/dividing'
else:
    data_dir = '/home/lefnui/CRISPRlysis_outputs/02_tables/nondividing'
    project_output_subdir = '/home/lefnui/CRISPRlysis_outputs/03_charts/nondividing'

# test if the directory exists, if not create it
if not os.path.exists(project_output_subdir):
    os.makedirs(project_output_subdir)
# Get current date and time
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')
##########################################################################################
def load_tsv_files_from_directory_to_dic(directory_path):
    tsv_files = glob.glob(os.path.join(directory_path, "02_*.tsv"))
    dataframes = {}
    for file in tsv_files:
        df_name = os.path.splitext(os.path.basename(file))[0]
        dataframes[df_name] = pd.read_csv(file, sep='\t')
    return dataframes

def summarize_dataframes(dataframes):
    summary = {}
    for name, df in dataframes.items():
        summary[name] = {
            'num_rows': df.shape[0],
            'num_columns': df.shape[1]
        }
    return summary

def add_metadata_column(df, column_name, value):
    df[column_name] = value
    return df

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
df = load_tsv_files_from_directory_to_dic(data_dir)
summary = summarize_dataframes(df)
# Join all dataframes into a single dataframe
df = pd.concat(df.values(), ignore_index=True)
##########################################################################################
#            I   N    N   DDDD    EEEEEE  L        SSSS   U   U   BBBB
#            I   N N  N   D   D   E       L       S       U   U   B   B
#            I   N  N N   D   D   EEEE    L        SSS    U   U   BBBB
#            I   N   NN   D   D   E       L           S   U   U   B   B
#            I   N    N   DDDD    EEEEEE  LLLLLL  SSSS     UUU    BBBB
##########################################################################################
file_number = '01'
counter = 1
# Define lists
default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'replicate']
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
no_of_samples = len(pivot_table_avg['gRNA_seq'].unique())
text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                    \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                    \nIN+DEL+SUB+OTHERS+no_mutation = 100%.\
                    \nNo of samples: {no_of_samples}'
name_to_save = f'{file_number}_{counter}_INDELSUB_summary_violinplot'
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='Event_around_cut', y='value'\
        , number_of_columns=3\
        , project_output_dir=project_output_subdir\
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
pivot_table_avg.to_csv(Path(str(project_output_subdir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
stat_info = pivot_table_avg.groupby('Event_around_cut').value.describe().unstack()
save_stat(stat_info, name_to_save, project_output_subdir)
