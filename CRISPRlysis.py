# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
This code will runn all the analysis based on the parameters.txt file which must be located in the same folder as this code CRISPRlysis.py
Output files will be stored in the 'project_output_dir' which is defined in the parameters.txt
Example: 'project_output_dir':'../../CRISPRlysis_outputs/Pribylova_et_al_2022'
"""
# %%
# TODO
# TODO do not use DEL_orientation_sober, but use the DEL-orientation from apply_return_cut_type_based_on_cut_categories
# %%
import ast
import csv
from datetime import datetime
from glob import glob
import os
from os.path import dirname, split, basename
from pathlib import Path
import shutil
import time
import sys

import itertools
import matplotlib
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.patches import PathPatch
import matplotlib.transforms as transforms
import more_itertools as mit
import numpy as np
import pandas as pd
import seaborn as sns
from statannot import add_stat_annotation
from scripts.general_functions import write_into_report, copy_files_with_parent_folder_to_another, find_files_paths_with_criteria, check_if_path_exist_create_if_not

from scripts.make_mastertables import filter_out_values_from_df, clasify_INDELSUB_from_columns, change_Event_based_on_mut_position, find_ranges, count_no_of_gaps, find_DEL_closest_to_cut, find_min_and_max_gap_ranges, find_DEL_around_cut_range, gap_span_analysed_window, find_microhomology, get_insertion_event, gap_in_aligned_seq, gap_in_range, find_difference, mutation_in_distance, make_original_seq_dic, make_seq_around_optimum_cut, size_of_del_around_the_cut, no_of_gRNA_samples_txt, return_shifted_DEL_range, return_category_of_DEL, return_cut_type_based_on_cut_categories, process_table, create_mastertable, add_columns_to_df, make_subtables_df

from scripts.plot_functions import save_stat, make_Ubuntu_mnt_path, make_reverse_complement, find_element_in_list, adjust_box_widths, plot_sns_violin_and_print, plot_general_stacked_horizontal_barchart, plot_barplot, check_if_all_events_are_present_after_pivot_reset_index

from scripts.dict_functions import make_color_dic_CutType, make_hatches_dic_CutType, make_color_dic_blunt_stagg_orientation

# %%
# Get current date and time
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')

# Load parameters.txt as dictionary
with open('./parameters.txt', 'r', encoding='utf-8') as file:
    file_content = file.read()
    parameters = ast.literal_eval(file_content)

# save parameters as variables
source_dir = parameters['source_dir']
project_output_dir = parameters['project_output_dir']
add_columns = parameters['add_columns']
data_source_original = parameters['data_source_original']
cut_position = parameters['cut_position']
cut_maximum = parameters['cut_maximum']
cut_minimum = parameters['cut_minimum']
nuclease = parameters['nuclease']
minimum_Reads = parameters['minimum_Reads']
output_format = parameters['output_format']
make_mastertable_00 = parameters['make_mastertable_00']
path_to_existing_mastertable_00 = parameters['path_to_existing_mastertable_00']
make_mastertable_01 = parameters['make_mastertable_01']
path_to_existing_mastertable_01 = parameters['path_to_existing_mastertable_01']
make_subtables = parameters['make_subtables']
subtables_param = parameters['subtables_param']
filter_out_values = parameters['filter_out_values']
show_plots = parameters['show_plots']
plot_INDELs = parameters['plot_INDELs']
plot_INs = parameters['plot_INs']
plot_DELs = parameters['plot_DELs']
plot_DEL_dinochart = parameters['plot_DEL_dinochart']
plot_MH = parameters['plot_MH']
plot_MH_dinochart = parameters['plot_MH_dinochart']
IN_greater_than = parameters['IN_greater_than']
sort_alphabetically = parameters['sort_alphabetically']

# Get the current path to the .py file
current_path = os.path.dirname(os.path.abspath(__file__))

# Check if project_output_dir exists, if not, create it
if not os.path.exists(project_output_dir):
    os.makedirs(project_output_dir)
    write_into_report(f'Report from {data_source_original} from date {date}.', destination_folder=project_output_dir)
    write_into_report(f'Project output dir was created and is {os.path.abspath(project_output_dir)}.', destination_folder=project_output_dir)
else:
    # Check if report.txt file exists and if yes, wipe it and make clear one
    if os.path.exists(Path(project_output_dir, 'report.txt')):
        with open(Path(project_output_dir, 'report.txt'), 'w') as file:
            pass
        # Make a new clean report.txt file
        write_into_report(f'Report from {data_source_original} from date {date}.', destination_folder=project_output_dir)
        write_into_report(f'Project output dir already exists and is {os.path.abspath(project_output_dir)}.', destination_folder=project_output_dir)
    else:
        # Make a new clean report.txt file
        write_into_report(f'Report from {data_source_original} from date {date}.', destination_folder=project_output_dir)
        write_into_report(f'Project output dir already exists and is {os.path.abspath(project_output_dir)}.', destination_folder=project_output_dir)

# TODO delete in the final version of the code
write_into_report(f'Cut position is set to location: {cut_position}', destination_folder=project_output_dir)

# Check if project output dir exists, create it if not
check_if_path_exist_create_if_not(project_output_dir)

# Make new mastertable from CRISPResso2 outputs or load existing one
if make_mastertable_00 is True:
    create_mastertable(source_dir, project_output_dir, parameters)
    # create_mastertable does not return the DataFrame, thus the table must be loaded later
    mastertable = pd.read_csv(Path(str(project_output_dir) + '/00_mastertable_' + date + '.tsv'), sep='\t')
    write_into_report(f'00_mastertable shape: {mastertable.shape}', destination_folder=project_output_dir)
else:
    mastertable = pd.read_csv(path_to_existing_mastertable_00, sep='\t')
    write_into_report(f'Reading 00_mastertable from: {path_to_existing_mastertable_00}. Shape {mastertable.shape}', destination_folder=project_output_dir)

# Filter out values from mastertable
if filter_out_values:
    mastertable = filter_out_values_from_df(mastertable, filter_out_values, project_output_dir)

# Make new mastertable containing additional columns from 00_mastertable or load existing one
if make_mastertable_01 is True:
    df = add_columns_to_df(mastertable, project_output_dir, parameters, date)
else:
    df = pd.read_csv(path_to_existing_mastertable_01, sep="\t")
    write_into_report(f'Reading 01_mastertable from: {path_to_existing_mastertable_01}', destination_folder=project_output_dir)

if make_subtables is True:
    tables_dic = make_subtables_df(df, project_output_dir, subtables_param, date)
else:
    tables_dic = {path_to_existing_mastertable_01:subtables_param}
# %%
tables_dic
# %%
# if tables_dic is not empty (imported empty from variables_for_analysis_linux.py) it will use the first half of the code.
# else it will use the bottom part of the code, where there are no sub-tables

for key, value in tables_dic.items():
    print(f'Processing {key} and {value}')
    df = pd.read_csv(key, sep='\t', index_col=None)
    write_into_report('###########################\n# Analysing next subtable #\n###########################', destination_folder=project_output_dir)
    write_into_report(f'Processing {key}', destination_folder=project_output_dir)
    write_into_report(f'{df["gRNA_seq"].unique()} are unique gRNAs in file.', destination_folder=project_output_dir)

    # TODO the part bellow might be deleted as the filtering should be done between mastertable 00 and 01 by 'filter_out_values':{} from the parameters.txt file
    # # filter out some gRNA
    # if filter_out_some_gRNA == True:
    #     df = filter_gRNA(df, list_of_gRNA_seq_to_filter)
    # write_into_report(f'After filtering {df["gRNA_seq"].unique()} are unique gRNAs in file.', destination_folder=project_output_dir)

    # Define new data_dource for individual key
    i=0
    data_source = f'{data_source_original} \n'
    for param in subtables_param:
        data_source += f'{param} is {value[i]}; '
        i+=1
    write_into_report(f'New data_source is {data_source}.', destination_folder=project_output_dir)
    # Change project_output_dir
    project_output_subdir = Path(project_output_dir, f'{value[0]}_{value[1]}_{date}')
    check_if_path_exist_create_if_not(project_output_subdir)
    write_into_report(f'New project_output_subdir is {project_output_subdir}', destination_folder=project_output_dir)
    # Sort Alphabetically
    if sort_alphabetically is True:
        order_list = sorted(df['gRNA_seq'].unique().tolist())
        # TODO delete if the above row is sorting properly
        # a = df['gRNA_seq'].unique()
        # order_list = a.tolist()
        # order_list.sort()
        write_into_report(f'df was sorted alphabetically. The order is: {order_list}', destination_folder=project_output_dir)
    else:
        order_list = df['gRNA_seq'].unique()
    if len(order_list) == 0:
        write_into_report('There is no gRNA to be plotted, moved to the next subtable.', destination_folder=project_output_dir)
        continue

    no_of_gRNA_samples = no_of_gRNA_samples_txt(df, order_list, project_output_dir)

    ##########################################################################################
    #            I   N    N   DDDD    EEEEEE  L        SSSS   U   U   BBBB
    #            I   N N  N   D   D   E       L       S       U   U   B   B
    #            I   N  N N   D   D   EEEE    L        SSS    U   U   BBBB
    #            I   N   NN   D   D   E       L           S   U   U   B   B
    #            I   N    N   DDDD    EEEEEE  LLLLLL  SSSS     UUU    BBBB
    ##########################################################################################
    if plot_INDELs is True:
        file_number = '01'
        counter = 1
        # Get data
        df_INDELSUB = df[['order', 'sample', 'gRNA_seq', 'Event_around_cut', '%Reads', 'name']]
        list_of_others = ['INDEL', 'INSUB', 'INDELSUB', 'DELSUB']

        pivot = pd.pivot_table(df_INDELSUB, index=['sample'], columns='Event_around_cut', values='%Reads', aggfunc='sum', fill_value=0)
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

        df_INDELSUB_summ = df_INDELSUB_summ[['sample', 'IN', 'DEL', 'SUB', 'OTHERS']]
        df_INDELSUB_summ = pd.melt(df_INDELSUB_summ, id_vars=['sample']\
                                , value_vars=['IN', 'DEL', 'SUB', 'OTHERS'])
        ##########################################################################################
        # plot VIOLINPLOT
        ##########################################################################################
        # Set charts properties
        title = 'INDELSUB representation in all samples'
        ylabel = 'Percentage'
        no_of_samples = len(df_INDELSUB_summ['sample'].unique())
        text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                            \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                            \nIN+DEL+SUB+OTHERS+no_mutation = 100%.\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_INDELSUB_summary_violinplot'
        plot = plot_sns_violin_and_print(data=df_INDELSUB_summ, x='Event_around_cut', y='value'\
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
        df_INDELSUB_summ.to_csv(Path(str(project_output_subdir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
        stat_info = df_INDELSUB_summ.groupby('Event_around_cut').value.describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1

        ##########################################################################################
        # plot HORIZONTAL BAR PLOT for individual samples
        ##########################################################################################
        title = 'INDELSUB analysis not normalised'
        xlabel = 'Average of %Reads from individual samples'
        name_to_save = f'{file_number}_{counter}_INDELSUB_horizontal_barchart_all_samples_not_normalised'
        text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                        \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                        \nIN+DEL+SUB+OTHERS+no_mutation = 100%.\
                        \n{no_of_gRNA_samples}'
        INDELSUB_color_dic=  {'IN':'#091f92'\
                            ,'DEL':'#ba8759'\
                            ,'SUB':'#ffcc33'\
                            ,'OTHERS':'#808080'}
        plot = plot_general_stacked_horizontal_barchart(df=df_INDELSUB_summ, x='Event_around_cut', y='sample', values='value', project_output_dir=project_output_subdir, data_source=data_source, title=title, name_to_save=name_to_save, output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False, palette=INDELSUB_color_dic, x_as_perc=True, text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default', xaxis_formater=100, decimals=0, xlim=(0,100), bars_order=['IN', 'DEL', 'SUB', 'OTHERS'], show_plots=show_plots)
        counter+=1

        # normalise the events for normalised plot form 0-100%
        pivot = pd.crosstab(index=df_INDELSUB_summ['sample'],\
                    columns=df_INDELSUB_summ['Event_around_cut'], \
                    values=df_INDELSUB_summ['value'],\
                    aggfunc='sum', normalize='index')
        df_INDELSUB_samples_norm = pivot.reset_index()
        df_INDELSUB_samples_norm = pd.melt(df_INDELSUB_samples_norm, id_vars=['sample'], value_vars=['IN', 'DEL', 'SUB', 'OTHERS'])
        title = 'INDELSUB analysis normalised'
        name_to_save = f'{file_number}_{counter}_INDELSUB_horizontal_barchart_all_samples_normalised'
        text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                            \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                            \nIN+DEL+SUB+OTHERS = 100%.\
                            \n{no_of_gRNA_samples}'
        # TODO solve y_legend_cooridnate
        # for 3 rows -1 is fine
        plot = plot_general_stacked_horizontal_barchart(df=df_INDELSUB_samples_norm, x='Event_around_cut', y='sample', values='value', project_output_dir=project_output_subdir, data_source=data_source, title=title, name_to_save=name_to_save, output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False, palette=INDELSUB_color_dic, x_as_perc=True, text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default', xaxis_formater=1, decimals=0, xlim=(0,1), bars_order=['IN', 'DEL', 'SUB', 'OTHERS'], show_plots=show_plots)
        stat_info = df_INDELSUB_samples_norm.groupby('Event_around_cut').value.describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1

    ##########################################################################################
    #            I   N    N
    #            I   N N  N
    #            I   N  N N
    #            I   N   NN
    #            I   N    N
    ##########################################################################################
    if plot_INs is True:
        file_number = '02'
        counter = 1
        # Get data
        df_IN = df
        pivot = pd.pivot_table(df_IN, index=['name', 'n_inserted', '4th_nt'], columns=(), values='%Reads', aggfunc='sum', fill_value=0)
        df_pivot = pivot.reset_index()
        df_pivot = df_pivot[df_pivot['n_inserted']!=0]

        ##########################################################################################
        # plot BARPLOT with number of inserted nucleotides
        ##########################################################################################
        name_to_save = f'{file_number}_{counter}_No_of_inserted_nucleotides'
        title = 'No of inserted nucleotides'

        # Plot average percentage of insertion and drop all despite 0 insertions
        plot = plot_barplot(data=df_pivot, x='n_inserted', y='%Reads', project_output_dir=project_output_subdir, name_to_save=name_to_save, data_source=data_source, hue=None, color='#091f92', title=title, x_lim=(-0.5,14.5), y_lim=(0,55), xlabel='No of nucleotides', yaxis_formater=100, decimals=1, errorbar='sd', fill_empty_x=True, output_format=output_format, dpi=100, show_plots=show_plots)
        stat_info = df_pivot.groupby('n_inserted').describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1

        ##########################################################################################
        # SUM up IN>3 nt to single column, to change the number, change grater_than
        ##########################################################################################
        # next two lines must be there, without it the former plot changes the df_pivot
        # and currently I do not why
        # TODO find, why plot_barplot when fill_empty_x=True changes the df_pivot
        df_pivot = pivot.reset_index()
        df_pivot = df_pivot[df_pivot['n_inserted']!=0]

        pivot = pd.crosstab(index=df_pivot['name'], columns=df_pivot['n_inserted'],\
                            values=df_pivot['%Reads'], aggfunc='sum', normalize='index')
        df_pivot = pivot.reset_index()

        # sum columns where n_inserted>IN_greater_than and discard n_inserted==0
        greater_than = IN_greater_than
        df_pivot[f'greater_than_{greater_than}'] = 0
        for column in df_pivot.columns:
            if isinstance(column, str):
                continue
            if column < (greater_than+1):
                continue
            df_pivot[f'greater_than_{greater_than}'] += df_pivot[column]
            df_pivot.drop(labels=column, axis=1, inplace=True)
        columns = df_pivot.iloc[:0,1:]
        df_pivot = df_pivot.melt(id_vars=['name'], value_vars=columns)
        ##########################################################################################
        # plot VIOLINPLOT
        ##########################################################################################
        # Set charts properties
        color_dic = {'0':'#ffffff', '1':'#A8B7FF', '2':'#536CE8', '3':'#466270', '4':'#514C7A', '5':'#4E8387'\
                    , '6':'#487D72', '7':'#4E9472'\
                    , 0:'#ffffff', 1:'#A8B7FF', 2:'#536CE8', 3:'#466270', 4:'#514C7A', 5:'#4E8387'\
                    , 6:'#487D72', 7:'#4E9472'}
        color_dic[f'greater_than_{IN_greater_than}'] = '#6F79A8'
        title='Normalised numbers od inserted nucleotides around the cut site'
        ylabel='Normalised percentage\n of inserted nucleotides'
        no_of_samples = len(df_pivot['name'].unique())
        text_bellow_chart = f'Only IN events around cut.\
                            \nSum of all IN in one sample = 100%\
                            \n{no_of_gRNA_samples}'
        name_to_save=f'{file_number}_{counter}_IN_no_of_IN_nucleotides_summary_violinplot'

        plot = plot_sns_violin_and_print(data=df_pivot, x='n_inserted', y='value', number_of_columns=len(df_pivot['n_inserted'].unique()), show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, order=None, project_output_dir=project_output_subdir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=title, xtickt_rotation=0, ylabel=ylabel, yaxis_formater=1, name_to_save=name_to_save, text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=1.02, output_format=output_format, show_plots=show_plots)
        stat_info = df_pivot.groupby('n_inserted').value.describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1

        title = 'Number od inserted nucleotides'
        xlabel='Normalised percentage of inserted nucleotides'
        name_to_save=f'{file_number}_{counter}_no_of_IN_nucleotides_stacked_barchart_name'
        # palette_dic = {0:'#ffffff', 1:'#A8B7FF', 2:'#536CE8', 3:'#424763', 'greater_than_2':'#6F79A8', 'greater_than_3':'#6F79A8'}
        plot = plot_general_stacked_horizontal_barchart(df=df_pivot, x='n_inserted', y='name', project_output_dir=project_output_subdir, values='value', sns_palette=False, palette=color_dic, data_source=data_source, title=title, name_to_save=name_to_save, output_format=output_format, dpi=200, xlabel=xlabel, x_as_perc=True, xaxis_formater=1, decimals=0, xlim=(0,1), text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default', y_order=None, bars_order=None, show_plots=show_plots)
        counter+=1
        #########################################################################################
        #########################################################################################
        #########################################################################################
        # Plot if the 4th nucleotide the inserted nucleotide
        # Get data
        df_IN = df[df['Event_around_cut'] == 'IN']
        df_IN = df[df['n_inserted'] == 1]

        pivot = pd.pivot_table(df_IN, index=['sample', '4th_nt_inserted', '4th_nt'], values='%Reads', aggfunc='sum', fill_value=0)
        df_pivot = pivot.reset_index()

        ##########################################################################################
        # plot STACKED HORIZONTAL BARCHART for 4th nt inserted True/False
        ##########################################################################################
        color_dic={True:'#06145E', False:'#B8BEDB'}
        text_bellow_chart = 'Only 1nt IN around the cut are plotted.\nSequences without DEL/SUB around the cut.'
        title='Is the 4th nucleotide the inserted nucleotide?'
        name_to_save=f'{file_number}_{counter}_4th_nt_inserted_not_normalised'
        xlabel='Sum of %Reads from individual samples'
        plot = plot_general_stacked_horizontal_barchart(df=df_pivot, y='sample'\
            , x='4th_nt_inserted', project_output_dir=project_output_subdir\
            , values='%Reads', sns_palette=False, palette=color_dic\
            , hatches_dic=None\
            , data_source=data_source\
            , title=title\
            , name_to_save=name_to_save, output_format=output_format\
            , dpi=200, xlabel=xlabel, x_as_perc=True, xaxis_formater=100\
            , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default'\
            , bars_order=[True, False], xlim=(0,100), show_plots=show_plots)
        stat_info = df_pivot.groupby('4th_nt_inserted').describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
        ##########################################################################################
        ## Normalise to have True+False = 100%
        pivot = pd.crosstab(index=[df_IN['sample'], df_IN['4th_nt']],\
                    columns=df_IN['4th_nt_inserted'], \
                    values=df_IN['%Reads'],\
                    aggfunc='sum', normalize='index')
        df_pivot = pivot.reset_index()

        list_to_check = [True, False]
        for item in list_to_check:
            if item not in df_pivot.columns.tolist():
                df_pivot[item] = int()

        df_pivot = pd.melt(df_pivot, id_vars=(['sample', '4th_nt']), value_vars=[True, False])

        title='Is the 4th nucleotide the inserted nucleotide? Normalised'
        name_to_save=f'{file_number}_{counter}_4th_nt_inserted_normalised'
        plot = plot_general_stacked_horizontal_barchart(df=df_pivot, y='sample'\
            , x='4th_nt_inserted', project_output_dir=project_output_subdir\
            , values='value', sns_palette=False, palette=color_dic\
            , hatches_dic=None\
            , data_source=data_source\
            , title=title\
            , name_to_save=name_to_save, output_format=output_format\
            , dpi=200, xlabel=(xlabel+' normalised')\
            , x_as_perc=True, xaxis_formater=1\
            , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default'\
            , bars_order=[True, False], xlim=(0,1), show_plots=show_plots)
        stat_info = df_pivot.groupby('4th_nt_inserted').describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1

        df_pivot = df_pivot[df_pivot['4th_nt_inserted']==True]
        # make hue_order according to present 4th nucleotides
        nucleotide_order_list = {'A':1, 'C':2, 'G':3, 'T':4}
        nucleotide_list = df_pivot['4th_nt'].unique()
        sort_list = []
        for nucleotide in nucleotide_list:
            sort_list.append(nucleotide_order_list[nucleotide])
        hue_order = [x for _, x in sorted(zip(sort_list, nucleotide_list))]        
        ##########################################################################################
        # plot VIOLINPLOT 1
        ##########################################################################################
        # Set charts properties
        # palette_dic=  {True:'#D3D3D3'
        #             , False:'#B8BEDB'
        #             , 'A':'#70c936'
        #             , 'C':'#ff0000'
        #             , 'G':'#ffd42a'
        #             , 'T':'#9955ff'}
        title='Is the 4th nucleotide the inserted nucleotide?'
        no_of_samples = len(df_pivot['sample'].unique())
        ylabel='Normalised percentage\n of inserted 4th nucleotides'
        text_bellow_chart = f'Only 1nt IN around the cut are plotted.\
                            \nSequences without DEL/SUB around the cut.\
                            \nIn individual sample True + False = 100%.\
                            \nHorizontal lines show medians for individual 4th nucleotide.\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_IN_4th_nt_inserted_normalised_violinplot'
        # TODO Maybe dodge should be True here
        plot = plot_sns_violin_and_print(data=df_pivot, x='4th_nt_inserted', y='value', number_of_columns=1.2\
            , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=False\
            , show_hue_medians=False\
            , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
            , clip_right=True\
            , title=title, xtickt_rotation=0, ylabel=ylabel\
            , order=None, plot_legend=True \
            , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.05, yaxis_formater=1\
            , x_max=1\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir \
            , data_source=data_source, output_format=output_format\
            , dpi=200\
            , y_min=0\
            , y_max=1.02, show_plots=show_plots
                )
        stat_info = df_pivot.groupby('4th_nt_inserted').describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
        ##########################################################################################
        # plot VIOLINPLOT 2
        ##########################################################################################
        name_to_save = f'{file_number}_{counter}_IN_4th_nt_inserted_normalised_violinplot'
        # TODO seems median line color is not correct and legend for stripplot is not shown
        plot = plot_sns_violin_and_print(data=df_pivot, x='4th_nt_inserted', y='value', number_of_columns=1\
            , show_datapoints=True, jitter=0.25, stripplot_offset=0, dodge=True\
            , show_hue_medians=True\
            , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
            , clip_right=False\
            , title=title, xtickt_rotation=0, ylabel=ylabel\
            , order=None, plot_legend=True \
            , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.05, yaxis_formater=1\
            , x_max=None\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir \
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
        stat_info = df_pivot.groupby('4th_nt').value.describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
    ##########################################################################################
    #               DDDD    EEEEEE  L
    #               D   D   E       L
    #               D   D   EEEE    L
    #               D   D   E       L
    #               DDDD    EEEEEE  LLLLLL
    ##########################################################################################
    # %%
    if plot_DELs is True:
        file_number = '03'
        counter = 1
        print('Plotting DELs')
        # boxplot from all individual samples showing DEL orientation
        df_orientation = df[['order', 'gRNA_seq', 'sample', 'name', 'DEL_orientation'\
                    , '%Reads', 'n_deleted_around_cut', 'Gap_around_cut']]
        # drop data without gap
        if 'X' in df_orientation['DEL_orientation'].unique():
            df_orientation = df_orientation[~df_orientation['DEL_orientation'].str.contains('X')]
        # filter out deletions which are not around the cut
        df_orientation = df_orientation[df_orientation['Gap_around_cut']==True]
        # Change ['U', 'B', 'D', 'C'] in DEL_orientation for ['upstream', 'both', 'downstream', 'cannot_be_determined']
        # Replace values in the DEL_orientation column
        df_orientation['DEL_orientation'] = df_orientation['DEL_orientation'].replace({
            'U': 'upstream',
            'B': 'both',
            'D': 'downstream',
            'C': 'cannot_be_determined'
        })
        # TODO delete order once it works
        # order = ['upstream', 'upstream_both', 'upstream_both_downstream', 'both', 'both_downstream', 'downstream']
        order = ['upstream', 'both', 'downstream', 'cannot_be_determined']
        ##########################################################################################
        # plot VIOLINPLOT
        ##########################################################################################
        # make pivot without normalisation
        pivot = pd.crosstab(index=[df_orientation['name'], df_orientation['gRNA_seq']], \
                columns=df_orientation['DEL_orientation'], values=df_orientation['%Reads'],\
                aggfunc='sum')
        df_pivot = pivot.reset_index()
        # check if all columns exists, if not make an empty one
        for orientation in order:
            if orientation not in df_pivot.columns.tolist():
                df_pivot[orientation] = int()

        columns = df_pivot.iloc[:0,2:]
        df_pivot = df_pivot.melt(id_vars=['name', 'gRNA_seq'], \
            value_vars=columns)

        # Set charts properties
        title = 'Not-normalised incidence of DEL according to the orientation'
        ylabel = 'Percentage of DEL around the cut'
        no_of_samples = df_pivot.groupby('DEL_orientation').value.count().unique()
        text_bellow_chart = f'Only DELs around the cut are plotted.\
                            \nNot-normalised\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_1_DEL_orientation_summary_violinplot'
        plot = plot_sns_violin_and_print(data=df_pivot, x='DEL_orientation', y='value'\
            , number_of_columns=len(df_pivot['DEL_orientation'].unique())/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=True\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=True, color_dic=None\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.4, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=True\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=100, show_plots=show_plots)
        stat_info = df_pivot.groupby('DEL_orientation').value.describe()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=0 # Should be 0 as next chart is just normalised
        ##########################################################################################
        # plot VIOLINPLOT normalised
        ##########################################################################################
        # normalise data as percentage from row total
        pivot = pd.crosstab(index=[df_orientation['name'], df_orientation['gRNA_seq']]\
            , columns=df_orientation['DEL_orientation'], values=df_orientation['%Reads']\
            , aggfunc='sum', normalize='index')
        df_pivot = pivot.reset_index()

        # check if all columns exists, if not make an empty one
        order = ['upstream', 'both', 'downstream', 'cannot_be_determined']
        # TODO delete order once it works
        # order = ['upstream', 'upstream_both', 'upstream_both_downstream', 'both', 'both_downstream', 'downstream']
        for orientation in order:
            if orientation not in df_pivot.columns.tolist():
                df_pivot[orientation] = int()

        columns = df_pivot.iloc[:0,2:]
        df_pivot = df_pivot.melt(id_vars=['name', 'gRNA_seq'], \
            value_vars=columns)

        # Set charts properties
        title = 'Normalised incidence of DEL according to the orientation'
        ylabel = 'Normalised percentage of DEL around the cut'
        no_of_samples = df_pivot.groupby('DEL_orientation').value.count().unique()
        text_bellow_chart = f'Only DELs around the cut are plotted.\
                            \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_2_DEL_orientation_summary_violinplot'
        plot = plot_sns_violin_and_print(data=df_pivot, x='DEL_orientation', y='value'\
            , number_of_columns=len(df_pivot['DEL_orientation'].unique())/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=True\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=True, color_dic=None\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.4, yaxis_formater=1\
            , x_max=None\
            , change_edge_color=True\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
        stat_info = df_pivot.groupby('DEL_orientation').value.describe()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        # %%
        # Get data
        df_CutType = df.copy()
        ## drop data without gap
        df_CutType = df_CutType[df_CutType['n_deleted_around_cut'] != 0]
        ## filter out deletions which are not around the cut
        df_CutType = df_CutType[df_CutType['Gap_around_cut'] == True]
        # make pivot, where all samples have individual row and values are sum of individual events
        pivot = pd.pivot_table(df_CutType, index=['name', 'gRNA_seq'], columns=['Cut_type', 'Likely_cut_type'], values='%Reads', aggfunc='sum', fill_value=0)
        df_Blunt_Stagg_samples = pivot.reset_index()
        # Flatten the index
        df_Blunt_Stagg_samples.columns = ['_'.join(map(str, col)).strip('_') for col in df_Blunt_Stagg_samples.columns]
        ## check if all events are present, if not, add missing column with 0
        df_Blunt_Stagg_samples, list_of_events = check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, 'BluntStag', cut_maximum, cut_minimum)
        # TODO ensure C_1, O_1 and X_1 does not exists in df and convert C_0 to C only etc.
        # make single column table
        df_Blunt_Stagg_samples = pd.melt(df_Blunt_Stagg_samples, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
        pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq']], columns=df_Blunt_Stagg_samples['variable'], values=df_Blunt_Stagg_samples['value'], aggfunc='sum', normalize='index')
        df_Blunt_Stagg_samples_norm = pivot.reset_index()
        df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
        # Set charts properties
        # TODO is order needed bellow?
        # order=list_of_events
        ##########################################################################################
        # Horizontal bar charts for Blunt/Staggered cuts
        ##########################################################################################
        # normalise the events for normalised plot form 0-100%
        title = 'Staggered/Blunt cut ratio for DEL adjacent or around the cut site'
        xlabel = 'Average of %Reads from individual samples'
        name_to_save = f'{file_number}_{counter}_Stagg_blunt_ratio_all_samples'

        # TODO solve color_dic and hatches_dic - new function was made in dict_functions.py
        color_dic = make_color_dic_CutType(list_of_events)
        # color_dic={\
        #     'stagg2nt':'#E890B8'\
        #     ,'stagg2nt_likely':'#E890B8'\
        #     ,'stagg1nt':'#ffaaaa'\
        #     ,'stagg1nt_likely':'#ffaaaa'\
        #     ,'blunt':'#6ea7ff'\
        #     ,'blunt_likely':'#6ea7ff'\
        #     ,'cannot_be_determined':'white'\
        #     ,'can_not_be_determined':'white'\
        #     ,'gap_out_of_cut':'gray'}
        # TODO solve type of hatches, especitally for blunt
        hatches_dic = make_hatches_dic_CutType(list_of_events)
        # hatches_dic = {'stagg2nt':''\
        #     ,'stagg2nt_likely':'//'\
        #     ,'stagg1nt':''\
        #     ,'stagg1nt_likely':'////'\
        #     ,'blunt':''\
        #     ,'blunt_likely':'\\\\'\
        #     ,'cannot_be_determined':''\
        #     ,'can_not_be_determined':''}

        text_bellow_chart = f'Only DELs around the cut are plotted.\
                        \nNormalised to: SUM of all DEL categories around the cut\n in individual sample = 100%.\
                        \n{no_of_gRNA_samples}'
        # TODO solve bars_order, color_dic, hatches_dic
        plot = plot_general_stacked_horizontal_barchart(df=df_Blunt_Stagg_samples_norm, x='variable'\
        , y='name', values='value', project_output_dir=project_output_subdir\
        , data_source=data_source, title=title, name_to_save=name_to_save\
        , output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False\
        , palette=color_dic, hatches_dic=hatches_dic, x_as_perc=True\
        , xaxis_formater=1, decimals=0, xlim=(0,1)\
        , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default'\
        , bars_order=list_of_events, show_plots=show_plots)
        stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
        # %%
        ##########################################################################################
        # plot VIOLINPLOT without PAM proximal/distal orientation
        ##########################################################################################
        title='Normalised Staggered/Blunt representation'
        ylabel='Normalised percentage of DELs \naround cuts events'
        text_bellow_chart = f'Only DELs around the cut are plotted.\
                            \nNormalised to: SUM of all DEL categories around the cut\n in individual sample = 100%.\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_DEL_Stagg_blunt_summary_violinplot'
        plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='variable', y='value'\
            , number_of_columns=len(df_Blunt_Stagg_samples_norm['variable'])/6\
            , show_datapoints=True, jitter=0.10, stripplot_offset=6, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel\
            , order=list_of_events\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=color_dic\
            , hatches=True, hatches_list=hatches_dic\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.15, yaxis_formater=1\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
        stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
        ##########################################################################################
        # plot VIOLINPLOTs with PAM proximal/distal orientation
        ##########################################################################################
        # %%
        # order = ["distal_stagg2nt", "distal_stagg2nt_likely", "distal_stagg1nt"\
        #     , "distal_stagg1nt_likely", "distal_blunt", "distal_blunt_likely"\
        #     , "proximal_blunt_likely", "proximal_blunt", "proximal_stagg1nt_likely"\
        #     , "proximal_stagg1nt", "proximal_stagg2nt_likely", "proximal_stagg2nt"\
        #     , "cannot_be_determined"]

        # make pivot, where all samples have individual row and values are sum of individual events
        # TODO for columns add the orientation
        pivot = pd.pivot_table(df_CutType, index=['name', 'gRNA_seq'], columns=['DEL_orientation', 'Cut_type', 'Likely_cut_type'], values='%Reads', aggfunc='sum', fill_value=0)
        # make order for pivot table based on the upstream/downstream, cut maximum/minimum and likeliness of the cut
        order = list(itertools.product(['U', 'D'], range(cut_maximum, cut_minimum+1), [0, 1]))
        order.append(('C', 0, '0'))
        order.append(('B', 'C', '0'))
        order = [(item[0], str(item[1]), item[2]) for item in order]
        pivot = pivot.reindex(columns=pd.MultiIndex.from_tuples(order), fill_value=0)
        # Flatten the order
        order = ['_'.join(map(str, col)).strip('_') for col in pivot.columns]	
        # Reset the index
        df_Blunt_Stagg_samples = pivot.reset_index()
        # Flatten the index
        df_Blunt_Stagg_samples.columns = ['_'.join(map(str, col)).strip('_') for col in df_Blunt_Stagg_samples.columns]

        # check if all events are present, if not, add missing column with 0
        # df_Blunt_Stagg_samples = check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, order, cut_maximum, cut_minimum)
        # # make single column table
        df_Blunt_Stagg_samples = pd.melt(df_Blunt_Stagg_samples, id_vars=['name', 'gRNA_seq'], value_vars=order)
        # %%
        df_Blunt_Stagg_samples
        # %%
        ##########################################################################################
        # plot VIOLINPLOT with PAM proximal/distal orientation with cannot_be_determined NOT-normalised
        ##########################################################################################
        # Make pivot crosstab without normalisation
        pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq']]\
                    , columns=df_Blunt_Stagg_samples['variable']\
                    , values=df_Blunt_Stagg_samples['value']\
                    , aggfunc='sum')
        df_Blunt_Stagg_samples_norm = pivot.reset_index()
        df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq'], value_vars=order)
        # %%
        # Set charts properties
        title='Not-normalised Staggered/Blunt representation, cut site oriented'
        ylabel='Percentage of DELs \naround cuts events'
        text_bellow_chart = f'Only DELs around the cut are plotted.\
                    \nNot-normalised.\
                    \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_1_DEL_Stagg_blunt_oriented_summary_violinplot'
        color_dic = make_color_dic_blunt_stagg_orientation(order)
        plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples, x='variable', y='value'\
            , number_of_columns=len(order)/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=6, dodge=True\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=True\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=color_dic, text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.20, yaxis_formater=100\
            , x_max=None\
            , plot_vertical_line_middle=True, middle_line_position=len(order)/2-1, change_edge_color=True\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
        # plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='DEL_orientation_cut_type', y='value'\
        #     , number_of_columns=6\
        #     , show_datapoints=True, jitter=0.10, stripplot_offset=6, dodge=True\
        #     , show_hue_medians=False\
        #     , hue=None, hue_only_for_data_points=False \
        #     , hatches=True\
        #     , clip_right=True\
        #     , title=title, xtickt_rotation=-55, ylabel=ylabel\
        #     , order=order\
        #     , plot_legend=False\
        #     , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
        #     , y_legend_coordinate=-0.40, yaxis_formater=100\
        #     , x_max=None\
        #     , plot_vertical_line_middle=True, middle_line_position=6, change_edge_color=True\
        #     , name_to_save=name_to_save\
        #     , project_output_dir=project_output_subdir\
        #     , data_source=data_source, output_format=output_format\
        #     , dpi=200, y_min=0, y_max=102)
        stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=0 # should be 0 as next chart is the same, just normalised
        ##########################################################################################
        # plot VIOLINPLOT with PAM proximal/distal orientation with cannot_be_determined
        ##########################################################################################
        # Make pivot crosstab where values will be normalised to sum of row = 100%
        pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq']]\
                    , columns=df_Blunt_Stagg_samples['variable']\
                    , values=df_Blunt_Stagg_samples['value']\
                    , aggfunc='sum', normalize='index')
        df_Blunt_Stagg_samples_norm = pivot.reset_index()
        df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq']\
                                    , value_vars=order)
        # Set charts properties
        title='Normalised Staggered/Blunt representation, cut site oriented'
        ylabel='Normalised percentage of DELs \naround cuts events'
        text_bellow_chart = f'Only DELs around the cut are plotted.\
                            \nNormalised to: SUM of all DEL categories around the cut\n in individual sample = 100%.\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'{file_number}_{counter}_2_DEL_Stagg_blunt_oriented_summary_violinplot_normalised'
        plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='variable'\
            , y='value', number_of_columns=len(order)/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=6, dodge=True\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=True, clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel\
            , order=order, plot_legend=False\
            , integrated_properties_dic=color_dic, text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.20, yaxis_formater=1\
            , x_max=None\
            , plot_vertical_line_middle=True, middle_line_position=len(order)/2-1, change_edge_color=True\
            , name_to_save=name_to_save\
            , project_output_dir=project_output_subdir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
        stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
        save_stat(stat_info, name_to_save, project_output_subdir)
        counter+=1
        # %%
# TODO make color_dic for all charts
color_dic = {'stagg2nt':{'color':'#E890B8','hatch_color':'gray', 'hatch':''}\
        ,'stagg2nt_likely':{'color':'#E890B8','hatch_color':'gray', 'hatch':'//'}\
        ,'stagg1nt':{'color':'#ffaaaa','hatch_color':'gray', 'hatch':''}\
        ,'stagg1nt_likely':{'color':'#ffaaaa','hatch_color':'gray', 'hatch':'////'}\
        ,'blunt':{'color':'#6ea7ff','hatch_color':'gray', 'hatch':''}\
        ,'blunt_likely':{'color':'#6ea7ff','hatch_color':'gray', 'hatch':'\\\\'}\
        ,'cannot_be_determined':{'color':'white','hatch_color':'gray', 'hatch':''}\
        ,'can_not_be_determined':{'color':'white','hatch_color':'gray', 'hatch':''}\
        ,'gap_out_of_cut':{'color':'gray','hatch_color':'', 'hatch':''}\
        ,'IN':{'color':'#091f92','hatch_color':'', 'hatch':''}\
        ,'INs':{'color':'#091f92','hatch_color':'', 'hatch':''}\
        ,'DEL':{'color':'#ba8759','hatch_color':'', 'hatch':''}\
        ,'DELs':{'color':'#ba8759','hatch_color':'', 'hatch':''}\
        ,'SUB':{'color':'#ffcc33','hatch_color':'', 'hatch':''}\
        ,'OTHERS':{'color':'#808080','hatch_color':'', 'hatch':''}\
        ,'no_mutation':{'color':'#111111','hatch_color':'', 'hatch':''}\
        ,'distal':{'color':'#f0d271','hatch_color':'', 'hatch':''}\
        ,'distal_both':{'color':'#C4C4C4','hatch_color':'', 'hatch':''}\
        ,'distal_both_proximal':{'color':'#979797','hatch_color':'', 'hatch':''}\
        ,'both':{'color':'#808080','hatch_color':'', 'hatch':''}\
        ,'both_proximal':{'color':'#C4C4C4','hatch_color':'', 'hatch':''}\
        ,'proximal':{'color':'#f5ad62','hatch_color':'', 'hatch':''}\
        ,'PAM_distal_nt_in_MH':{'color':'#52D8FF','hatch_color':'', 'hatch':''}\
        ,'PAM_proximal_nt_in_MH':{'color':'#FADA52','hatch_color':'', 'hatch':''}\
        ,'no_MH_likely':{'color':'#A9A9A9','hatch_color':'', 'hatch':''}\
        ,'no_MH':{'color':'#D3D3D3','hatch_color':'', 'hatch':''}\
        ,'A':{'color':'#70c936','hatch_color':'', 'hatch':''}\
        ,'C':{'color':'#ff0000','hatch_color':'', 'hatch':''}\
        ,'G':{'color':'#ffd42a','hatch_color':'', 'hatch':''}\
        ,'T':{'color':'#9955ff','hatch_color':'', 'hatch':''}}

# %%

