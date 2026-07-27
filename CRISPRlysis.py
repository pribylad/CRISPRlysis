# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
This code will runn all the analysis based on the parameters.txt file which must be located in the same folder as this code CRISPRlysis.py
Output files will be stored in the 'project_output_dir' which is defined in the parameters.txt
Example: 'project_output_dir':'../../CRISPRlysis_outputs/Pribylova_et_al_2022'
"""

import ast
import csv
from datetime import datetime
from glob import glob
import os
from os.path import dirname, split, basename
from pathlib import Path
from random import randint
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
# TODO remove statannot once all works on statannotations bellow
# from statannot import add_stat_annotation
# from statannotations.Annotator import Annotator

from scripts.general_functions import write_into_report, copy_files_with_parent_folder_to_another, find_files_paths_with_criteria, check_if_path_exist_create_if_not

from scripts.make_mastertables import filter_out_values_from_df, clasify_INDELSUB_from_columns, change_Event_based_on_mut_position, find_ranges, count_no_of_gaps, find_DEL_closest_to_cut, find_min_and_max_gap_ranges, find_DEL_around_cut_range, gap_span_analysed_window, find_microhomology, get_insertion_event, gap_in_aligned_seq, gap_in_range, find_difference, mutation_in_distance, make_original_seq_dic, make_seq_around_optimum_cut, size_of_del_around_the_cut, no_of_gRNA_samples_txt, return_shifted_DEL_range, return_category_of_DEL, return_cut_type_based_on_cut_categories, process_table, create_mastertable, add_columns_to_df, make_subtables_df

from scripts.plot_functions import save_stat, make_Ubuntu_mnt_path, make_reverse_complement, find_element_in_list, adjust_box_widths, plot_sns_violin_and_print, plot_general_stacked_horizontal_barchart, plot_barplot, check_if_all_events_are_present_after_pivot_reset_index, plot_general_stacked_vertical_barchart

from scripts.dict_functions import make_color_dic_CutType, make_hatches_dic_CutType, make_color_dic_blunt_stagg_orientation, make_color_dic_Dino, make_hatches_dic_Dino, make_color_dic_Dino_end, make_hatches_dic_Dino_end, get_color_dic_35_violin

# Do not show future warnings in the output
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore", FutureWarning)

# Get start time
start_time = time.time()
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
if 'list_of_columns' in parameters:
    list_of_columns = parameters['list_of_columns']
else:
    list_of_columns = None
path_to_existing_mastertable_01 = parameters['path_to_existing_mastertable_01']
make_subtables = parameters['make_subtables']
path_to_existing_subtable_02 = parameters['path_to_existing_subtable_02']
subtables_param = parameters['subtables_param']
filter_out_values = parameters['filter_out_values']
show_plots = parameters['show_plots']
plot_INDELs = parameters['plot_INDELs']
plot_INs = parameters['plot_INs']
plot_DELs = parameters['plot_DELs']
plot_DEL_individual_gRNAs = parameters['plot_DEL_individual_gRNAs']
plot_DEL_dinochart = parameters['plot_DEL_dinochart']
plot_participating_end = parameters['plot_participating_end']
plot_participating_end_dinochart = parameters['plot_participating_end_dinochart']
# plot_MH = parameters['plot_MH']
# plot_MH_dinochart = parameters['plot_MH_dinochart'] 
IN_greater_than = parameters['IN_greater_than']
sort_alphabetically = parameters['sort_alphabetically']
# try if parameter 'random_gen_data' is in parameters, if not, set it to False
if 'random_gen_data' in parameters:
    random_gen_data = parameters['random_gen_data']
else:
    random_gen_data = False

# Get the current path to the .py file
current_path = os.path.dirname(os.path.abspath(__file__))

# Check if source_dir exists, if not, raise an error
if not os.path.exists(source_dir):
    raise FileNotFoundError(f'Source directory {source_dir} does not exist. Please check the path in parameters.txt and make sure it is correct.')

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

# Check if project output dir exists, create it if not
check_if_path_exist_create_if_not(project_output_dir)

# Make new mastertable from CRISPResso2 outputs or load existing one
if make_mastertable_00 is True:
    create_mastertable(source_dir, project_output_dir, parameters, random_gen_data)
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
    write_into_report('Making 01_mastertable.', destination_folder=project_output_dir)
    df = add_columns_to_df(mastertable, project_output_dir, parameters, date, sort_alphabetically, list_of_columns)
else:
    df = pd.read_csv(path_to_existing_mastertable_01, sep="\t")
    write_into_report(f'Reading 01_mastertable from: {path_to_existing_mastertable_01}', destination_folder=project_output_dir)

if make_subtables is True:
    tables_dic = make_subtables_df(df, project_output_dir, subtables_param, date)
else:
    if path_to_existing_subtable_02 != '':
        tables_dic = {path_to_existing_subtable_02:subtables_param}
        print(f'path_to_existing_subtable_02 is {path_to_existing_subtable_02}')
    else:
        tables_dic = {path_to_existing_mastertable_01:subtables_param}
# if tables_dic is not empty (imported empty from variables_for_analysis_linux.py) it will use the first half of the code.
# else it will use the bottom part of the code, where there are no sub-tables

for key, value in tables_dic.items():
    print(f'Processing {key} and {value}')
    df = pd.read_csv(key, sep='\t', index_col=None)
    write_into_report('###########################\n# Analysing next subtable #\n###########################', destination_folder=project_output_dir)
    write_into_report(f'Processing {key}', destination_folder=project_output_dir)

    # Define new data_dource for individual key
    i=0
    data_source = f'{data_source_original} \n'
    # if subtables_param is not empty list
    # TODO solve for filtering to other than 2 stuff
    if len(subtables_param) > 0:
        subpath = ''
        for param in subtables_param:
            subpath = subpath + str(value[i]) + '_'
            data_source += f'{param} is {value[i]}; '
            i+=1
        write_into_report(f'New data_source is {subpath}.', destination_folder=project_output_dir)
        # Change project_output_dir
        project_output_subdir = Path(project_output_dir, subpath + f'{date}')
        print(f'New project_output_subdir is {project_output_subdir}')
    else:
        project_output_subdir = Path(project_output_dir, f'{date}')
    check_if_path_exist_create_if_not(project_output_subdir)
    write_into_report(f'New project_output_subdir is {project_output_subdir}', destination_folder=project_output_dir)
    # Sort Alphabetically
    if sort_alphabetically is True:
        order_list = sorted(df['gRNA_seq'].unique().tolist())
        write_into_report(f'df was sorted alphabetically. The order is: {order_list}', destination_folder=project_output_dir)
    else:
        order_list = df['gRNA_seq'].unique()
    if len(order_list) == 0:
        write_into_report('There is no gRNA to be plotted, moved to the next subtable.', destination_folder=project_output_dir)
        continue

    no_of_gRNA_samples = no_of_gRNA_samples_txt(df, order_list, project_output_dir)

    # If column 'Event' has only rows with 'no_mutation' skip this cycle
    if df['Event'].nunique() == 1 and df['Event'].unique()[0] == 'no_mutation':
        write_into_report(f'There is no mutation in {key}, moved to the next subtable.', destination_folder=project_output_dir)
        continue

    ##########################################################################################
    #            I   N    N   DDDD    EEEEEE  L        SSSS   U   U   BBBB
    #            I   N N  N   D   D   E       L       S       U   U   B   B
    #            I   N  N N   D   D   EEEE    L        SSS    U   U   BBBB
    #            I   N   NN   D   D   E       L           S   U   U   B   B
    #            I   N    N   DDDD    EEEEEE  LLLLLL  SSSS     UUU    BBBB
    ##########################################################################################
    if plot_INDELs is True:
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
        ylabel = 'Percentage of reads'
        no_of_samples = len(df_INDELSUB_summ['sample'].unique())
        text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                            \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                            \nIN+DEL+SUB+OTHERS+no_mutation = 100%.\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'01_01_INDELSUB_summary_violinplot'
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

        ##########################################################################################
        # plot HORIZONTAL BAR PLOT for individual samples
        ##########################################################################################
        # If more then 40 samples are present, divide the plot into as many plots as needed to have maximum 40 samples in one plot, if not, make one plot with all samples
        # count df_INDELSUB_summ['sample'].unique() and if more than 40, make a list of lists with maximum 40 unique samples in each list
        unique_samples = df_INDELSUB_summ['sample'].unique()
        if len(unique_samples) > 40:
            sample_lists = [unique_samples[i:i+40] for i in range(0, len(unique_samples), 40)]
        else:
            sample_lists = [unique_samples]
        i=1
        for sample_list in sample_lists:
            df_INDELSUB_summ_sample = df_INDELSUB_summ[df_INDELSUB_summ['sample'].isin(sample_list)]
            title = 'INDELSUB analysis not normalised'
            xlabel = 'Average of %Reads from individual samples'
            name_to_save = f'01_02_{i}_INDELSUB_horizontal_barchart_all_samples_not_normalised'
            text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                            \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                            \nIN+DEL+SUB+OTHERS+no_mutation = 100%.'
            INDELSUB_color_dic=  {'IN':'#091f92'\
                                ,'DEL':'#ba8759'\
                                ,'SUB':'#ffcc33'\
                                ,'OTHERS':'#808080'}
            plot = plot_general_stacked_horizontal_barchart(df=df_INDELSUB_summ_sample, x='Event_around_cut', y='sample', values='value', project_output_dir=project_output_subdir, data_source=data_source, title=title, name_to_save=name_to_save, output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False, palette=INDELSUB_color_dic, x_as_perc=True, text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default', xaxis_formater=100, decimals=0, xlim=(0,100), bars_order=['IN', 'DEL', 'SUB', 'OTHERS'], show_plots=show_plots)

            # normalise the events for normalised plot form 0-100%
            pivot = pd.crosstab(index=df_INDELSUB_summ_sample['sample'],\
                        columns=df_INDELSUB_summ_sample['Event_around_cut'], \
                        values=df_INDELSUB_summ_sample['value'],\
                        aggfunc='sum', normalize='index')
            df_INDELSUB_samples_norm = pivot.reset_index()
            df_INDELSUB_samples_norm = pd.melt(df_INDELSUB_samples_norm, id_vars=['sample'], value_vars=['IN', 'DEL', 'SUB', 'OTHERS'])
            title = 'INDELSUB analysis normalised'
            name_to_save = f'01_03_{i}_INDELSUB_horizontal_barchart_all_samples_normalised'
            text_bellow_chart = f'INDELSUBOTHERS category clasified only for events around the cut.\
                                \nOTHERS are combinations of INs, DELs and SUBs (INDELs etc.).\
                                \nIN+DEL+SUB+OTHERS = 100%.'
            plot = plot_general_stacked_horizontal_barchart(df=df_INDELSUB_samples_norm, x='Event_around_cut', y='sample', values='value', project_output_dir=project_output_subdir, data_source=data_source, title=title, name_to_save=name_to_save, output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False, palette=INDELSUB_color_dic, x_as_perc=True, text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default', xaxis_formater=1, decimals=0, xlim=(0,1), bars_order=['IN', 'DEL', 'SUB', 'OTHERS'], show_plots=show_plots)
            stat_info = df_INDELSUB_samples_norm.groupby('Event_around_cut').value.describe().unstack()
            save_stat(stat_info, name_to_save, project_output_subdir)
            i+=1
    ##########################################################################################
    #            I   N    N
    #            I   N N  N
    #            I   N  N N
    #            I   N   NN
    #            I   N    N
    ##########################################################################################
    if plot_INs is True:
        # Get data
        df_IN = df.copy()
        pivot = pd.pivot_table(df_IN, index=['name', 'n_inserted', '4th_nt'], columns=(), values='%Reads', aggfunc='sum', fill_value=0)
        # Find largest number in n_inserted and make in each sample all numbers from 0 to max present. If not present fill with 0.
        max_n_inserted = df_IN['n_inserted'].max()
        all_n_inserted_list = list(range(0, max_n_inserted+1))
        index = pd.MultiIndex.from_product([df_IN['name'].unique(), all_n_inserted_list, df_IN['4th_nt'].unique()], names=['name', 'n_inserted', '4th_nt'])
        pivot = pivot.reindex(index, fill_value=0)
        df_pivot = pivot.reset_index()
        df_pivot = df_pivot[df_pivot['n_inserted']!=0]
        ##########################################################################################
        # SUM up IN>IN_greater_than nt to single column, to change the number, change grater_than
        ##########################################################################################
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
        name_to_save=f'02_01_IN_no_of_IN_nucleotides_summary_violinplot'
        plot = plot_sns_violin_and_print(data=df_pivot, x='n_inserted', y='value', number_of_columns=len(df_pivot['n_inserted'].unique()), show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, order=None, project_output_dir=project_output_subdir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=title, xtickt_rotation=0, ylabel=ylabel, yaxis_formater=1, xlabel="No of inserted nucleotides", name_to_save=name_to_save, text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=1.02, output_format=output_format, show_plots=show_plots)
        stat_info = df_pivot.groupby('n_inserted').value.describe().unstack()
        save_stat(stat_info, name_to_save, project_output_subdir)

        title = 'Number od inserted nucleotides'
        xlabel='Normalised percentage of inserted nucleotides'
        name_to_save=f'02_02_IN_no_of_IN_nucleotides_stacked_barchart_name'
        # palette_dic = {0:'#ffffff', 1:'#A8B7FF', 2:'#536CE8', 3:'#424763', 'greater_than_2':'#6F79A8', 'greater_than_3':'#6F79A8'}
        plot = plot_general_stacked_horizontal_barchart(df=df_pivot, x='n_inserted', y='name', project_output_dir=project_output_subdir, values='value', sns_palette=False, palette=color_dic, data_source=data_source, title=title, name_to_save=name_to_save, output_format=output_format, dpi=200, xlabel=xlabel, x_as_perc=True, xaxis_formater=1, decimals=0, xlim=(0,1), text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default', y_order=None, bars_order=None, show_plots=show_plots)
        #########################################################################################
        #########################################################################################
        #########################################################################################
        # Plot if the 4th nucleotide the inserted nucleotide
        # Get data
        df_IN = df.copy()
        df_IN = df_IN[df_IN['Event_around_cut'] == 'IN']
        df_IN = df_IN[df_IN['n_inserted'] == 1]
        pivot = pd.pivot_table(df_IN, index=['sample', '4th_nt'], columns='4th_nt_inserted', values='%Reads', aggfunc='sum', fill_value=0)
        # check any sample is present, if not, make dummy sample with 0 values to be able to plot
        # check if columns True and False are present, if not make them with 0 values
        list_to_check = [True, False]
        for item in list_to_check:
            if item not in pivot.columns.tolist():
                pivot[item] = int()
        if pivot.shape[0] == 0:
            pivot = pd.DataFrame({'sample':['dummy_sample'], '4th_nt':['A'], True:[0], False:[0]})
            pivot = pivot.set_index(['sample', '4th_nt'])
        df_pivot = pivot.reset_index()
        # for each row make two rows and make one column named '4th_nt_inserted' and fill there 'False' or 'True' and in the last column named '%Reads' fill the value from the original row where '4th_nt_inserted' is False or True
        df_pivot = pd.melt(df_pivot, id_vars=['sample', '4th_nt'], value_vars=[True, False], var_name='4th_nt_inserted', value_name='%Reads')
        ##########################################################################################
        # plot STACKED HORIZONTAL BARCHART for 4th nt inserted True/False
        ##########################################################################################
        # If more then 40 samples are present, divide the plot into as many plots as needed to have maximum 40 samples in one plot, if not, make one plot with all samples
        # count df_pivot['sample'].unique() and if more than 40, make a list of lists with maximum 40 unique samples in each list
        unique_samples = df_pivot['sample'].unique()
        if len(unique_samples) > 40:
            sample_lists = [unique_samples[i:i+40] for i in range(0, len(unique_samples), 40)]
        else:
            sample_lists = [unique_samples]
        i=1
        for sample_list in sample_lists:
            df_pivot_sample = df_pivot[df_pivot['sample'].isin(sample_list)]
            color_dic={True:'#06145E', False:'#B8BEDB'}
            text_bellow_chart = 'Only 1nt IN around the cut are plotted.\nSequences without DEL/SUB around the cut.'
            title='Is the 4th nucleotide the inserted nucleotide?'
            name_to_save=f'02_03_{i}_4th_nt_inserted_not_normalised'
            xlabel='Sum of %Reads from individual samples'
            plot = plot_general_stacked_horizontal_barchart(df=df_pivot_sample, y='sample'\
                , x='4th_nt_inserted', project_output_dir=project_output_subdir\
                , values='%Reads', sns_palette=False, palette=color_dic\
                , hatches_dic=None\
                , data_source=data_source\
                , title=title\
                , name_to_save=name_to_save, output_format=output_format\
                , dpi=200, xlabel=xlabel, x_as_perc=True, xaxis_formater=100\
                , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default'\
                , bars_order=[True, False], xlim=(0,100), show_plots=show_plots)
            stat_info = df_pivot_sample.groupby('4th_nt_inserted').describe().unstack()
            print(stat_info)
            save_stat(stat_info, name_to_save, project_output_subdir)
            ##########################################################################################
            # Normalise to have True+False = 100%
            pivot_sample = pd.crosstab(index=[df_pivot_sample['sample'], df_pivot_sample['4th_nt']],\
                        columns=df_pivot_sample['4th_nt_inserted'], \
                        values=df_pivot_sample['%Reads'],\
                        aggfunc='sum', normalize='index')
            df_pivot_sample = pivot_sample.reset_index()

            list_to_check = [True, False]
            for item in list_to_check:
                if item not in df_pivot_sample.columns.tolist():
                        df_pivot_sample[item] = int()

            df_pivot_sample = pd.melt(df_pivot_sample, id_vars=(['sample', '4th_nt']), value_vars=[True, False])

            title='Is the 4th nucleotide the inserted nucleotide? Normalised'
            name_to_save=f'02_04_{i}_4th_nt_inserted_normalised'
            plot = plot_general_stacked_horizontal_barchart(df=df_pivot_sample, y='sample'\
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
            stat_info = df_pivot_sample.groupby('4th_nt_inserted').describe().unstack()
            save_stat(stat_info, name_to_save, project_output_subdir)
            i+=1
        # Normalise to have True+False = 100%
        pivot = pd.crosstab(index=[df_pivot['sample'], df_pivot['4th_nt']],\
                    columns=df_pivot['4th_nt_inserted'], \
                    values=df_pivot['%Reads'],\
                    aggfunc='sum', normalize='index')
        df_pivot = pivot.reset_index()
        list_to_check = [True, False]
        for item in list_to_check:
            if item not in df_pivot.columns.tolist():
                    df_pivot[item] = int()
                    
        df_pivot = pd.melt(df_pivot, id_vars=(['sample', '4th_nt']), value_vars=[True, False])
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
        title='Is the 4th nucleotide the inserted nucleotide?'
        no_of_samples = len(df_pivot['sample'].unique())
        ylabel='Normalised percentage\n of inserted 4th nucleotides'
        text_bellow_chart = f'Only 1nt IN around the cut are plotted.\
                            \nSequences without DEL/SUB around the cut.\
                            \nIn individual sample True + False = 100%.\
                            \n{no_of_gRNA_samples}'
        name_to_save = f'02_05_IN_4th_nt_inserted_normalised_violinplot'
        plot = plot_sns_violin_and_print(data=df_pivot, x='4th_nt_inserted', y='value', number_of_columns=1.2\
            , show_datapoints=True, jitter=0.2, stripplot_offset=25, dodge=True\
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
        ##########################################################################################
        # plot VIOLINPLOT 2
        ##########################################################################################
        name_to_save = f'02_06_IN_4th_nt_inserted_normalised_violinplot'
        text_bellow_chart = f'Only 1nt IN around the cut are plotted.\
                    \nSequences without DEL/SUB around the cut.\
                    \nIn individual sample True + False = 100%.\
                    \nHorizontal lines show medians for individual 4th nucleotide.\
                    \n{no_of_gRNA_samples}'
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

    ##########################################################################################
    #               DDDD    EEEEEE  L
    #               D   D   E       L
    #               D   D   EEEE    L
    #               D   D   E       L
    #               DDDD    EEEEEE  LLLLLL
    ##########################################################################################
    file_number = '03'
    counter = 1
    if plot_DELs is True:
        print('Plotting DELs')
        # boxplot from all individual samples showing DEL orientation
        df_orientation = df[['order', 'gRNA_seq', 'sample', 'name', 'DEL_orientation'\
                    , '%Reads', 'n_deleted_around_cut', 'Gap_around_cut', 'Event_around_cut']]
        # drop data without gap
        if 'X' in df_orientation['DEL_orientation'].unique():
            df_orientation = df_orientation[~df_orientation['DEL_orientation'].str.contains('X')]
        # keep oly DELs Event_around_cut
        df_orientation = df_orientation[df_orientation['Event_around_cut']=='DEL']
        # filter out deletions which are not around the cut
        df_orientation = df_orientation[df_orientation['Gap_around_cut']==True]
        # Change ['U', 'B', 'D', 'C'] in DEL_orientation for ['upstream', 'bidirectional', 'downstream', 'cannot be determined']
        # Replace values in the DEL_orientation column
        df_orientation['DEL_orientation'] = df_orientation['DEL_orientation'].replace({
            'U': 'upstream',
            'B': 'bidirectional',
            'D': 'downstream',
            'C': 'cannot be determined'
        })
        order = ['upstream', 'bidirectional', 'downstream', 'cannot be determined']
        ##########################################################################################
        # plot VIOLINPLOT
        ##########################################################################################
        # make pivot without normalisation
        pivot = pd.crosstab(index=[df_orientation['name'], df_orientation['gRNA_seq']], \
                columns=df_orientation['DEL_orientation'], values=df_orientation['%Reads'],\
                aggfunc='sum')
        # reset index and NaN values to 0
        df_pivot = pivot.fillna(0).reset_index()
        ##########################################################################################
        # If df_pivot is Empty skip the next plot
        if df_pivot.empty:
            write_into_report('There is no DEL around the cut, skipping Not-normalised and normalised incidence of DEL according to the orientation plotting.', destination_folder=project_output_dir)
        else:
            ##########################################################################################
            # plot VIOLINPLOT not normalised
            ##########################################################################################
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
            name_to_save = f'03_01_DEL_orientation_summary_violinplot'
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
            ##########################################################################################
            # plot VIOLINPLOT normalised
            ##########################################################################################
            # normalise data as percentage from row total
            pivot = pd.crosstab(index=[df_orientation['name'], df_orientation['gRNA_seq']]\
                , columns=df_orientation['DEL_orientation'], values=df_orientation['%Reads']\
                , aggfunc='sum', normalize='index')
            df_pivot = pivot.reset_index()

            # check if all columns exists, if not make an empty one
            order = ['upstream', 'bidirectional', 'downstream', 'cannot be determined']
            for orientation in order:
                if orientation not in df_pivot.columns.tolist():
                    df_pivot[orientation] = int()

            columns = df_pivot.iloc[:0,2:]
            df_pivot = df_pivot.melt(id_vars=['name', 'gRNA_seq'], \
                value_vars=columns)

            # Set charts properties
            title = 'Normalised incidence of DEL according to the orientation'
            ylabel = 'Normalised percentage of DEL\naround the cut'
            no_of_samples = df_pivot.groupby('DEL_orientation').value.count().unique()
            text_bellow_chart = f'Only DELs around the cut are plotted.\
                                \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                                \n{no_of_gRNA_samples}'
            name_to_save = f'03_02_DEL_orientation_summary_violinplot'
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
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        # Get data
        df_CutType = df.copy()
        ## drop data without gap
        df_CutType = df_CutType[df_CutType['n_deleted_around_cut'] != 0]
        ## filter out deletions which are not around the cut
        df_CutType = df_CutType[df_CutType['Gap_around_cut'] == True]
        # Keep only Event_around_cut for DEL
        # df_CutType = df_CutType[df_CutType['Event_around_cut'] == 'DEL']
        # Make new column MH_around_cut, '0' is False, 'gap_out_of_cut' is 'gap span analysed window' and all other values are True
        df_CutType['MH_around_cut'] = df_CutType['Maximal_MH_len'].apply(lambda x: 'gap spans analysed window' if x == 'gap_out_of_cut' else ('False' if x == '0' else 'True'))
        df_CutType['MH_around_cut'].unique()
        # make pivot, where all samples have individual row and values are sum of individual events
        pivot = pd.pivot_table(df_CutType, index=['name', 'gRNA_seq'], columns=['Cut_type', 'MH_around_cut'], values='%Reads', aggfunc='sum', fill_value=0)
        # If df_pivot is Empty skip the next plot
        if pivot.empty:
            write_into_report('There is no DEL for blunt/stagg, skipping Blunt/Staggered cut ratio plotting.', destination_folder=project_output_dir)
        else:
            df_Blunt_Stagg_samples = pivot.reset_index()
            # Flatten the index
            df_Blunt_Stagg_samples.columns = ['_'.join(map(str, col)).strip('_') for col in df_Blunt_Stagg_samples.columns]
            ## check if all events are present, if not, add missing column with 0
            df_Blunt_Stagg_samples, list_of_events = check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, 'BluntStag', cut_maximum, cut_minimum)
            # make single column table
            df_Blunt_Stagg_samples = pd.melt(df_Blunt_Stagg_samples, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
            pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq']], columns=df_Blunt_Stagg_samples['variable'], values=df_Blunt_Stagg_samples['value'], aggfunc='sum', normalize='index')
            df_Blunt_Stagg_samples_norm = pivot.reset_index()
            df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
            # In column variable rename values
            df_Blunt_Stagg_samples_norm['variable'] = df_Blunt_Stagg_samples_norm['variable'].replace({
                '-2_False': '-2nt staggered cut',
                '-2_True': '-2nt staggered cut with complementarity',
                '-1_False': '-1nt staggered cut',
                '-1_True': '-1nt staggered cut with complementarity',
                '0_False': 'blunt cut',
                '0_True': 'blunt cut with complementarity',
                'C_False': 'type of cut cannot be determined',
                'C_True': 'type of cut cannot be determined with complementarity',
                'X_gap spans analysed window': 'DEL span analysed window'
            })
            # Change list_of_events accordingly
            list_of_events = ['-2nt staggered cut', '-2nt staggered cut with complementarity', '-1nt staggered cut', '-1nt staggered cut with complementarity', 'blunt cut', 'blunt cut with complementarity', 'type of cut cannot be determined', 'type of cut cannot be determined with complementarity', 'DEL span analysed window']
            # check if all events are present in list_of_events, if not add based on the df_Blunt_Stagg_samples_norm['variable'].unique()
            for event in df_Blunt_Stagg_samples_norm['variable'].unique():
                if event not in list_of_events:
                    list_of_events.append(event)
            ##########################################################################################
            # Horizontal bar charts for Blunt/Staggered cuts
            ##########################################################################################
            color_dic = make_color_dic_CutType(list_of_events)
            hatches_dic = make_hatches_dic_CutType(list_of_events)
            if plot_DEL_individual_gRNAs is True:
                # If more then 40 samples are present, divide the plot into as many plots as needed to have maximum 40 samples in one plot, if not, make one plot with all samples
                unique_samples = df_Blunt_Stagg_samples_norm['name'].unique()
                if len(unique_samples) > 400:
                    # shorten the list to unique samples with more than 400 samples, to be able to make plots, if more than 400 samples are present, make a list of lists with maximum 400 unique samples in each list
                    write_into_report(f'There are {len(unique_samples)} samples, which is more than 400, only first 400 are plotted.', destination_folder=project_output_dir)
                    unique_samples = unique_samples[:400]
                if len(unique_samples) > 40:
                    sample_lists = [unique_samples[i:i+40] for i in range(0, len(unique_samples), 40)]
                else:
                    sample_lists = [unique_samples]
                i=1
                for sample_list in sample_lists:
                    df_Blunt_Stagg_samples_norm_sample = df_Blunt_Stagg_samples_norm[df_Blunt_Stagg_samples_norm['name'].isin(sample_list)]
                    # normalise the events for normalised plot form 0-100%
                    title = 'Staggered/Blunt cut ratio for DEL adjacent or around the cut site'
                    xlabel = 'Average of %Reads from individual samples'
                    name_to_save = f'03_03_{i}_Stagg_blunt_ratio_all_samples'
                    text_bellow_chart = 'Only DELs around the cut are plotted.\
                                    \nNormalised to: SUM of all DEL categories around the cut in individual sample = 100%.'

                    plot = plot_general_stacked_horizontal_barchart(df=df_Blunt_Stagg_samples_norm_sample, x='variable'\
                    , y='name', values='value', project_output_dir=project_output_subdir\
                    , data_source=data_source, title=title, name_to_save=name_to_save\
                    , output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False\
                    , palette=color_dic, hatches_dic=hatches_dic, x_as_perc=True\
                    , xaxis_formater=1, decimals=0, xlim=(0,1)\
                    , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default'\
                    , bars_order=list_of_events, show_plots=show_plots)
                    stat_info = df_Blunt_Stagg_samples_norm_sample.groupby('variable').value.describe()
                    save_stat(stat_info, name_to_save, project_output_subdir)
                    i+=1
            ##########################################################################################
            # plot VIOLINPLOT without PAM downstream/upstream orientation
            ##########################################################################################
            title='Normalised Staggered/Blunt representation'
            ylabel='Normalised percentage of DELs \naround cuts events'
            text_bellow_chart = f'Only DELs around the cut are plotted.\
                                \nNormalised to: SUM of all DEL categories around the cut\n in individual sample = 100%.\
                                \n{no_of_gRNA_samples}'
            name_to_save = f'03_04_DEL_Stagg_blunt_summary_violinplot'
            plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='variable', y='value'\
                , number_of_columns=len(df_Blunt_Stagg_samples_norm['variable'].unique())/2\
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
                , y_legend_coordinate=-0.85, yaxis_formater=1\
                , x_max=None\
                , change_edge_color=False\
                , name_to_save=name_to_save\
                , project_output_dir=project_output_subdir\
                , data_source=data_source, output_format=output_format\
                , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
            stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)
        ##########################################################################################
        # plot VIOLINPLOTs with PAM downstream/upstream orientation
        ##########################################################################################
        # make pivot, where all samples have individual row and values are sum of individual events
        pivot = pd.pivot_table(df_CutType, index=['name', 'gRNA_seq'], columns=['DEL_orientation', 'Cut_type', 'MH_around_cut'], values='%Reads', aggfunc='sum', fill_value=0)
        # make order for pivot table based on the upstream/downstream, cut maximum/minimum and presence of MH around the cut
        # If pivot is Empty skip the next plot
        if pivot.empty:
            write_into_report('There is no DEL for PAM downstream/upstream orientation, skipping plotting.', destination_folder=project_output_dir)
        else:
            order = list(itertools.product(['U', 'D'], range(cut_maximum, cut_minimum+1), ['False', 'True']))
            order.append(('C', 0, 'True'))
            order.append(('B', 'C', 'False'))
            order.append(('B', 'C', 'True'))
            order.append(('X', 'X', 'gap spans analysed window'))
            order = [(item[0], str(item[1]), item[2]) for item in order]
            pivot = pivot.reindex(columns=pd.MultiIndex.from_tuples(order), fill_value=0)
            # Flatten the order
            order = ['_'.join(map(str, col)).strip('_') for col in pivot.columns]	
            # Reset the index
            df_Blunt_Stagg_samples = pivot.reset_index()
            # Flatten the index
            df_Blunt_Stagg_samples.columns = ['_'.join(map(str, col)).strip('_') for col in df_Blunt_Stagg_samples.columns]
            # # make single column table
            df_Blunt_Stagg_samples = pd.melt(df_Blunt_Stagg_samples, id_vars=['name', 'gRNA_seq'], value_vars=order)
            ###############################################################################
            # plot VIOLINPLOT with PAM downstream/upstream orientation with cannot_be_determined NOT-normalised
            ##########################################################################################
            # Make pivot crosstab without normalisation
            pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq']]\
                        , columns=df_Blunt_Stagg_samples['variable']\
                        , values=df_Blunt_Stagg_samples['value']\
                        , aggfunc='sum')
            df_Blunt_Stagg_samples_norm = pivot.reset_index()
            df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq'], value_vars=order)
            # In column variable change 'U_-2_0' to 2nt upstream staggered cut, 'U_-2_1' to 2nt upstream staggered cut with complementarity, 'U_-1_0' to 1nt upstream staggered cut, 'U_-1_1' to 1nt upstream staggered cut with complementarity, 'U_0_0' to blunt cut, 'U_0_1' to blunt cut with complementarity, 'C_0_0' to type of cut cannot be determined, 'B_C_0' to bidirectional
            df_Blunt_Stagg_samples['variable'] = df_Blunt_Stagg_samples_norm['variable'].replace({
                'U_-2_False': 'upstream -2nt staggered cut',
                'U_-2_True': 'upstream -2nt staggered cut with complementarity',
                'U_-1_False': 'upstream -1nt staggered cut',
                'U_-1_True': 'upstream -1nt staggered cut with complementarity',
                'U_0_False': 'upstream blunt cut',
                'U_0_True': 'upstream blunt cut with complementarity',
                'D_-2_False': 'downstream -2nt staggered cut',
                'D_-2_True': 'downstream -2nt staggered cut with complementarity',
                'D_-1_False': 'downstream -1nt staggered cut',
                'D_-1_True': 'downstream -1nt staggered cut with complementarity',
                'D_0_False': 'downstream blunt cut',
                'D_0_True': 'downstream blunt cut with complementarity',
                'C_0_True': 'orientation cannot be determined blunt cut with complementarity',
                'B_C_False': 'bidirectional',
                'B_C_True': 'bidirectional with complementarity',
                'X_X_gap spans analysed window': 'DEL span analysed window'
            })
            df_Blunt_Stagg_samples_norm['variable'] = df_Blunt_Stagg_samples_norm['variable'].replace({
                'U_-2_False': 'upstream -2nt staggered cut',
                'U_-2_True': 'upstream -2nt staggered cut with complementarity',
                'U_-1_False': 'upstream -1nt staggered cut',
                'U_-1_True': 'upstream -1nt staggered cut with complementarity',
                'U_0_False': 'upstream blunt cut',
                'U_0_True': 'upstream blunt cut with complementarity',
                'D_-2_False': 'downstream -2nt staggered cut',
                'D_-2_True': 'downstream -2nt staggered cut with complementarity',
                'D_-1_False': 'downstream -1nt staggered cut',
                'D_-1_True': 'downstream -1nt staggered cut with complementarity',
                'D_0_False': 'downstream blunt cut',
                'D_0_True': 'downstream blunt cut with complementarity',
                'C_0_True': 'orientation cannot be determined blunt cut with complementarity',
                'B_C_False': 'bidirectional',
                'B_C_True': 'bidirectional with complementarity',
                'X_X_gap spans analysed window': 'DEL span analysed window'
            })
            # Change list_of_events accordingly
            list_of_events = ['upstream -2nt staggered cut', 'upstream -2nt staggered cut with complementarity', 'upstream -1nt staggered cut', 'upstream -1nt staggered cut with complementarity', 'upstream blunt cut', 'upstream blunt cut with complementarity', 'bidirectional', 'bidirectional with complementarity', 'downstream -2nt staggered cut', 'downstream -2nt staggered cut with complementarity', 'downstream -1nt staggered cut', 'downstream -1nt staggered cut with complementarity', 'downstream blunt cut', 'downstream blunt cut with complementarity', 'orientation cannot be determined blunt cut with complementarity', 'DEL span analysed window']
            # check if all events are present in list_of_events, if not add based on the df_Blunt_Stagg_samples_norm['variable'].unique()
            for event in df_Blunt_Stagg_samples_norm['variable'].unique():
                if event not in list_of_events:
                    list_of_events.append(event)
            # Set charts properties
            title='Not-normalised Staggered/Blunt representation, cut site oriented'
            ylabel='Percentage of DELs \naround cut events'
            text_bellow_chart = f'Only DELs around the cut are plotted.\
                        \nNot-normalised.\
                        \n{no_of_gRNA_samples}'
            name_to_save = '03_05_DEL_Stagg_blunt_oriented_summary_violinplot'
            color_dic = make_color_dic_blunt_stagg_orientation(list_of_events)
            plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples, x='variable', y='value'\
                , number_of_columns=len(list_of_events)/2\
                , show_datapoints=True, jitter=0.10, stripplot_offset=6, dodge=False\
                , show_hue_medians=False\
                , hue=None, hue_only_for_data_points=False \
                , hatches=True\
                , clip_right=True\
                , title=title, xtickt_rotation=-55, ylabel=ylabel\
                , order=list_of_events\
                , plot_legend=False\
                , integrated_properties_dic=color_dic, text_bellow_chart=text_bellow_chart\
                , y_legend_coordinate=-0.80, yaxis_formater=100\
                , x_max=None\
                , plot_vertical_line_middle=False, change_edge_color=True\
                , name_to_save=name_to_save\
                , project_output_dir=project_output_subdir\
                , data_source=data_source, output_format=output_format\
                , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
            stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)
            ##########################################################################################
            # plot VIOLINPLOT with PAM downstream/upstream orientation with cannot_be_determined
            ##########################################################################################
            # Make pivot crosstab where values will be normalised to sum of row = 100%
            pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq']]\
                        , columns=df_Blunt_Stagg_samples['variable']\
                        , values=df_Blunt_Stagg_samples['value']\
                        , aggfunc='sum', normalize='index')
            df_Blunt_Stagg_samples_norm = pivot.reset_index()
            df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq']\
                                        , value_vars=list_of_events)
            # Set charts properties
            title='Normalised Staggered/Blunt representation, cut site oriented'
            ylabel='Normalised percentage of DELs \naround cut events'
            text_bellow_chart = f'Only DELs around the cut are plotted.\
                                \nNormalised to: SUM of all DEL categories around the cut\n in individual sample = 100%.\
                                \n{no_of_gRNA_samples}'
            name_to_save = '03_06_DEL_Stagg_blunt_oriented_summary_violinplot_normalised'
            plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='variable'\
                , y='value', number_of_columns=len(list_of_events)/2\
                , show_datapoints=True, jitter=0.10, stripplot_offset=6, dodge=False\
                , show_hue_medians=False\
                , hue=None, hue_only_for_data_points=False \
                , hatches=True, clip_right=True\
                , title=title, xtickt_rotation=-55, ylabel=ylabel\
                , order=list_of_events, plot_legend=False\
                , integrated_properties_dic=color_dic, text_bellow_chart=text_bellow_chart\
                , y_legend_coordinate=-0.80, yaxis_formater=1\
                , x_max=None\
                , plot_vertical_line_middle=False, change_edge_color=True\
                , name_to_save=name_to_save\
                , project_output_dir=project_output_subdir\
                , data_source=data_source, output_format=output_format\
                , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
            stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        # make DataFrame for ramining DEL_orientations
        # take columns name, n_deleted_around_cut, Gap_around_cut, Event_around_cut, DEL_orientation, Original_Sequence_around_cut, Maximal_MH_len, Cut_type, sample, %Reads, #Reads
        df_CutType = df[['name', 'n_deleted_around_cut', 'Gap_around_cut', 'Event_around_cut', 'DEL_orientation', 'Original_Sequence_around_cut', 'Maximal_MH_len', 'Cut_type', 'sample', '%Reads', '#Reads']].copy()
        # drop data without gap
        df_CutType = df_CutType[df_CutType['n_deleted_around_cut'] != 0]
        # filter out deletions which are not around the cut
        df_CutType = df_CutType[df_CutType['Gap_around_cut'] == True]
        # filter out all events which are not DEL
        df_CutType = df_CutType[df_CutType['Event_around_cut'] == 'DEL']
        # filter our bidirectional DELs
        df_CutType = df_CutType[df_CutType['DEL_orientation'] != 'B']
        df_CutType = df_CutType[df_CutType['DEL_orientation'] != 'X']
        # make DataFrame where samples where DEL_orientation cannot_be_determined is splitted between downstream and upstream
        df_CutType_cannot_be_determined = df_CutType[df_CutType['DEL_orientation']=='C'].copy()
        df_CutType_cannot_be_determined['%Reads'] = df_CutType_cannot_be_determined['%Reads']/2
        df_CutType_cannot_be_determined['#Reads'] = df_CutType_cannot_be_determined['#Reads']/2
        df_CutType_cannot_be_determined_upstream = df_CutType_cannot_be_determined.copy()
        df_CutType_cannot_be_determined_upstream['DEL_orientation'] = 'sU'
        df_CutType_cannot_be_determined_downstream = df_CutType_cannot_be_determined.copy()
        df_CutType_cannot_be_determined_downstream['DEL_orientation'] = 'sD'
        df_CutType_cannot_be_determined = pd.concat([df_CutType_cannot_be_determined_downstream,\
                                                        df_CutType_cannot_be_determined_upstream],\
                                                        ignore_index=True).copy()
        # Remove category C as it was splitted between downstream and upstream
        df_CutType = df_CutType[df_CutType['DEL_orientation'] != 'C']

        df_CutType = pd.concat([df_CutType, df_CutType_cannot_be_determined], ignore_index=True)
        # If df_CutType has zero rows, skip the next plotting
        if df_CutType.empty:
            write_into_report('There is no DEL for Dino chart plotting, skipping plotting.', destination_folder=project_output_dir)  
            plot_DEL_dinochart = False
        else:
            # Make new column MH_around_cut, '0' is False, 'gap_out_of_cut' is 'gap span analysed window' and all other values are True
            df_CutType['MH_around_cut'] = df_CutType['Maximal_MH_len'].apply(lambda x: 'gap spans analysed window' if x == 'gap_out_of_cut' else ('False' if x == '0' else 'True'))
            #  make column combining DEL_orientation with Cut_type
            df_CutType['DEL_orientation_cut_type'] = df_CutType.apply(lambda row: row['DEL_orientation'] + '_' + row['Cut_type'] + '_' + str(row['MH_around_cut']), axis=1)
            df_CutType['n_deleted_around_cut'] = df_CutType.apply(lambda row: -int(row['n_deleted_around_cut']) if (row['DEL_orientation']=='U' or row['DEL_orientation']=='sU') else row['n_deleted_around_cut'], axis=1)
            order_list = df_CutType['sample'].unique()
            df_CutType['n_deleted_around_cut'] = df_CutType['n_deleted_around_cut'].astype(int)
        if plot_DEL_dinochart is True:
            subcounter=1
            for sample in order_list:
                print(sample)
                # filter individual gRNA for exact match with sample name
                df_Cut_gRNA = df_CutType[df_CutType['sample'] == sample]
                seq_around = df_Cut_gRNA['Original_Sequence_around_cut'].unique()
                # print(seq_around)
                name_to_save= f'03_07_{subcounter:03d}_Stagg_blunt_DEL_orientation_dino_30nt_NotNormalised_' + sample
                title = 'Stagg/blunt + DEL orientation analysis ' + str(sample) + '\n' + str(seq_around)
                xlabel='Number of deleted nucleotides (minus = upstream from the probable cut, plus = downstream from the probable cut)'
                ylabel='Sum of %Reads'
                text_bellow_chart = 'Data from deletions upstream and downstream from the probable cut only. Without bidirectional deletions.\nCannot be determined DELs were split between upstream and downstream.\nNot normalised.'
                pivot = pd.crosstab(index=[df_Cut_gRNA['DEL_orientation_cut_type']\
                                            , df_Cut_gRNA['n_deleted_around_cut']]\
                        , columns=df_Cut_gRNA['name'], values=df_Cut_gRNA['%Reads']\
                        , aggfunc='sum', normalize=False)
                df_pivot = pivot.reset_index()
                # add column with gRNA name, where will be average from all samples
                df_pivot[sample] = df_pivot.iloc[:,2:].to_numpy().mean(axis=1)

                # Order df_pivot according to the DEL_orientation_cut_type to have blunt first
                df_pivot['Orientation'] = df_pivot['DEL_orientation_cut_type'].str.split('_').str[0]
                df_pivot['Cut_type'] = df_pivot['DEL_orientation_cut_type'].str.split('_').str[1]
                df_pivot['MH_around_cut'] = df_pivot['DEL_orientation_cut_type'].str.split('_').str[2]
                df_pivot = df_pivot.sort_values(by=['Orientation', 'Cut_type', 'MH_around_cut'], ascending=[True, False, False])

                df_pivot['DEL_orientation_cut_type'] = df_pivot['DEL_orientation_cut_type'].replace({
                    'U_-2_False': 'upstream -2nt staggered cut',
                    'U_-2_True': 'upstream -2nt staggered cut with complementarity',
                    'U_-1_False': 'upstream -1nt staggered cut',
                    'U_-1_True': 'upstream -1nt staggered cut with complementarity',
                    'U_0_False': 'upstream blunt cut',
                    'U_0_True': 'upstream blunt cut with complementarity',
                    'D_-2_False': 'downstream -2nt staggered cut',
                    'D_-2_True': 'downstream -2nt staggered cut with complementarity',
                    'D_-1_False': 'downstream -1nt staggered cut',
                    'D_-1_True': 'downstream -1nt staggered cut with complementarity',
                    'D_0_False': 'downstream blunt cut',
                    'D_0_True': 'downstream blunt cut with complementarity',
                    'sU_0_True': 'upstream cannot be determined blunt cut with complementarity',
                    'sD_0_True': 'downstream cannot be determined blunt cut with complementarity',
                })

                # make order bars starting with upstream and descending cut type, then sU and sD, then downstream and descending cut type
                order = ['downstream blunt cut with complementarity',
                        'downstream blunt cut',
                        'downstream -1nt staggered cut with complementarity',
                        'downstream -1nt staggered cut',
                        'downstream -2nt staggered cut with complementarity',
                        'downstream -2nt staggered cut',
                        'downstream cannot be determined blunt cut with complementarity',
                        'upstream cannot be determined blunt cut with complementarity',
                        'upstream blunt cut with complementarity',
                        'upstream blunt cut',
                        'upstream -1nt staggered cut with complementarity',
                        'upstream -1nt staggered cut',
                        'upstream -2nt staggered cut with complementarity',
                        'upstream -2nt staggered cut']
                list_for_order = df_pivot['DEL_orientation_cut_type'].unique().tolist()
                # sort list for order based on the order list above
                list_for_order = sorted(list_for_order, key=lambda x: order.index(x) if x in order else len(order))
                color_dic = make_color_dic_Dino(df_pivot['DEL_orientation_cut_type'].unique())
                hatches_dic = make_hatches_dic_Dino(df_pivot['DEL_orientation_cut_type'].unique())
                plot = plot_general_stacked_vertical_barchart(df=df_pivot, x='n_deleted_around_cut'\
                        , y='DEL_orientation_cut_type', project_output_dir=project_output_subdir\
                        , values=sample, sns_palette=False, palette=color_dic\
                        , hatches_dic=hatches_dic\
                        , data_source=data_source, plot_vertical_line_middle=True\
                        , title=title, name_to_save=name_to_save, output_format=output_format\
                        , dpi=200, xlabel=xlabel, ylabel=ylabel, y_as_perc=True, yaxis_formater=100\
                        , xlim=(-30,30), ylim=None, text_bellow_chart=text_bellow_chart\
                        , fill_empty_x=True, order_y=list_for_order\
                        , y_legend_coordinate=-0.2, show_plots=show_plots)
                subcounter+=1
        else:
            print('Dinocharts are not going to be plotted.')
        ##########################################################################################
        ##########################################################################################
        ##########################################################################################
        # Get data
        df_CutType = df.copy()
        # drop data without gap
        df_CutType = df_CutType[df_CutType['n_deleted_around_cut'] != 0]
        # filter out deletions which are not around the cut
        df_CutType = df_CutType[df_CutType['Gap_around_cut'] == True]
        # filter out all events which are not DEL
        df_CutType = df_CutType[df_CutType['Event_around_cut'] == 'DEL']
        # Make new column MH_around_cut, '0' is False, 'gap_out_of_cut' is 'gap span analysed window' and all other values are True
        df_CutType['MH_around_cut'] = df_CutType['Maximal_MH_len'].apply(lambda x: 'gap spans analysed window' if x == 'gap_out_of_cut' else ('False' if x == '0' else 'True'))
        # make pivot, where all samples have individual row and values are sum of individual events
        pivot = pd.pivot_table(df_CutType, index=['name', 'gRNA_seq', '4th_nt'], columns=['Cut_type', 'MH_around_cut'], values='%Reads', aggfunc='sum', fill_value=0)
        # if pivot is Empty skip the next plot
        if pivot.empty:
            write_into_report('There is no DEL for plotting, skipping Normalised Staggered/Blunt representation plotting.', destination_folder=project_output_dir)
        else:
            df_Blunt_Stagg_samples = pivot.reset_index()
            # Flatten the index
            df_Blunt_Stagg_samples.columns = ['_'.join(map(str, col)).strip('_') for col in df_Blunt_Stagg_samples.columns]
            ## check if all events are present, if not, add missing column with 0
            df_Blunt_Stagg_samples, list_of_events = check_if_all_events_are_present_after_pivot_reset_index(df_Blunt_Stagg_samples, 'BluntStag', cut_maximum, cut_minimum)
            # make single column table
            df_Blunt_Stagg_samples = pd.melt(df_Blunt_Stagg_samples, id_vars=['name', 'gRNA_seq', '4th_nt'], value_vars=list_of_events)
            # rename variable     
            df_Blunt_Stagg_samples['variable'] = df_Blunt_Stagg_samples['variable'].replace({
                '-2_False': '-2nt staggered cut',
                '-2_True': '-2nt staggered cut with complementarity',
                '-1_False': '-1nt staggered cut',
                '-1_True': '-1nt staggered cut with complementarity',
                '0_False': 'blunt cut',
                '0_True': 'blunt cut with complementarity',
                'C_False': 'type of cut cannot be determined',
                'C_True': 'type of cut cannot be determined with complementarity',
                'X_gap spans analysed window': 'DEL span analysed window'
            })
            # Make list of events based on the variable column
            list_of_events = df_Blunt_Stagg_samples['variable'].unique().tolist()
            pivot = pd.crosstab(index=[df_Blunt_Stagg_samples['name'],df_Blunt_Stagg_samples['gRNA_seq'],df_Blunt_Stagg_samples['4th_nt']], columns=df_Blunt_Stagg_samples['variable'], values=df_Blunt_Stagg_samples['value'], aggfunc='sum', normalize='index')
            df_Blunt_Stagg_samples_norm = pivot.reset_index()
            df_Blunt_Stagg_samples_norm = pd.melt(df_Blunt_Stagg_samples_norm, id_vars=['name', 'gRNA_seq', '4th_nt'], value_vars=list_of_events)
            ##########################################################################################
            # plot STACKED BARPLOT
            ##########################################################################################
            color_dic=make_color_dic_CutType(list_of_events)
            hatches_dic = make_hatches_dic_CutType(list_of_events)

            title = 'Normalised Staggered/Blunt representation\
                    \nin samples divided by 4th nucleotide upstream from PAM'
            xlabel = '4th nucleotide upstream from PAM'
            ylabel = 'Normalised percentage of\nDELs around cut events'
            text_bellow_chart = f'Only DELs around the cut are plotted.\
                                \nNormalised to SUM of all DEL categories around the cut in individual sample = 100%,\
                                \n       then samples divided in 4th nt groups and mean was calculated for each DEL categry and plotted.\
                                \n{no_of_gRNA_samples}'
            name_to_save = f'03_08_DEL_Stagg_blunt_summary_barchart_4th_nt'
            plot = plot_general_stacked_vertical_barchart(df=df_Blunt_Stagg_samples_norm, x='4th_nt'\
                , y='variable', values='value', order_y=list_of_events\
                , project_output_dir=project_output_subdir\
                , sns_palette=False, palette=color_dic\
                , hatches_dic=hatches_dic\
                , title=title, xlabel=xlabel, ylabel=ylabel\
                , text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.1\
                , fill_empty_x = False, xlim=None\
                , y_as_perc=True, ylim=(0,1)\
                , name_to_save=name_to_save, data_source=data_source, output_format=output_format, show_plots=show_plots)
            stat_info = df_Blunt_Stagg_samples_norm.groupby('4th_nt').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)
            ##########################################################################################
            # plot VIOLINPLOT
            ##########################################################################################
            title='Normalised Staggered/Blunt representation based on 4th nucleotide upstream from PAM'
            ylabel='Normalised percentage of DELs \naround cuts events'
            text_bellow_chart = f'Only DELs around the cut are plotted.\
                                \nNormalised to: SUM of all DEL categories around the cut in individual sample = 100%.\
                                \n{no_of_gRNA_samples}'
            name_to_save = f'03_09_DEL_Stagg_blunt_summary_violinplot_4th_nt'
            plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='variable', y='value'\
                    , number_of_columns=len(list_of_events)*1.5\
                    , show_datapoints=False, jitter=0.05, stripplot_offset=0, dodge=False\
                    , show_hue_medians=True\
                    , hue='4th_nt', hue_only_for_data_points=False \
                    , hatches=False\
                    , clip_right=False\
                    , title=title, xtickt_rotation=-55, ylabel=ylabel\
                    , order=list_of_events\
                    , plot_legend=True \
                    , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
                    , y_legend_coordinate=-0.75, yaxis_formater=1\
                    , boxplot_distances_with=[0.8, 0.3]\
                    , x_max=None\
                    , name_to_save=name_to_save\
                    , project_output_dir=project_output_subdir \
                    , data_source=data_source, output_format=output_format\
                    , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
            stat_info = df_Blunt_Stagg_samples_norm.groupby('variable').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)

            name_to_save = f'03_10_DEL_Stagg_blunt_summary_violinplot_4th_nt'
            color_dic = make_color_dic_CutType(list_of_events)
            hatches_dic = make_hatches_dic_CutType(list_of_events)
            plot = plot_sns_violin_and_print(data=df_Blunt_Stagg_samples_norm, x='4th_nt', y='value'\
                    , number_of_columns=16\
                    , show_datapoints=False, jitter=0.05, stripplot_offset=0, dodge=False\
                    , show_hue_medians=False\
                    , hue='variable', hue_order=list_of_events, hue_only_for_data_points=False \
                    , hatches=hatches_dic\
                    , clip_right=False\
                    , title=title, xtickt_rotation=0, ylabel=ylabel\
                    , order=None\
                    , plot_legend=True\
                    , integrated_properties_dic=False, color_dic=color_dic, text_bellow_chart=text_bellow_chart\
                    , y_legend_coordinate=-0.08, yaxis_formater=1\
                    , boxplot_distances_with=[0.8, 0.3]\
                    , x_max=None\
                    , name_to_save=name_to_save\
                    , project_output_dir=project_output_subdir \
                    , data_source=data_source, output_format=output_format\
                    , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
            stat_info = df_Blunt_Stagg_samples_norm.groupby('4th_nt').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    if plot_participating_end is True:
        df_end = df.copy()
        # drop data without gap
        df_end = df_end[df_end['n_deleted_around_cut'] != 0]
        # filter out deletions which are not around the cut
        df_end = df_end[df_end['Gap_around_cut'] == True]
        # Keep only DELs
        df_end = df_end[df_end['Event_around_cut'] == 'DEL']
        # filter out deletions which does not have MH around the cut site
        df_end = df_end[df_end['End_participating_on_MH'] != 'X']
        # Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
        df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'either 3\' or 5\''})
        pivot = pd.pivot_table(df_end, index=['name', 'gRNA_seq'], columns=['End_3_5', 'blunt_cut_with_MH_possible'], values='%Reads', aggfunc='sum', fill_value=0)
        # If pitvot is Empty skip the next plot
        if pivot.empty:
            write_into_report('There is no DEL for plotting, skipping Normalised nad Not-normalised ends participating on the repair.', destination_folder=project_output_dir)
        else:
            df_end_samples = pivot.reset_index()
            # Flattern the MultiIndex columns
            df_end_samples.columns = ['_'.join(map(str, col)).strip() for col in df_end_samples.columns.values]
            # Rename columns
            df_end_samples.rename(columns={'name_': 'name', 'gRNA_seq_':'gRNA_seq', '3\'_True':'3\' possible blunt', '5\'_True':'5\' possible blunt', 'either 3\' or 5\'_True':'either 3\' or 5\' possible blunt', '3\'_False':'3\' NOT blunt', '5\'_False':'5\' NOT blunt', 'either 3\' or 5\'_False':'either 3\' or 5\' NOT blunt'}, inplace=True)

            df_end_samples, list_of_events = check_if_all_events_are_present_after_pivot_reset_index(df_end_samples, ['3\' possible blunt', '5\' possible blunt', 'either 3\' or 5\' possible blunt', '3\' NOT blunt', '5\' NOT blunt', 'either 3\' or 5\' NOT blunt'], cut_maximum, cut_minimum)
            df_end_samples = pd.melt(df_end_samples, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
            pivot_normalised = pd.crosstab(index=[df_end_samples['name'],df_end_samples['gRNA_seq']], columns=df_end_samples['variable'], values=df_end_samples['value'], aggfunc='sum', normalize='index')
            pivot_NOTnorm = pd.crosstab(index=[df_end_samples['name'],df_end_samples['gRNA_seq']], columns=df_end_samples['variable'], values=df_end_samples['value'], aggfunc='sum')
            df_end_samples_norm = pivot_normalised.reset_index()
            df_end_samples_norm = pd.melt(df_end_samples_norm, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
            df_end_samples_NOTnorm = pivot_NOTnorm.reset_index()
            df_end_samples_NOTnorm = pd.melt(df_end_samples_NOTnorm, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
            ##########################################################################################
            # Horizontal bar charts for participating ends
            ##########################################################################################
            # NORMALISED
            # If more then 40 samples are present, divide the plot into as many plots as needed to have maximum 40 samples in one plot, if not, make one plot with all samples
            # count df_end_samples_norm['sample'].unique() and if more than 40, make a list of lists with maximum 40 unique samples in each list
            unique_samples = df_end_samples_norm['name'].unique()
            if len(unique_samples) > 400:
                # shorten the list to unique samples with more than 400 samples, to be able to make plots, if more than 400 samples are present, make a list of lists with maximum 400 unique samples in each list
                write_into_report(f'There are {len(unique_samples)} samples, which is more than 400, only first 400 are plotted.', destination_folder=project_output_dir)
                unique_samples = unique_samples[:400]
            if len(unique_samples) > 40:
                sample_lists = [unique_samples[i:i+40] for i in range(0, len(unique_samples), 40)]
            else:
                sample_lists = [unique_samples]
            i=1
            for sample_list in sample_lists:
                df_end_samples_norm_part = df_end_samples_norm[df_end_samples_norm['name'].isin(sample_list)]
                title = 'Normalised ends participating on the repair'
                xlabel = 'Average of %Reads from individual samples'
                text_bellow_chart = f'Only DELs with complementarity are plotted.\
                                    \nNormalised to: SUM of all ends in individual sample = 100%.\
                                    \n{no_of_gRNA_samples}'
                name_to_save = f'04_01_{i}_DEL_end_summary_norm_horizontal_barchart'
                color_dic = get_color_dic_35_violin()
                plot = plot_general_stacked_horizontal_barchart(df=df_end_samples_norm_part, x='variable'\
                    , y='name', values='value', project_output_dir=project_output_subdir\
                    , data_source=data_source, title=title, name_to_save=name_to_save\
                    , output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False\
                    , palette=color_dic, x_as_perc=True\
                    , xaxis_formater=1, decimals=0, xlim=(0,1)\
                    , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default' \
                    , bars_order=['5\' NOT blunt', '5\' possible blunt', 'either 3\' or 5\' NOT blunt', 'either 3\' or 5\' possible blunt', '3\' NOT blunt', '3\' possible blunt'], show_plots=show_plots)
                stat_info = df_end_samples_norm_part.groupby('variable').value.describe()
                save_stat(stat_info, name_to_save, project_output_subdir)

                # NOT-NORMALISED
                df_end_samples_NOTnorm_part = df_end_samples_NOTnorm[df_end_samples_NOTnorm['name'].isin(sample_list)]
                title = 'Ends participating on the repair'
                xlabel = 'Average of %Reads from individual samples'
                text_bellow_chart = f'Only DELs with complementarity are plotted.\
                                    \n{no_of_gRNA_samples}'
                name_to_save = f'04_02_{i}_DEL_end_summary_NOTnorm_horizontal_barchart'
                plot = plot_general_stacked_horizontal_barchart(df=df_end_samples_NOTnorm_part, x='variable'\
                    , y='name', values='value', project_output_dir=project_output_subdir\
                    , data_source=data_source, title=title, name_to_save=name_to_save\
                    , output_format=output_format, dpi=200, xlabel=xlabel, sns_palette=False\
                    , palette=color_dic, x_as_perc=True\
                    , xaxis_formater=100, decimals=0, xlim=(0,100)\
                    , text_bellow_chart=text_bellow_chart, y_legend_coordinate='Default' \
                    , bars_order=['5\' NOT blunt', '5\' possible blunt', 'either 3\' or 5\' NOT blunt', 'either 3\' or 5\' possible blunt', '3\' NOT blunt', '3\' possible blunt'], show_plots=show_plots)
                stat_info = df_end_samples_NOTnorm_part.groupby('variable').value.describe()
                save_stat(stat_info, name_to_save, project_output_subdir)
                i+=1
            #########################################################################################
            # plot VIOLINPLOT
            ##########################################################################################
            # NORMALISED
            # Set charts properties for violinplot
            title = 'Normalised ends participating on the repair'
            ylabel = 'Normalised percentage of ends \n around cut events'
            text_bellow_chart = f'Only DELs with complementarity are plotted.\
                                \nNormalised to: SUM of all ends in individual sample = 100%\
                                \n{no_of_gRNA_samples}'
            name_to_save = f'04_03_DEL_end_summary_norm_violinplot'
            color_dic = get_color_dic_35_violin()
            plot = plot_sns_violin_and_print(data=df_end_samples_norm, x='variable', y='value'\
                , number_of_columns=3\
                , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
                , show_hue_medians=False\
                , hue=None, hue_only_for_data_points=False \
                , hatches=False\
                , clip_right=True\
                , title=title, xtickt_rotation=-55, ylabel=ylabel\
                , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', '3\' possible blunt', 'either 3\' or 5\' NOT blunt', 'either 3\' or 5\' possible blunt']\
                , plot_legend=False\
                , integrated_properties_dic=False, color_dic=color_dic\
                , text_bellow_chart=text_bellow_chart\
                , y_legend_coordinate=-0.45, yaxis_formater=1\
                , x_max=None\
                , change_edge_color=False\
                , name_to_save=name_to_save\
                , project_output_dir=project_output_subdir\
                , data_source=data_source, output_format=output_format\
                , dpi=200, y_min=0, y_max=1.02, show_plots=show_plots)
            stat_info = df_end_samples_norm.groupby('variable').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)

            # NOT- NORMALISED
            # Set charts properties for violinplot
            title = 'Ends participating on the repair'
            ylabel = 'Percentage from all events'
            text_bellow_chart = f'Only DELs with complementarity are plotted.\
                                \n{no_of_gRNA_samples}'
            name_to_save = f'04_04_DEL_end_summary_NOTnorm_violinplot'
            plot = plot_sns_violin_and_print(data=df_end_samples_NOTnorm, x='variable', y='value'\
                , number_of_columns=3\
                , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
                , show_hue_medians=False\
                , hue=None, hue_only_for_data_points=False \
                , hatches=False\
                , clip_right=True\
                , title=title, xtickt_rotation=-55, ylabel=ylabel\
                , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', '3\' possible blunt', 'either 3\' or 5\' NOT blunt', 'either 3\' or 5\' possible blunt']\
                , plot_legend=False\
                , integrated_properties_dic=False, color_dic=color_dic\
                , text_bellow_chart=text_bellow_chart\
                , y_legend_coordinate=-0.45, yaxis_formater=100\
                , x_max=None\
                , change_edge_color=False\
                , name_to_save=name_to_save\
                , project_output_dir=project_output_subdir\
                , data_source=data_source, output_format=output_format\
                , dpi=200, y_min=0, show_plots=show_plots)
            stat_info = df_end_samples_NOTnorm.groupby('variable').value.describe()
            save_stat(stat_info, name_to_save, project_output_subdir)
        ##########################################################################################
        if plot_participating_end_dinochart is True:
            # make DataFrame for ramining DEL_orientations
            df_end = df.copy()
            # drop data without gap
            df_end = df_end[df_end['n_deleted_around_cut'] != 0]
            # filter out deletions which are not around the cut
            df_end = df_end[df_end['Gap_around_cut'] == True]
            # Keep only DELs
            df_end = df_end[df_end['Event_around_cut'] == 'DEL']
            # if df_end has zero rows, skip the next plotting
            if df_end.empty:
                write_into_report('There is no DEL for Dino chart plotting, skipping plotting.', destination_folder=project_output_subdir)  
                # continue
            # make DataFrame where samples where DEL_orientation cannot_be_determined is splitted between downstream and upstream
            df_end_cannot_be_determined = df_end[df_end['DEL_orientation']=='C'].copy()
            df_end_cannot_be_determined['%Reads'] = df_end_cannot_be_determined['%Reads']/2
            df_end_cannot_be_determined['#Reads'] = df_end_cannot_be_determined['#Reads']/2
            df_end_cannot_be_determined_upstream = df_end_cannot_be_determined.copy()
            df_end_cannot_be_determined_upstream['DEL_orientation'] = 'sU'
            df_end_cannot_be_determined_downstream = df_end_cannot_be_determined.copy()
            df_end_cannot_be_determined_downstream['DEL_orientation'] = 'sD'
            df_end_cannot_be_determined = pd.concat([df_end_cannot_be_determined_downstream,\
                                                            df_end_cannot_be_determined_upstream],\
                                                            ignore_index=True).copy()

            # filter out deletions which are on both sites of cut, so downstream and upstream will remain
            df_end = df_end[df_end['DEL_orientation'] != 'B']
            df_end = df_end[df_end['DEL_orientation'] != 'C']
            df_end = df_end[df_end['DEL_orientation'] != 'X']
            df_end = pd.concat([df_end, df_end_cannot_be_determined], ignore_index=True)

            df_end['DEL_orientation_participating_end'] = df_end.apply(lambda row: row['DEL_orientation'] + '_' + str(row['End_3_5']) + '_' + str(row['blunt_cut_with_MH_possible']), axis=1)
            df_end['n_deleted_around_cut'] = df_end.apply(lambda row: -int(row['n_deleted_around_cut']) if (row['DEL_orientation']=='U' or row['DEL_orientation']=='sU') else row['n_deleted_around_cut'], axis=1)
            df_end['n_deleted_around_cut'] = df_end['n_deleted_around_cut'].astype(int)

            df_end['DEL_orientation_participating_end'].unique()
            # rename DEL_orientation_participating_end values to be more descriptive for the plot
            df_end['DEL_orientation_participating_end'] = df_end['DEL_orientation_participating_end'].replace({
                'U_5_True': 'upstream 5\' possible blunt',
                'U_5_False': 'upstream 5\' NOT blunt',
                'U_0_True': 'upstream either 5\' or 3\' possible blunt',
                'U_0_False': 'upstream no MH from predicted cut sites',
                'D_3_True': 'downstream 3\' possible blunt',
                'D_3_False': 'downstream 3\' NOT blunt',
                'D_0_True': 'downstream either 5\' or 3\' possible blunt',
                'D_0_False': 'downstream no MH from predicted cut sites',
                'X_0_unknown': 'Unknown, No MH',
                'sD_0_True': 'downstream cannot be determined blunt cut with complementarity',
                'sU_0_True': 'upstream cannot be determined blunt cut with complementarity'
            })
            subcounter=1
            for sample in df_end['sample'].unique():
                print(sample)
                # filter individual gRNA
                df_end_gRNA = df_end[df_end['sample'] == sample]
                seq_around = df_end_gRNA['Original_Sequence_around_cut'].unique()
                name_to_save= f'04_05_{subcounter:03d}_end_participating_DEL_orientation_dino_30nt_NotNormalised_' + sample
                title = 'End participating on the repair + DEL orientation analysis ' + str(sample) + '\n' + str(seq_around)
                xlabel='Number of deleted nucleotides (minus = upstream from the probable cut, plus = downstream from the probable cut)'
                ylabel='Sum of %Reads'
                text_bellow_chart = 'Data from deletions around cut only.\nWithout bidirectional deletions.\nCannot be determined was split between upstream and downstream.\nNot normalised data.'
                pivot = pd.crosstab(index=[df_end_gRNA['DEL_orientation_participating_end']\
                                            , df_end_gRNA['n_deleted_around_cut']]\
                        , columns=df_end_gRNA['name'], values=df_end_gRNA['%Reads']\
                        , aggfunc='sum', normalize=False)
                df_pivot = pivot.reset_index()
                # add column with gRNA name, where will be average from all samples
                # if more samples share the same name, make average from them
                df_pivot[sample] = df_pivot.iloc[:,2:].to_numpy().mean(axis=1)
                df_pivot = df_pivot.sort_values(by=['DEL_orientation_participating_end'], ascending=[False])
                order = list(df_pivot['DEL_orientation_participating_end'].unique())
                color_dic = make_color_dic_Dino_end(df_pivot['DEL_orientation_participating_end'].unique())
                hatches_dic = make_hatches_dic_Dino_end(df_pivot['DEL_orientation_participating_end'].unique())
                plot = plot_general_stacked_vertical_barchart(df=df_pivot, x='n_deleted_around_cut'\
                        , y='DEL_orientation_participating_end', project_output_dir=project_output_subdir\
                        , values=sample, sns_palette=False, palette=color_dic\
                        , hatches_dic=hatches_dic\
                        , data_source=data_source, plot_vertical_line_middle=True\
                        , title=title, name_to_save=name_to_save, output_format=output_format\
                        , dpi=200, xlabel=xlabel, ylabel=ylabel, y_as_perc=True, yaxis_formater=100\
                        , xlim=(-30,30), ylim=None, text_bellow_chart=text_bellow_chart\
                        , fill_empty_x=True, order_y=order\
                        , y_legend_coordinate=-0.2, show_plots=show_plots)
                subcounter+=1
    else:
        print('Plots for ends participating on the repair are not going to be plotted.')
# Write into the report how long the script run took
end_time = time.time()
total_time = end_time - start_time
write_into_report(f'\n\nScript run time: {total_time/60:.2f} minutes.', destination_folder=project_output_subdir)