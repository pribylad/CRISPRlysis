# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""
This code will runn all the analysis based on the parameters.txt file which must be located in the same folder as this code CRISPRlysis.py
Output files will be stored in the 'project_output_dir' which is defined in the parameters.txt
Example: 'project_output_dir':'../../CRISPRlysis_outputs/Pribylova_et_al_2022'
"""
# %%
import ast
import csv
from datetime import datetime
from pathlib import Path
import os
from os.path import dirname, split, basename
import sys
import more_itertools as mit
import pandas as pd
from scripts.general_functions import write_into_report, copy_files_with_parent_folder_to_another, find_files_paths_with_criteria, check_if_path_exist_create_if_not

from scripts.make_mastertables import filter_out_values_from_df, clasify_INDELSUB_from_columns, change_Event_based_on_mut_position, find_ranges, count_no_of_gaps, find_DEL_closest_to_cut, find_min_and_max_gap_ranges, find_DEL_around_cut_range, gap_span_analysed_window, find_microhomology, get_insertion_event, gap_in_aligned_seq, gap_in_range, find_difference, mutation_in_distance, make_original_seq_dic, make_seq_around_optimum_cut, size_of_del_around_the_cut, no_of_gRNA_samples_txt, return_shifted_DEL_range, return_category_of_DEL, return_cut_type_based_on_cut_categories, process_table, create_mastertable, add_columns_to_df, make_subtables_df

from scripts.plot_functions import save_stat, make_Ubuntu_mnt_path, make_reverse_complement, find_element_in_list, adjust_box_widths, plot_sns_violin_and_print

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
else:
    mastertable = pd.read_csv(path_to_existing_mastertable_00, sep='\t')
    write_into_report(f'Reading 00_mastertable from: {path_to_existing_mastertable_00}', destination_folder=project_output_dir)

# Filter out values from mastertable
if filter_out_values:
    mastertable = filter_out_values_from_df(mastertable, filter_out_values, project_output_dir)

# Make new mastertable containing additional columns from 00_mastertable or load existing one
if make_mastertable_01 is True:
    df = add_columns_to_df(mastertable, project_output_dir, parameters, date)
else:
    df = pd.read_csv(path_to_existing_mastertable_01, sep="\t")
    write_into_report(f'Reading 01_mastertable from: {path_to_existing_mastertable_01}', destination_folder=project_output_dir)

# %%
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
    print('Here')
    if plot_INDELs is True:
        file_number = '01'
        counter = 1
        # Get data
        df_INDELSUB = df[['order', 'sample', 'gRNA_seq', 'Event_around_cut', '%Reads', 'name']]
        list_of_others = ['INDEL', 'INSUB', 'INDELSUB', 'DELSUB']

        pivot = pd.pivot_table(df_INDELSUB, index=['sample'], columns='Event_around_cut', values='%Reads',\
                aggfunc='sum', fill_value=0)
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
                , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=True, show_hue_medians=False\
                , hue=None, hue_only_for_data_points=True, transparency=0.7 \
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

