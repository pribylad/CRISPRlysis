# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name

"""This code will runn all the analysis based on the parameters.txt file which must be located in the same folder as this code CRISPRlysis.py"""

# %%
import ast
import csv
from datetime import datetime
from pathlib import Path
import os
from os.path import dirname, split, basename
import more_itertools as mit
import pandas as pd
from functions.general_functions import write_into_report, copy_files_with_parent_folder_to_another, find_files_paths_with_criteria, check_if_path_exist_create_if_not

from functions.make_mastertables import filter_out_values_from_df, clasify_INDELSUB_from_columns, change_Event_based_on_mut_position, find_ranges, count_no_of_gaps, find_DEL_closest_to_cut, find_min_and_max_gap_ranges, find_DEL_around_cut_range, gap_span_analysed_window, find_microhomology, get_insertion_event, gap_in_aligned_seq, gap_in_range, find_difference, mutation_in_distance, make_original_seq_dic, make_seq_around_optimum_cut, size_of_del_around_the_cut, no_of_gRNA_samples_txt, return_shifted_DEL_range, return_category_of_DEL, return_cut_type_based_on_cut_categories, process_table, create_mastertable, add_columns_to_df, make_subtables_df


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

# Make a report.txt file
write_into_report(f'Report from {data_source_original} from date {date}.')
write_into_report(f'Report from {data_source_original} from date {date}.')
write_into_report(f'Project output dir is {project_output_dir}.')
write_into_report(f'Cut position is set to location: {cut_position}')

# Check if project output dir exists, create it if not
check_if_path_exist_create_if_not(project_output_dir)

# Make new mastertable from CRISPResso2 outputs or load existing one
if make_mastertable_00 is True:
    create_mastertable(source_dir, project_output_dir, parameters)
    # create_mastertable does not return the DataFrame, thus the table must be loaded later
    mastertable = pd.read_csv(Path(str(project_output_dir) + '/00_mastertable_' + date + '.tsv'), sep='\t')
else:
    mastertable = pd.read_csv(path_to_existing_mastertable_00, sep='\t')
    write_into_report(f'Reading 00_mastertable from: {path_to_existing_mastertable_00}')

# Filter out values from mastertable
if filter_out_values is True:
    mastertable = filter_out_values_from_df(mastertable, filter_out_values, project_output_dir)

# Make new mastertable containing additional columns from 00_mastertable or load existing one
if make_mastertable_01 is True:
    df = add_columns_to_df(mastertable, project_output_dir, parameters, date)
else:
    df = pd.read_csv(path_to_existing_mastertable_01, sep="\t")
    write_into_report(f'Reading 01_mastertable from: {path_to_existing_mastertable_01}')

if make_subtables is True:
    tables_dic = make_subtables_df(df, project_output_dir, subtables_param, date)
# %%

