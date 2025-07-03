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

from scripts.plot_functions import save_stat, make_Ubuntu_mnt_path, make_reverse_complement, find_element_in_list, adjust_box_widths, plot_sns_violin_and_print, plot_general_stacked_horizontal_barchart, plot_barplot, check_if_all_events_are_present_after_pivot_reset_index, plot_general_stacked_vertical_barchart

from scripts.dict_functions import make_color_dic_CutType, make_hatches_dic_CutType, make_color_dic_blunt_stagg_orientation, make_color_dic_Dino, make_hatches_dic_Dino, make_color_dic_Dino_end
# %%
cut_maximum = -2
cut_minimum = 0
# %%
df = pd.read_csv('/home/lefnui/CRISPRlysis_outputs/Pribylova_et_al_2022/02_mastertable_P35S_WT_250702.tsv', sep='\t')
df
# %%
df_end = df
# drop data without gap
df_end = df_end[df_end['n_deleted_around_cut'] != 0]
# filter out deletions which are not around the cut
df_end = df_end[df_end['Gap_around_cut'] == True]
# filter out deletionw which does not have MH around the cut site
df_end = df_end[df_end['End_participating_on_MH'] != 'X']
# Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'eider 3\' or 5\''})
pivot = pd.pivot_table(df_end, index=['name', 'gRNA_seq'], columns=['End_3_5', 'blunt_cut_with_MH_possible'], values='%Reads', aggfunc='sum', fill_value=0)
df_end_samples = pivot.reset_index()
# Flattern the MultiIndex columns
df_end_samples.columns = ['_'.join(map(str, col)).strip() for col in df_end_samples.columns.values]
# Rename columns
df_end_samples.rename(columns={'name_': 'name', 'gRNA_seq_':'gRNA_seq', '3\'_True':'3\' possible blunt', '5\'_True':'5\' possible blunt', 'eider 3\' or 5\'_True':'eider 3\' or 5\' possible blunt', '3\'_False':'3\' NOT blunt', '5\'_False':'5\' NOT blunt', 'eider 3\' or 5\'_False':'eider 3\' or 5\' NOT blunt'}, inplace=True)

df_end_samples, list_of_events = check_if_all_events_are_present_after_pivot_reset_index(df_end_samples, ['3\' possible blunt', '5\' possible blunt', 'eider 3\' or 5\' possible blunt', '3\' NOT blunt', '5\' NOT blunt', 'eider 3\' or 5\' NOT blunt'], cut_maximum, cut_minimum)
df_end_samples = pd.melt(df_end_samples, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
pivot_normalised = pd.crosstab(index=[df_end_samples['name'],df_end_samples['gRNA_seq']], columns=df_end_samples['variable'], values=df_end_samples['value'], aggfunc='sum', normalize='index')
pivot_normalised
pivot_NOTnorm = pd.crosstab(index=[df_end_samples['name'],df_end_samples['gRNA_seq']], columns=df_end_samples['variable'], values=df_end_samples['value'], aggfunc='sum')
df_end_samples_norm = pivot_normalised.reset_index()
df_end_samples_norm = pd.melt(df_end_samples_norm, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
df_end_samples_NOTnorm = pivot_NOTnorm.reset_index()
df_end_samples_NOTnorm = pd.melt(df_end_samples_NOTnorm, id_vars=['name', 'gRNA_seq'], value_vars=list_of_events)
df_end_samples_norm
df_end_samples_NOTnorm