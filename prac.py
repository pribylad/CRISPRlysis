# pylint: disable=too-many-lines
# pylint: disable=line-too-long
# pylint: disable=invalid-name
# %%
import time
import numpy as np
import pandas as pd
import seaborn as sns


from .make_mastertables import filter_out_values_from_df, clasify_INDELSUB_from_columns, change_Event_based_on_mut_position, find_ranges, count_no_of_gaps, find_DEL_closest_to_cut, find_min_and_max_gap_ranges, find_DEL_around_cut_range, gap_span_analysed_window, find_microhomology, get_insertion_event, gap_in_aligned_seq, gap_in_range, find_difference, mutation_in_distance, make_original_seq_dic, make_seq_around_optimum_cut, size_of_del_around_the_cut, no_of_gRNA_samples_txt, return_shifted_DEL_range, return_category_of_DEL, return_cut_type_based_on_cut_categories, process_table, create_mastertable, add_columns_to_df, make_subtables_df

from plot_functions import save_stat, make_Ubuntu_mnt_path, make_reverse_complement, find_element_in_list, adjust_box_widths, plot_sns_violin_and_print, plot_general_stacked_horizontal_barchart, plot_barplot, check_if_all_events_are_present_after_pivot_reset_index, plot_general_stacked_vertical_barchart

from dict_functions import make_color_dic_CutType, make_hatches_dic_CutType, make_color_dic_blunt_stagg_orientation, make_color_dic_Dino, make_hatches_dic_Dino, make_color_dic_Dino_end, get_color_dic_35_violin

# Get current date and time
now = datetime.now()
# Format the date to YYMMDD
date = now.strftime('%y%m%d')

# %%
grouped_df_norm = pd.read_csv('random_categories_grouped_norm_-2_50_0_250524.tsv', sep='\t')
grouped_df_norm[['Gap_category_1', 'Gap_category_2', 'Cut_type', 'End_participating_on_MH', 'DEL_orientation', 'End_3_5', 'Likely_cut_type']] = grouped_df['Combination'].str.split('_', expand=True)
grouped_df_norm
