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
import scipy.stats as stats
import seaborn as sns
import statannotations
from statannotations.Annotator import Annotator

from scripts.plot_functions import plot_sns_violin_and_print, save_stat
from scripts.dict_functions import get_color_dic_35_violin
##########################################################################################
output_format = 'png'
show_plots = True
plot_4th_nucleotide_chart = False
x = 95  # efficacy threshold percentage

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
# filter_criterion = 'HCT116'
# filter_value = '48'
# df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit{100-x}efficacy_{filter_criterion}_RNP_{filter_value}_{datetime.now().strftime("%y%m%d")}.tsv'
# data_source = f'{filter_criterion}_RNP_{filter_value}'
# project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/time_series/{filter_criterion}_{filter_value}')

filter_criterion = 'Kamoen'
filter_value = 'WT'
df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_95/03_filtered_efficacy_samples_limit5efficacy_{filter_criterion}_{filter_value}_{datetime.now().strftime("%y%m%d")}.tsv'
data_source = f'{filter_criterion}_{filter_value}'
project_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/time_series/{filter_criterion}_{filter_value}')

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
##################################################################################
#       Ends participating on the repair
##################################################################################
# %%
file_number = '003_Supp'
counter = 2

# Define list of indexes for pivot table
default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism', 'replicate']
list_of_events = ['3\' NOT blunt', '5\' possible blunt', '5\' NOT blunt', 'eider 3\' or 5\' possible blunt', 'eider 3\' or 5\' NOT blunt']
df_end = df[default_list + ['Gap_around_cut', 'Maximal_MH_len', 'End_3_5', 'blunt_cut_with_MH_possible', '%Reads']]
# Remove rows where Maximal_MH_len is 'no_gap', '0' or 'gap_out_of_cut'
df_end = df_end[~df_end['Maximal_MH_len'].isin(['no_gap', '0', 'gap_out_of_cut'])]
# filter out deletions which are not around the cut
df_end = df_end[df_end['Gap_around_cut'] == True]
# Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'eider 3\' or 5\''})
# Change True and False in  blunt_cut_with_MH_possible to 'possible blunt' and 'NOT possible blunt'
df_end['blunt_cut_with_MH_possible'] = df_end['blunt_cut_with_MH_possible'].replace({'TRUE':'possible blunt', 'FALSE':'NOT blunt', True:'possible blunt', False:'NOT blunt', 'True':'possible blunt', 'False':'NOT blunt'})

# Make column where both End_3_5 and blunt_cut_with_MH_possible are combined
df_end['End_3_5_blunt'] = df_end['End_3_5'] + ' ' + df_end['blunt_cut_with_MH_possible']
# remove columns End_3_5, blunt_cut_with_MH_possible
df_end = df_end.drop(columns=['End_3_5', 'blunt_cut_with_MH_possible'], inplace=False)

# Get data
pivot = pd.pivot_table(df_end, index=default_list, columns='End_3_5_blunt', values='%Reads', aggfunc='sum', fill_value=0)
# reset pivot and remove index
df_pivot = pivot.reset_index()
# remove index name
df_pivot.columns.name = None
# Normalise rows to 100%
df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
df_pivot.head()
# Melt dataframe
df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
# remove 'replicate' from default_list
default_list.remove('replicate')
# Make pivot table with all columns without column replicate and make average of values
pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
# Melt pivot table_avg reset index
pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
pivot_table_avg
#########################################################################################
# plot VIOLINPLOT
##########################################################################################
# NORMALISED
# Set charts properties for violinplot
title = 'Normalised ends participating on the repair, mean of replicates'
ylabel = 'Normalised percentage of ends \n around cut events'
no_of_samples =  len(pivot_table_avg[['source', 'gRNA_seq', 'organism']].drop_duplicates())
text_bellow_chart = f'Only DELs around the cut are plotted.\
                    \nNormalised to: SUM of all ends in individual sample = 100%\
                    \n3\' possible blunt can not be in -2,0,0 setup -> it is eider 3\' or 5\' NOT blunt\
                    \nNo of samples: {no_of_samples}'

name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_violinplot'
color_dic = get_color_dic_35_violin()
# TODO split it to True/False for can be blunt cut with MH
plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
    , number_of_columns=3\
    , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
    , show_hue_medians=False\
    , hue=None, hue_only_for_data_points=False \
    , hatches=False\
    , clip_right=True\
    , title=title, xtickt_rotation=-55, ylabel=ylabel\
    , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']\
    , plot_legend=False\
    , integrated_properties_dic=False, color_dic=color_dic\
    , text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.45, yaxis_formater=100\
    , x_max=None\
    , change_edge_color=False\
    , name_to_save=name_to_save\
    , project_output_dir=project_output_dir\
    , data_source=data_source, output_format=output_format\
    , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
# Save pivot table as tsv
pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
# Save statistics
stat_info = pivot_table_avg.groupby('variable').value.describe()
save_stat(stat_info, name_to_save, project_output_dir)
##########################################################################################
# Statistical testing between 5' NOT blunt and 3' NOT blunt
##########################################################################################
file_number = '003_Supp'
counter = 1

# Count statistics for 5' NOT blunt and 3' NOT blunt, those data are compositional, so must be transformed before statistical testing as ALR=log(x1/x2)
df_stats = pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])]
# Pivot the dataframe to have separate columns for each variable
df_stats = df_stats.pivot_table(index=default_list, columns='variable', values='value').reset_index()
# Calculate ALR transformation log(x1/x2) and avoid division by zero by adding a small constant 1e-12
df_stats['ALR_5_NOTblunt_vs_3_NOTblunt'] = np.log((df_stats["5' NOT blunt"] + 1e-12) / (df_stats["3' NOT blunt"] + 1e-12))

# Perform Wilcoxon signed-rank test on ALR values
statistic, p_value = stats.wilcoxon(df_stats['ALR_5_NOTblunt_vs_3_NOTblunt'])
print(f'Wilcoxon signed-rank test statistic: {statistic}, p-value: {p_value}')
median_lf = np.median(df_stats['ALR_5_NOTblunt_vs_3_NOTblunt'])
if median_lf > 0:
    interpretation = "5’ microhomology staggered cut is significantly more frequent than 3’ microhomology staggered cut."
elif median_lf < 0:
    interpretation = "3’ microhomology staggered cut is significantly more frequent than 5’ microhomology staggered cut."
else:
    interpretation = "No difference in frequency between 5’ microhomology staggered cut and 3’ microhomology staggered cut."
print(f'Median is {median_lf}. Interpretation based on median log-fold change: {interpretation}')
# Save statistics
stat_info = pd.DataFrame({'statistic': [statistic], 'p_value': [p_value]})
save_stat(stat_info, f'{file_number}_{counter}_DEL_end_summary_norm_violinplot_filtered_stats', project_output_dir)

plt.figure(figsize=(7, 4)) 
plt.hist(df_stats['ALR_5_NOTblunt_vs_3_NOTblunt'], bins=40, color="steelblue", edgecolor="black", alpha=0.7)
plt.axvline(0, color="gray", linestyle="--", label="5’ = 3’")
plt.axvline(median_lf, color="green", linestyle="--", label="Median log-ratio")
plt.xlabel("log(5’ microhomology staggered cut/3’ microhomology staggered cut)")
plt.ylabel("Count of samples")
plt.title(f"{filter_value}\nDistribution of log-ratios: 5’ microhomology staggered cut vs 3’ microhomology staggered cut")
plt.legend()
plt.tight_layout()
plt.savefig(Path(str(project_output_dir) + f'/{file_number}_{counter}_DEL_end_summary_stats_histogram_' + date + f'.{output_format}'), dpi=200)
plt.show()
# %%
# Plot only 5' NOT blunt and 3' NOT blunt from the pivot_table_avg as x variables
file_number = '003_Fig'
counter = 1
name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_violinplot_filtered'
color_dic = get_color_dic_35_violin()
text_bellow_chart = f'Only DELs around the cut are plotted.\
                    \nNormalised to: SUM of all DELs with MH in individual sample = 100%\
                    \nInterpretation based on median log-fold change:\n    {interpretation}\
                    \n    Wilcoxon signed-rank test statistic: {statistic}, p-value: {p_value}\
                    \nNo of samples: {no_of_samples}'

plot = plot_sns_violin_and_print(data=pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])], x='variable', y='value'\
    , number_of_columns=2\
    , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
    , show_hue_medians=False\
    , hue=None, hue_only_for_data_points=False \
    , hatches=False\
    , clip_right=True\
    , title=title, xtickt_rotation=-55, ylabel=ylabel\
    , order=['5\' NOT blunt', '3\' NOT blunt']\
    , plot_legend=False\
    , integrated_properties_dic=False, color_dic=color_dic\
    , text_bellow_chart=text_bellow_chart\
    , y_legend_coordinate=-0.25, yaxis_formater=100\
    , x_max=None\
    , change_edge_color=False\
    # , add_stat=True, stat_pairs=[("5' NOT blunt", "3' NOT blunt")], stat_test='t-test_ind', stat_format='full', stat_correction=None\
    , name_to_save=name_to_save\
    , project_output_dir=project_output_dir\
    , data_source=data_source, output_format=output_format\
    , dpi=200, y_min=0, y_max=119, show_plots=show_plots)

# Plot End participating on the repair with respect of the 4th nucleotid
if plot_4th_nucleotide_chart is True:
    file_number = '003_Supp'
    counter = 3
    # Define list of indexes for pivot table
    default_list = ['source', 'gene', 'gRNA_seq', '4th_nt', 'mutant', 'replicate']
    list_of_events = ['3\' NOT blunt', '5\' possible blunt', '5\' NOT blunt', 'eider 3\' or 5\' possible blunt', 'eider 3\' or 5\' NOT blunt']
    df_end = df[default_list + ['Gap_around_cut', 'Maximal_MH_len', 'End_3_5', 'blunt_cut_with_MH_possible', '%Reads']]
    # Remove rows where Maximal_MH_len is 'no_gap', '0' or 'gap_out_of_cut'
    df_end = df_end[~df_end['Maximal_MH_len'].isin(['no_gap', '0', 'gap_out_of_cut'])]
    # filter out deletions which are not around the cut
    df_end = df_end[df_end['Gap_around_cut'] == True]
    # Change values in End_3_5 to be 3', 5' and cannot_be_determined needed for plotting with color_dic
    df_end['End_3_5'] = df_end['End_3_5'].replace({3:'3\'', 5:'5\'', 0:'eider 3\' or 5\''})
    # Change True and False in  blunt_cut_with_MH_possible to 'possible blunt' and 'NOT possible blunt'
    df_end['blunt_cut_with_MH_possible'] = df_end['blunt_cut_with_MH_possible'].replace({'TRUE':'possible blunt', 'FALSE':'NOT blunt'})

    # Make column where both End_3_5 and blunt_cut_with_MH_possible are combined
    df_end['End_3_5_blunt'] = df_end['End_3_5'] + ' ' + df_end['blunt_cut_with_MH_possible']
    # remove columns End_3_5, blunt_cut_with_MH_possible
    df_end = df_end.drop(columns=['End_3_5', 'blunt_cut_with_MH_possible'], inplace=False)

    # Get data
    pivot = pd.pivot_table(df_end, index=default_list, columns='End_3_5_blunt', values='%Reads', aggfunc='sum', fill_value=0)
    # reset pivot and remove index
    df_pivot = pivot.reset_index()
    # remove index name
    df_pivot.columns.name = None
    # Normalise rows to 100%
    df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
    df_pivot.head()

    # Melt dataframe
    df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
    # remove 'replicate' from default_list
    default_list.remove('replicate')
    # Make pivot table with all columns without column replicate and make average of values
    pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
    # Melt pivot table_avg reset index
    pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
    pivot_table_avg

    #########################################################################################
    # plot VIOLINPLOT
    ##########################################################################################
    # NORMALISED 4th as hue
    # Set charts properties for violinplot
    title = 'Normalised ends participating on the repair, mean of replicates'
    ylabel = 'Normalised percentage of ends \n around cut events'
    text_bellow_chart = f'Only DELs around the cut are plotted.\
                        \nNormalised to: SUM of all ends in individual sample = 100%\
                        \n3\' possible blunt can not be in -2,0,0 setup -> it is eider 3\' or 5\' NOT blunt\
                        \nNo of samples: {no_of_samples}'
    name_to_save = f'{file_number}_{counter}_DEL_end_summary_norm_4thnt_violinplot'
    color_dic = get_color_dic_35_violin()
    # TODO split it to True/False for can be blunt cut with MH
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
        , number_of_columns=15\
        , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=True\
        , show_hue_medians=False\
        , boxplot_distances_with = [0.8, 0.18]\
        , hue='4th_nt', hue_only_for_data_points=False \
        , hatches=False\
        , clip_right=True\
        , title=title, xtickt_rotation=-55, ylabel=ylabel\
        , order=['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']\
        , hue_order=['A', 'C','G','T']\
        , plot_legend=True\
        , integrated_properties_dic=True, color_dic=None\
        , text_bellow_chart=text_bellow_chart\
        , y_legend_coordinate=-0.45, yaxis_formater=100\
        , x_max=None\
        , change_edge_color=False\
        , name_to_save=name_to_save\
        , project_output_dir=project_output_dir\
        , data_source=data_source, output_format=output_format\
        , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
    # Save pivot table as tsv
    pivot_table_avg.to_csv(Path(str(project_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
    # Save statistics
    stat_info = pivot_table_avg.groupby('variable').value.describe()
    save_stat(stat_info, name_to_save, project_output_dir)

# %%
