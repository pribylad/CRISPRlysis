# pylint: disable=line-too-long
# pylint: disable=invalid-name

# This script is for analysis of DELs around the cut site, their orientation and ends participating on the repair.

# %%
# import necessary libraries
from datetime import datetime
from pathlib import Path
import matplotlib
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import seaborn as sns
import statannotations
from statannotations.Annotator import Annotator
import statsmodels

from scripts.plot_functions import plot_sns_violin_and_print, save_stat
from scripts.dict_functions import get_color_dic_35_violin
##########################################################################################
output_format = 'png'
# output_format = 'svg'
show_plots = False
individual_dataset = True    # If False (for settings 1, 2, and 3), Kingdom and ID_no_reps is used in the default_list and additional plots for them are made.
use_3_5_possible = False

# # 1 no filter, analysis on all samples after merging 02_mastertable*.tsv tables from CRISPRlysis outputs
# filter_criterion = 'None'
# filter_value = 'None'
# data_source = 'all_samples_no_filters'
# individual_dataset = False
# no_of_runs = 1
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_All_260409.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_all_samples')

# # 2 with data filtered for overall efficacy >=5%, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy'
# filter_value = '95'
# individual_dataset = False
# no_of_runs = 1
# data_source = 'all_samples_filtered_95'
# # df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
# # charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/01_all_samples')
# # table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/01_all_samples')
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
# data_source = 'all_samples_filtered_95'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/01_all_samples')

# # 3 with data filtered for overall efficacy >=5% and kingdom = plant, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy_and_kingdom'
# filter_value = '95_and_plant'
# individual_dataset = True
# no_of_runs = 1
# data_source = 'plants_filtered_95'
# # df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/04_split_kingdom_plant.tsv'
# # charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/02_plants')
# # table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/02_plants')
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/04_split_kingdom_plant.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/02_plants')

# # 4 with data filtered for overall efficacy >=5% and kingdom = animal, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy_and_kingdom'
# filter_value = '95_and_animal'
# individual_dataset = True
# no_of_runs = 1
# data_source = 'animals_filtered_95'
# # df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/04_split_kingdom_animal.tsv'
# # charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/02_animals')
# # table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/02_animals')
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/04_split_kingdom_animal.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/02_animals')


# filter_criterion = 'Random_dataset'
# use_3_5_possible = True
# # filter_value = '-2_50_0'
# # filter_value = '-2_50_2'
# filter_value = '0_50_2'
# individual_dataset = True
# data_source = filter_criterion + '_' + filter_value
# df_dir = f'/home/lefnui/CRISPRlysis_outputs/random_2025/251020_PAM/01_mastertable_PAM_{filter_value}_251020.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/04_{filter_criterion}_{filter_value}_random_PAM')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/04_{filter_criterion}_{filter_value}_random_PAM')

# df_dir = f'/home/lefnui/CRISPRlysis_outputs/random_2025/251019_noPAM/01_mastertable_251019_{filter_value}_noPAM.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/04_{filter_criterion}_{filter_value}_random_noPAM')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/04_{filter_criterion}_{filter_value}_random_noPAM')

# 5 with data filtered for overall efficacy >=5% and individual dataset, input table is result from 02_filter_low_efficacy.py script
# Open list_of_datasets.txt and read it as list of datasets to be analysed individually
# list_of_datasets_file = '/home/lefnui/git/CRISPRlysis/list_of_individual_samples_1.txt'
# list_of_datasets_file = '/home/lefnui/git/CRISPRlysis/list_of_individual_samples_2.txt'
# mutants = True
# list_of_datasets_file = '/home/lefnui/git/CRISPRlysis/list_of_individual_samples_3_mutants.txt'
list_of_datasets_file = '/home/lefnui/git/CRISPRlysis/list_of_individual_samples_4_with_Ponce.txt'
list_of_datasets = []
with open(list_of_datasets_file, 'r') as f:
    for line in f:
        dataset = line.strip()
        if dataset:  # Check if the line is not empty
            list_of_datasets.append(dataset)
print(list_of_datasets)
no_of_runs = len(list_of_datasets)
print(f"No of individual datasets to be analysed: {no_of_runs}")

# # 6 with data from Cristina Ponce Lilli unpublished, withour filtering for efficacy, input table is result from 01_merge_filter_02_tables.py script
# filter_criterion = 'efficacy_and_author'
# filter_value = '100_and_Ponce'
# no_of_runs = 1
# data_source = 'Ponce_unfiltered'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_Ponce_260423.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples/05_Ponce')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/05_Ponce')

for run in range(no_of_runs):
    if no_of_runs > 1:
        if mutants is False:
            individual_dataset = True
            filter_criterion = 'efficacy_and_individual_dataset'
            filter_value = list_of_datasets[run]
            # df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/{filter_value}'
            df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/{filter_value}'
            # Remove 04_split_index_ from data_source and also .tsv, replace " _ " with "_" and spaces with "_" in data_source
            data_source = f'{filter_value}'
            data_source = data_source.replace('split_', '').replace('index_', '').replace('.tsv', '').replace(' - ', '_').replace(' ', '_').replace('A.', 'A').replace('.,','_').replace(',', '').replace('__', '_')
            # charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/{data_source}')
            # table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/{data_source}')
            charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/{data_source}')
            table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/{data_source}')
        else:
            filter_criterion = 'mutant'
            filter_value = list_of_datasets[run]
            df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/mutants/{filter_value}'
            # Remove 04_split_index_ from data_source and also .tsv, replace " _ " with "_" and spaces with "_" in data_source
            data_source = f'{filter_value}'
            data_source = data_source.replace('split_', '').replace('index_', '').replace('.tsv', '').replace(' - ', '_').replace(' ', '_').replace('A.', 'A').replace('.,','_').replace(',', '').replace('__', '_')
            print(f"Data source: {data_source}")
            charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/03_mutants/{data_source}')
            table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/mutants/{data_source}')
    print(f"Run {run+1}/{no_of_runs} sample being analysed: {filter_value}")
    ##########################################################################################
    now = datetime.now()
    # Format the date to YYMMDD
    date = now.strftime('%y%m%d')

    if charts_output_dir.exists() is False:
        os.makedirs(charts_output_dir)

    if table_output_dir.exists() is False:
        os.makedirs(table_output_dir)
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
    df.columns
    
    df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)
    df['blunt_cut_with_MH_possible'] = df['blunt_cut_with_MH_possible'].astype(str)
    # change TRUE to True and FALSE to False in blunt_cut_with_MH_possible column
    df['blunt_cut_with_MH_possible'] = df['blunt_cut_with_MH_possible'].replace({'TRUE':'True', 'FALSE':'False'})
    print(f"Loaded dataframe with shape {df.shape} from {df_dir}")
    
    df['blunt_cut_with_MH_possible'].unique()
    
    df['Maximal_MH_len'].unique()
    
    ##########################################################################################
    #               DDDD    EEEEEE  L
    #               D   D   E       L
    #               D   D   EEEE    L
    #               D   D   E       L
    #               DDDD    EEEEEE  LLLLLL
    ##########################################################################################

    # Define list of indexes for pivot table
    # default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism', 'replicate']
    if individual_dataset == False:
        default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep', 'gRNA_num', 'replicate']
    else:
        default_list = ['source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation','replicate']
    # list_of_events = ['upstream', 'both', 'downstream', 'cannot_be_determined', 'X']
    list_of_events = ['upstream_ambiguous', 'both_ambiguous', 'downstream_ambiguous', 'upstream_NO_ambi', 'both_NO_ambi', 'downstream_NO_ambi', 'cannot_be_determined_ambiguous', 'cannot_be_determined_NO_ambi', 'X_gap_out_of_cut']
    df_DEL_orientation = df[default_list + ['DEL_orientation', 'Gap_around_cut', 'Maximal_MH_len', '%Reads']]
    # Filter out rows where Gap_around_cut is False
    df_DEL_orientation = df_DEL_orientation[df_DEL_orientation['Gap_around_cut'] == True]
    # Replace values in the DEL_orientation column
    df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'].replace({
        'U': 'upstream',
        'B': 'both',
        'D': 'downstream',
        'C': 'cannot_be_determined'
    })
    order = ['upstream_ambiguous', 'both_ambiguous', 'downstream_ambiguous', 'upstream_NO_ambi', 'both_NO_ambi', 'downstream_NO_ambi', 'cannot_be_determined_ambiguous', 'cannot_be_determined_NO_ambi', 'X_gap_out_of_cut']
    print(df_DEL_orientation['Maximal_MH_len'].unique())
    # Replace in Maximal_MH_len '0' for 'NO_ambi', 'gap_out_of_cut' keep and anything else as 'ambiguous'
    df_DEL_orientation['Maximal_MH_len'] = df_DEL_orientation['Maximal_MH_len'].replace({'0': 'NO_ambi', 'gap_out_of_cut': 'gap_out_of_cut'})
    df_DEL_orientation['Maximal_MH_len'] = df_DEL_orientation['Maximal_MH_len'].apply(lambda x: 'ambiguous' if x not in ['NO_ambi', 'gap_out_of_cut'] else x)
    # Replace column 'DEL_orientation' with merged column of DEL_orientation and Maximal_MH_len
    df_DEL_orientation['DEL_orientation'] = df_DEL_orientation['DEL_orientation'] + '_' + df_DEL_orientation['Maximal_MH_len']
    print(df_DEL_orientation['DEL_orientation'].unique())
    
    df_DEL_orientation
    # Remove columns Gap_around_cut and Maximal_MH_len as they are not needed anymore
    df_DEL_orientation = df_DEL_orientation.drop(columns=['Gap_around_cut', 'Maximal_MH_len'], inplace=False)
    df_DEL_orientation
    
    ##########################################################################################
    #        Deletion orientation around the cut site
    ##########################################################################################
    # Get data
    pivot = pd.pivot_table(df_DEL_orientation, index=default_list, columns='DEL_orientation', values='%Reads', aggfunc='sum', fill_value=0)
    # reset pivot and remove index
    df_pivot = pivot.reset_index()
    df_pivot
    
    # remove index name
    df_pivot.columns.name = None
    # Save pivot table as tsv
    df_pivot.to_csv(Path(str(table_output_dir) + '/' + '05' + '_DEL_orientation_NOTnormalised_pivot' + '.tsv'), sep='\t')
    # Normalise rows to 100%
    df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
    df_pivot
    
    # Melt dataframe
    # Check if all list_of_events items are in columns, if not create them with 0 values
    for item in list_of_events:
        if item not in df_pivot.columns.tolist():
            df_pivot[item] = 0
    df_pivot = df_pivot.melt(id_vars=default_list, value_vars=list_of_events)
    # Remove 'replicate' from default_list
    default_list.remove('replicate')
    # Make pivot table with all columns without column replicate and make average of values
    pivot_table_avg = create_pivot_table(df_pivot, default_list, columns='variable', values='value', aggfunc='mean')
    # Met pivot table_avg reset index
    pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=default_list, value_vars=list_of_events)
    # Save pivot table as tsv
    pivot_table_avg.to_csv(Path(str(table_output_dir) + '/' + '05_pivot_DEL_orientation_normalised_avg' + '.tsv'), sep='\t')
    # Make column combining organism and cell_line for better visualisation in the plot
    pivot_table_avg['organism_cell_line_hpi'] = pivot_table_avg['organism'] + ' - ' + pivot_table_avg['cell_line'] + ' - ' + pivot_table_avg['hours_post_inoculation']
    pivot_table_avg
    
    # Set charts properties
    title = 'Normalised incidence of DEL according to the orientation,\n mean of replicates'
    ylabel = 'Normalised percentage \nof DEL around the cut'
    no_of_samples =  len(pivot_table_avg) / len(pivot_table_avg['variable'].unique())
    text_bellow_chart = f'Only DELs around the cut are plotted. Without gaps spanning analysed window\
                        \n - see table with source data.\
                        \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                        \nNo of samples: {no_of_samples}'
    name_to_save = f'003_1_DEL_orientation_summary_violinplot'
    # hue order
    hue_order=[
                'Arabidopsis thaliana - xxx - 504'
                ,'Arabidopsis thaliana - root culture - 504'
                ,'Arabidopsis thaliana - xxx - unknown'
                ,'A. thaliana - embryo - unknown'
                ,'A. thaliana - embryo - T1'
                ,'Nicotiana benthamiana - xxx - 120'
                ,'Nicotiana benthamiana - leaf - 120'
                ,'Nicotiana benthamiana - xxx - 336'
                ,'Petunia hybrida - xxx - 24'
                ,'Silene latifolia - xxx - xx'
                ,'Solanum lycopersicum - xxx - 24'
                ,'Solanum lycopersicum - M82 - 48'
                ,'Homo sapiens - HCT116 - 4'
                ,'Homo sapiens - HEK293 - 4'
                ,'Homo sapiens - K562 - 4'
                ,'Homo sapiens - HCT116 - 8'
                ,'Homo sapiens - HEK293 - 8'
                ,'Homo sapiens - K562 - 8'
                ,'Homo sapiens - HCT116 - 16'
                ,'Homo sapiens - HEK293 - 16'
                ,'Homo sapiens - K562 - 16'
                ,'Homo sapiens - HCT116 - 24'
                ,'Homo sapiens - HEK293 - 24'
                ,'Homo sapiens - K562 - 24'
                ,'Homo sapiens - HCT116 - 48'
                ,'Homo sapiens - HEK293 - 48'
                ,'Homo sapiens - K562 - 48'
                ,'Homo sapiens - HAP1 - 48to72'
                ,'Homo sapiens - Hela - 72'
                ,'Homo sapiens - Jurkat - 48to72'
                ,'UniCorn - random_dataset - 0'
                ,'Xenopus tropicalis - xxx - xxx'
                ,'Pleurodeles waltl - xxx - xxx'
                ,'Chlamydomonas reinhardtii - xxx - 48'
    ]
    # Test if all organism_cell_line_hpi.unique() are present in the hue order
    unique_hue_values = pivot_table_avg['organism_cell_line_hpi'].unique()
    missing_hue_values = [value for value in unique_hue_values if value not in hue_order]
    if missing_hue_values:
        raise ValueError(f"Missing values in hue_order: {missing_hue_values}")
    # If there are more values in hue_order and they are missing in organism_cell_line_hpi.unique() remove them from hue_order
    hue_order = [value for value in hue_order if value in unique_hue_values]

    
    if individual_dataset == False:
        # color_dic = {'upstream':'#f0d271', 'both':'#808080', 'downstream':'#cccccc', 'cannot_be_determined':'#ffffff', 'plant':'#009d73', 'animal':'#cc79a6'}
        colour_dic = {'upstream_ambiguous':'#f0d271', 'upstream_NO_ambi':'#f0d271', 'both_ambiguous':'#808080', 'both_NO_ambi':'#808080', 'downstream_ambiguous':'#cccccc', 'downstream_NO_ambi':'#cccccc', 'cannot_be_determined_ambiguous':'#ffffff', 'cannot_be_determined_NO_ambi':'#ffffff', 'X_gap_out_of_cut':'#000000', 'plant':'#009d73', 'animal':'#cc79a6'}
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=15\
            , show_datapoints=True, jitter=0.0, stripplot_offset=0, dodge=True\
            , show_hue_medians=False\
            , hue='organism_cell_line_hpi', hue_only_for_data_points=True \
            , hue_order=hue_order
            , hatches=False\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
            , order=order\
            , plot_legend=True\
            , integrated_properties_dic=True, color_dic=None\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.4, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=True\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='png'\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)

        name_to_save = f'003_2_DEL_orientation_summary_violinplot'
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=6\
            , show_datapoints=True, jitter=0.05, stripplot_offset=0, dodge=True\
            , show_hue_medians=True\
            , hue='kingdom', hue_only_for_data_points=True \
            , hatches=False\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
            , order=order\
            , plot_legend=True\
            , integrated_properties_dic=False, color_dic=colour_dic\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.5, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='png'\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)

    
    # Set charts properties
    title = 'Normalised incidence of DEL according to the orientation,\n mean of replicates'
    ylabel = 'Normalised percentage \nof DEL around the cut'
    no_of_samples =  len(pivot_table_avg) / len(pivot_table_avg['variable'].unique())

    # colour_dic = {'upstream':'#f0d271', 'both':'#808080', 'downstream':'#cccccc', 'cannot_be_determined':'#ffffff'}
    colour_dic = {'upstream_ambiguous':'#f0d271', 'upstream_NO_ambi':'#f0d271', 'both_ambiguous':'#808080', 'both_NO_ambi':'#808080', 'downstream_ambiguous':'#cccccc', 'downstream_NO_ambi':'#cccccc', 'cannot_be_determined_ambiguous':'#ffffff', 'cannot_be_determined_NO_ambi':'#ffffff', 'X_gap_out_of_cut':'#000000', 'plant':'#009d73', 'animal':'#cc79a6'}

    ##########################################################################################
    # Statistical testing between upstream_abiguous and upstream_NO_ambi
    ##########################################################################################

    # Count statistics for upstream_ambiguous and upstream_NO_ambi, those data are compositional, so must be transformed before statistical testing as ALR=log(x1/x2)
    df_stats = pivot_table_avg[pivot_table_avg['variable'].isin(["upstream_ambiguous", "upstream_NO_ambi"])]
    # Pivot the dataframe to have separate columns for each variable
    df_stats = df_stats.pivot_table(index=default_list, columns='variable', values='value').reset_index()
    # export df_stats to tsv
    df_stats.to_csv(Path(str(table_output_dir) + '/' + f'05_DEL_end_summary_norm_violinplot_filtered_stats_{filter_criterion}_{filter_value}' + '_' + date + '.tsv'), sep='\t', index=False)

    # Calculate ALR transformation log(x1/x2) and avoid division by zero by adding a small constant 1e-12
    df_stats['ALR_upstream_ambiguous_vs_upstream_NO_ambi'] = np.log((df_stats["upstream_ambiguous"] + 1e-12) / (df_stats["upstream_NO_ambi"] + 1e-12))

  # Perform Wilcoxon signed-rank test on ALR values
  #   Try if statistic rise ValueError, if yes skip
    try:
        statistic, p_value = stats.wilcoxon(df_stats['ALR_upstream_ambiguous_vs_upstream_NO_ambi'])
        print(f'Wilcoxon signed-rank test statistic: {statistic}, p-value: {p_value}')
        median_lf = np.median(df_stats['ALR_upstream_ambiguous_vs_upstream_NO_ambi'])
        if median_lf > 0:
            interpretation = "Upstream deletions with ambiguos nucleoties are significantly more frequent than upstream deletions without ambiguity."
        elif median_lf < 0:
            interpretation = "Upstream deletions without ambiguity are significantly more frequent than upstream deletions with ambiguity."
        else:
            interpretation = "No difference in frequency between upstream deletions with and without ambiguity."
        print(f'Median is {median_lf}. Interpretation based on median log-fold change: {interpretation}')
        # Save statistics
        stat_info = pd.DataFrame({'statistic': [statistic], 'p_value': [p_value]})
        save_stat(stat_info, f'003_5_DEL_end_summary_norm_violinplot_filtered_stats_{filter_criterion}_{filter_value}', charts_output_dir)
    except ValueError as e:
        print(f"Error occurred: {e}")
        print("Skipping statistical test.")
        statistic, p_value = None, None
    
    plt.figure(figsize=(7, 4)) 
    plt.hist(df_stats['ALR_upstream_ambiguous_vs_upstream_NO_ambi'], bins=40, color="steelblue", edgecolor="black", alpha=0.7)
    plt.axvline(0, color="gray", linestyle="--", label="5’ = 3’")
    plt.axvline(median_lf, color="green", linestyle="--", label="Median log-ratio")
    plt.xlabel("log(upstream_ambiguous/upstream_NO_ambi)")
    plt.ylabel("Count of samples")
    plt.title(f"{filter_value}\nDistribution of log-ratios: upstream_ambiguous vs upstream_NO_ambi")
    plt.legend()
    plt.tight_layout()
    plt.savefig(Path(str(charts_output_dir) + f'/003_3_upstream_DEL_summary_stats_histogram_' + '_' + date + f'.{output_format}'), dpi=200)
    plt.show()

    # Full plot as png
    name_to_save = f'003_3_DEL_orientation_summary_violinplot'
    text_bellow_chart = f'Only DELs around the cut are plotted. Without gaps spanning analysed window\
                    \n - see table with source data.\
                    \nNormalised to: SUM of all DEL orientation in individual sample = 100%\
                    \nInterpretation based on median log-fold change:\n    {interpretation}\
                    \n    Wilcoxon signed-rank test statistic: {statistic}, p-value: {p_value}\
                    \nNo of samples: {no_of_samples}'
    
    try:
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=len(pivot_table_avg['variable'].unique())/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=colour_dic\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.4, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            # , add_stat=True, stat_pairs=[("upstream_ambiguous", "upstream_NO_ambi")], stat_test='Wilcoxon', stat_format='star', stat_correction='BH'\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=115, show_plots=show_plots)

        # Without legend etc as png
        name_to_save = f'003_4_DEL_orientation_summary_violinplot_small'
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=len(pivot_table_avg['variable'].unique())/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=None, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=colour_dic\
            , text_bellow_chart=None\
            , y_legend_coordinate=-0.4, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='png'\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)

        # Without legend etc as svg
        name_to_save = f'003_4_DEL_orientation_summary_violinplot_small'
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=len(pivot_table_avg['variable'].unique())/2\
            , show_datapoints=True, jitter=0.10, stripplot_offset=8, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=None, xtickt_rotation=-55, ylabel=ylabel, xlabel='DEL orientation'\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=colour_dic\
            , text_bellow_chart=None\
            , y_legend_coordinate=-0.4, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='svg'\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)

        stat_info = pivot_table_avg.groupby('variable').value.describe()
        save_stat(stat_info, name_to_save, charts_output_dir)
    except ValueError as e:
        print(f"Error occurred while plotting: {e}")
        print("Skipping plot.")

    ##########################################################################################
    #               DDDD    EEEEEE  L           PPPPP     A    RRRRR   TTTTTTT
    #               D   D   E       L           P    P   A A   R    R     T
    #               D   D   EEEE    L           PPPPP    A A   RRRRR      T
    #               D   D   E       L           P       AAAAA  R   R      T
    #               DDDD    EEEEEE  LLLLLL      p      A     A R    R     T
    ##########################################################################################
    ##################################################################################
    #       Ends participating on the repair
    ##################################################################################

    # Define list of indexes for pivot table
    # default_list = ['source', 'gene', 'gRNA_seq', 'mutant', 'organism', 'replicate']
    if individual_dataset == False:
        default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep', 'replicate']
    else:
        default_list = ['source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'replicate']

    # For random dataset also use 3' possible blunt
    if use_3_5_possible == True:
        list_of_events = ['3\' possible blunt', '3\' NOT blunt', '5\' possible blunt', '5\' NOT blunt', 'eider 3\' or 5\' possible blunt', 'eider 3\' or 5\' NOT blunt']
    else:
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
    df_end
    
    # Get data
    pivot = pd.pivot_table(df_end, index=default_list, columns='End_3_5_blunt', values='%Reads', aggfunc='sum', fill_value=0)
    # reset pivot and remove index
    df_pivot = pivot.reset_index()
    # remove index name
    df_pivot.columns.name = None
    # Normalise rows to 100%
    df_pivot.iloc[:, len(default_list):] = df_pivot.iloc[:, len(default_list):].div(df_pivot.iloc[:, len(default_list):].sum(axis=1), axis=0) * 100
    df_pivot
    
    # If all list_of_events items are not in columns, create them with 0 values
    for item in list_of_events:
        if item not in df_pivot.columns.tolist():
            df_pivot[item] = 0
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
    no_of_samples = len(pivot_table_avg) / len(pivot_table_avg['variable'].unique())
    text_bellow_chart = f'Only DELs around the cut are plotted.\
                        \nNormalised to: SUM of all ends in individual sample = 100%\
                        \n3\' possible blunt can not be in -2,0,0 setup -> it is eider 3\' or 5\' NOT blunt\
                        \nNo of samples: {no_of_samples}'

    color_dic = get_color_dic_35_violin()
    if use_3_5_possible == True:
        order = ['5\' NOT blunt', '5\' possible blunt', '3\' NOT blunt', '3\' possible blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']
    else:
        order = ['5\' NOT blunt','5\' possible blunt', '3\' NOT blunt', 'eider 3\' or 5\' NOT blunt', 'eider 3\' or 5\' possible blunt']

    try:
        # Full plot as png
        name_to_save = f'003_4_DEL_end_summary_norm_violinplot'
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=3\
            , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=title, xtickt_rotation=-55, ylabel=ylabel\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=color_dic\
            , text_bellow_chart=text_bellow_chart\
            , y_legend_coordinate=-0.45, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)

        # Without legend etc as png
        name_to_save = f'003_4_DEL_end_summary_norm_violinplot_small'
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=3\
            , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=None, xtickt_rotation=-55, ylabel=ylabel\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=color_dic\
            , text_bellow_chart=None\
            , y_legend_coordinate=-0.45, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='png'\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)

        # Without legend etc as svg
        plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value'\
            , number_of_columns=3\
            , show_datapoints=True, jitter=0.1, stripplot_offset=10, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=None, xtickt_rotation=-55, ylabel=ylabel\
            , order=order\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=color_dic\
            , text_bellow_chart=None\
            , y_legend_coordinate=-0.45, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='svg'\
            , dpi=200, y_min=0, y_max=102, show_plots=show_plots)
    except ValueError as e:
        print(f"Error occurred: {e}")
        print("Skipping plot.")
    # Save pivot table as tsv
    pivot_table_avg.to_csv(Path(str(table_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
    # Save statistics
    stat_info = pivot_table_avg.groupby('variable').value.describe()
    save_stat(stat_info, name_to_save, charts_output_dir)
    
    ##########################################################################################
    # Statistical testing between 5' NOT blunt and 3' NOT blunt
    ##########################################################################################

    # Count statistics for 5' NOT blunt and 3' NOT blunt, those data are compositional, so must be transformed before statistical testing as ALR=log(x1/x2)
    df_stats = pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])]
    # Pivot the dataframe to have separate columns for each variable
    df_stats = df_stats.pivot_table(index=default_list, columns='variable', values='value').reset_index()
    # export df_stats to tsv
    df_stats.to_csv(Path(str(table_output_dir) + '/' + f'05_DEL_end_summary_norm_violinplot_filtered_stats_{filter_criterion}_{filter_value}' + '_' + date + '.tsv'), sep='\t', index=False)

    # Calculate ALR transformation log(x1/x2) and avoid division by zero by adding a small constant 1e-12
    df_stats['ALR_5_NOTblunt_vs_3_NOTblunt'] = np.log((df_stats["5' NOT blunt"] + 1e-12) / (df_stats["3' NOT blunt"] + 1e-12))


    # Perform Wilcoxon signed-rank test on ALR values
    try:
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
        save_stat(stat_info, f'003_5_DEL_end_summary_norm_violinplot_filtered_stats_{filter_criterion}_{filter_value}', charts_output_dir)
        stat_info
    except ValueError as e:
        print(f"Error occurred: {e}")
        statistic, p_value = None, None

    plt.figure(figsize=(7, 4)) 
    plt.hist(df_stats['ALR_5_NOTblunt_vs_3_NOTblunt'], bins=40, color="steelblue", edgecolor="black", alpha=0.7)
    plt.axvline(0, color="gray", linestyle="--", label="5’ = 3’")
    plt.axvline(median_lf, color="green", linestyle="--", label="Median log-ratio")
    plt.xlabel("log(5’ microhomology staggered cut/3’ microhomology staggered cut)")
    plt.ylabel("Count of samples")
    plt.title(f"{filter_value}\nDistribution of log-ratios: 5’ microhomology staggered cut vs 3’ microhomology staggered cut")
    plt.legend()
    plt.tight_layout()
    plt.savefig(Path(str(charts_output_dir) + f'/003_6_DEL_end_summary_stats_histogram_' + '_' + date + f'.{output_format}'), dpi=200)
    plt.show()
    
    # Plot only 5' NOT blunt and 3' NOT blunt from the pivot_table_avg as x variables
    color_dic = get_color_dic_35_violin()
    text_bellow_chart = f'Only DELs around the cut are plotted.\
                        \nNormalised to: SUM of all DELs with MH in individual sample = 100%\
                        \nInterpretation based on median log-fold change:\n    {interpretation}\
                        \n    Wilcoxon signed-rank test statistic: {statistic}, p-value: {p_value}\
                        \nNo of samples: {no_of_samples}'
    # Full plot as png
    name_to_save = f'003_7_DEL_end_summary_norm_violinplot_filtered'
    try:
        plot = plot_sns_violin_and_print(data=pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])], x='variable', y='value'\
            , number_of_columns=2\
            , show_datapoints=True, jitter=0.05, stripplot_offset=10, dodge=False\
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
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format=output_format\
            , dpi=200, y_min=0, y_max=119, show_plots=show_plots)
    except ValueError as e:
        print(f"Error occurred: {e}")
    # Without legend etc as png
    name_to_save = f'003_7_DEL_end_summary_norm_violinplot_filtered_small'
    try:
        plot = plot_sns_violin_and_print(data=pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])], x='variable', y='value'\
            , number_of_columns=2\
            , show_datapoints=True, jitter=0.05, stripplot_offset=10, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=None, xtickt_rotation=-55, ylabel=ylabel\
            , order=['5\' NOT blunt', '3\' NOT blunt']\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=color_dic\
            , text_bellow_chart=None\
            , y_legend_coordinate=-0.25, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , add_stat=True, stat_pairs=[("5' NOT blunt", "3' NOT blunt")], stat_test='Wilcoxon', stat_format='star', stat_correction='BH'\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='png'\
            , dpi=200, y_min=0, y_max=119, show_plots=show_plots)
    except ValueError as e:
        print(f"Error occurred: {e}")
    # Without legend etc as svg
    try:
        plot = plot_sns_violin_and_print(data=pivot_table_avg[pivot_table_avg['variable'].isin(["5' NOT blunt", "3' NOT blunt"])], x='variable', y='value'\
            , number_of_columns=2\
            , show_datapoints=True, jitter=0.05, stripplot_offset=10, dodge=False\
            , show_hue_medians=False\
            , hue=None, hue_only_for_data_points=False \
            , hatches=False\
            , clip_right=True\
            , title=None, xtickt_rotation=-55, ylabel=ylabel\
            , order=['5\' NOT blunt', '3\' NOT blunt']\
            , plot_legend=False\
            , integrated_properties_dic=False, color_dic=color_dic\
            , text_bellow_chart=None\
            , y_legend_coordinate=-0.25, yaxis_formater=100\
            , x_max=None\
            , change_edge_color=False\
            , add_stat=True, stat_pairs=[("5' NOT blunt", "3' NOT blunt")], stat_test='Wilcoxon', stat_format='star', stat_correction='BH'\
            , name_to_save=name_to_save\
            , project_output_dir=charts_output_dir\
            , data_source=data_source, output_format='svg'\
            , dpi=200, y_min=0, y_max=119, show_plots=show_plots)
    except ValueError as e:
        print(f"Error occurred: {e}")

print('Tadaaaa, all finished :o)')
# %%
