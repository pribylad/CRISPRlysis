# pylint: disable=line-too-long
# pylint: disable=invalid-name

# This script will make Figure 2A -> alanysis of insertion events
# Also, it will filter samples with 1nt insertion >= 80% and save the filtered dataframe for further analysis of these samples in the next script.

# %%
import os
import pandas as pd
from pathlib import Path
from datetime import date, datetime
import matplotlib.pyplot as plt
from scripts.plot_functions import plot_sns_violin_and_print, save_stat
##########################################################################################
output_format = 'png'
show_plots = False
maximum_to_view = 2    # max lenght of insertion to be views as inependent subset, all insertions longer than this will be grouped as >{maximum_to_view}

# # 1 no filter, analysis on all samples after merging 02_mastertable*.tsv tables from CRISPRlysis outputs
# filter_criterion = 'None'
# filter_value = 'None'
# data_source = 'all_samples_no_filters'
# no_of_runs = 1
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_All_260409.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/01_no_filter_all_samples/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_all_samples')

# # 2 with data filtered for overall efficacy >=5%, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy'
# filter_value = '95'
# no_of_runs = 1
# data_source = 'all_samples_filtered_95'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/01_all_samples')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/01_all_samples')

# # 3 with data filtered for overall efficacy >=5% and kingdom = plant, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy_and_kingdom'
# filter_value = '95_and_plant'
# no_of_runs = 1
# data_source = 'plants_filtered_95'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/04_split_kingdom_plant.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/02_plants')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/02_plants')

# # 4 with data filtered for overall efficacy >=5% and kingdom = animal, input table is result from 02_filter_low_efficacy.py script
# filter_criterion = 'efficacy_and_kingdom'
# filter_value = '95_and_animal'
# no_of_runs = 1
# data_source = 'animals_filtered_95'
# df_dir = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/04_split_kingdom_animal.tsv'
# charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/02_animals')
# table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/02_animals')

# 5 with data filtered for overall efficacy >=5% and individual dataset, input table is result from 02_filter_low_efficacy.py script
# Open list_of_datasets.txt and read it as list of datasets to be analysed individually
list_of_datasets_file = '/home/lefnui/git/CRISPRlysis/list_of_individual_samples_1.txt'
# list_of_datasets_file = '/home/lefnui/git/CRISPRlysis/list_of_individual_samples_2.txt'
list_of_datasets = []
with open(list_of_datasets_file, 'r') as f:
    for line in f:
        dataset = line.strip()
        if dataset:  # Check if the line is not empty
            list_of_datasets.append(dataset)
print(list_of_datasets)
no_of_runs = len(list_of_datasets)
print(f"No of individual datasets to be analysed: {no_of_runs}")

for run in range(no_of_runs):
    if no_of_runs > 1:
        filter_criterion = 'efficacy_and_individual_dataset'
        filter_value = list_of_datasets[run]
        df_dir = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/{filter_value}'
        # Remove 04_split_index_ from data_source and also .tsv, replace " _ " with "_" and spaces with "_" in data_source
        data_source = f'{filter_value}'
        data_source = data_source.replace('split_', '').replace('index_', '').replace('.tsv', '').replace(' - ', '_').replace(' ', '_').replace('A.', 'A').replace('.,','_').replace(',', '').replace('__', '_')
        charts_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/{data_source}')
        table_output_dir = Path(f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/{data_source}')
    print(f"Run {run+1}/{no_of_runs} sample being analysed: {filter_value}")
    ##########################################################################################
    now = datetime.now()
    date = now.strftime('%y%m%d')
    # Create output directory if it does not exist
    if not os.path.exists(charts_output_dir):
        os.makedirs(charts_output_dir)
    if not os.path.exists(table_output_dir):
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

    # Make a list of events with numbers from 1 to maximum_to_view and >maximum_to_view
    list_of_events = [str(i) for i in range(1, maximum_to_view+1)]
    list_of_events.append(f'>{maximum_to_view}')
    ##########################################################################################
    # Load datafreame
    df = pd.read_csv(df_dir, sep='\t')
    df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)
    df['blunt_cut_with_MH_possible'] = df['blunt_cut_with_MH_possible'].astype(str)
    # change TRUE to True and FALSE to False in blunt_cut_with_MH_possible column
    df['blunt_cut_with_MH_possible'] = df['blunt_cut_with_MH_possible'].replace({'TRUE':'True', 'FALSE':'False'})
    print(f'{df.shape} -> Original dataframe shape.')
    print(f"{len(df['name'].unique())} -> No of unique samples where each replicate is one sample.")
    print(df[['sample']].drop_duplicates().shape)
    print(len(df['ID_no_rep'].unique()))
    # Define lists
    default_list = ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation', 'ID_no_rep', 'replicate']
    # default_list = ['ID_no_rep', 'replicate']
    df_IN = df[default_list + ['n_inserted', '4th_nt', '4th_nt_inserted', '%Reads']]
    print(f'{df_IN.shape} -> After filtering for needed columns (no of rows should be the same)')

    print(len(df_IN[default_list].drop_duplicates()))

    # Get data
    pivot = pd.pivot_table(df_IN, index=default_list, columns='n_inserted', values='%Reads', aggfunc='sum', fill_value=0)
    print(f'{pivot.shape} -> After making pivot table with n_inserted as columns and %Reads as values.')
    # reset pivot and remove index
    df_pivot = pivot.reset_index()
    # remove index name
    df_pivot.columns.name = None
    # remove column 0
    df_pivot = df_pivot.drop(columns=[0], inplace=False)
    print(f'{df_pivot.shape} -> After dropping column 0')
    # Make sum from the column maximum_to_view+1 to the end as >{maximum_to_view}
    df_pivot[f'>{maximum_to_view}'] = df_pivot.iloc[:, len(default_list)+len(list_of_events):].sum(axis=1)
    # Convert columns to string
    df_pivot.columns = df_pivot.columns.map(str)
    # Check if columns from 1 to maximum_to_view are in the columns of df_pivot, if not add them with value 0
    for i in range(1, maximum_to_view+1):
        if str(i) not in df_pivot.columns:
            df_pivot[str(i)] = 0
    # Keep only default_list + list_of_events columns
    df_pivot = df_pivot[default_list + list_of_events]
    # Normalise rows to 100%
    df_pivot[list_of_events] = df_pivot[list_of_events].div(df_pivot[list_of_events].sum(axis=1), axis=0) * 100
    df_pivot = pd.melt(df_pivot, id_vars=default_list, value_vars=list_of_events)
    print(f'{df_pivot.shape} -> After melting pivot table, when variable column has number of inserted nucleotides and value column has normalised percentage of reads with that number of inserted nucleotides.')
    #####################################################################################
    # MAKE AVERAGE OF REPLICATES
    #####################################################################################
    # Make pivot table with all columns without column replicate and name and make average of values
    pivot_table_avg = create_pivot_table(df_pivot, default_list[:-1], columns='variable', values='value', aggfunc='mean')
    # Met pivot table_avg reset index
    pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars='ID_no_rep', value_vars=list_of_events)
    # Save pivot_table_avg as tsv in table_output_dir with name 03_pivot_table_avg_{date}.tsv
    pivot_table_avg.to_csv(Path(str(table_output_dir) + '/' + f'03_pivot_table_avg_{date}.tsv'), sep='\t', index=False)
    #####################################################################################
    # MAKE VIOLIN PLOT OF IN EVENTS
    #####################################################################################
    color_dic = {'1':'#A8B7FF', '2':'#536CE8', f'>{maximum_to_view}':'#AAAAAA'}
    title='Normalised numbers of inserted nucleotides around the cut site, mean of replicates'
    ylabel='Normalised percentage\n of inserted nucleotides'
    no_of_samples = len(pivot_table_avg) / len(pivot_table_avg['variable'].unique())
    text_bellow_chart = f'Only IN events around cut.\
                        \nSum of all IN in one sample = 100%\
                        \nNo of samples: {no_of_samples}'
    # Full plot as png
    name_to_save=f'002_1_IN_no_of_IN_nucleotides_summary_violinplot'
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value', number_of_columns=len(pivot_table_avg['variable'].unique()), show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, 
    order=None, project_output_dir=charts_output_dir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=title, xtickt_rotation=0, ylabel=ylabel, xlabel='No. of inserted nucleotides', yaxis_formater=100, name_to_save=name_to_save, text_bellow_chart=text_bellow_chart, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=102, output_format=output_format, show_plots=show_plots)

    # Without legend etc as png
    name_to_save=f'002_1_IN_no_of_IN_nucleotides_summary_violinplot_small'
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value', number_of_columns=len(pivot_table_avg['variable'].unique()), show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, 
    order=None, project_output_dir=charts_output_dir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=None, xtickt_rotation=0, ylabel=ylabel, xlabel='No. of inserted nucleotides', yaxis_formater=100, name_to_save=name_to_save, text_bellow_chart=None, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=102, output_format='png', show_plots=show_plots)

    # Without legend etc as svg
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='variable', y='value', number_of_columns=len(pivot_table_avg['variable'].unique()), show_datapoints=True, jitter=0.1, stripplot_offset=12, dodge=False, clip_right=True, 
    order=None, project_output_dir=charts_output_dir, hue=None, integrated_properties_dic=False, color_dic=color_dic, show_hue_medians=False, title=None, xtickt_rotation=0, ylabel=ylabel, xlabel='No. of inserted nucleotides', yaxis_formater=100, name_to_save=name_to_save, text_bellow_chart=None, y_legend_coordinate=-0.05, data_source=data_source, y_min=0, y_max=102, output_format='svg', show_plots=show_plots)

    # Save pivot table as tsv
    pivot_table_avg.to_csv(Path(str(charts_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
    stat_info = pivot_table_avg.groupby('variable').value.describe().unstack()
    save_stat(stat_info, name_to_save, charts_output_dir)
    ####################################################################################,
    # FILTER ORIGINAL DATAFRAME FOR SAMPLES WHERE 1NT INSERTION IS >= 80%
    #####################################################################################
    # Filter pivot_table_avg for variable == '1' and value >= 80
    pivot_table_avg_1 = pivot_table_avg[(pivot_table_avg['variable'] == '1') & (pivot_table_avg['value'] >= 80)]
    # Make list of unique combinations of ['kingdom', 'source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation'] separated by '_' from pivot_table_avg_1 and save it in a list called filter_list
    # print(min(pivot_table_avg_1['value']))
    filter_list = pivot_table_avg_1['ID_no_rep'].tolist()
    print(len(filter_list))
    print(f'{len(filter_list)} -> Number of unique combinations of filter in pivot_table_avg_1, which will be used to filter original dataframe for samples with 1 inserted nucleotide and value >= 80%.')
    # Save filter_list as tsv
    pd.DataFrame(filter_list, columns=['ID_no_rep']).to_csv(Path(str(table_output_dir) + '/' + f'04_1_filter_list_1ntIN_above80_{date}.tsv'), sep='\t', index=False)
    print(df.columns)
    # Filter df for rows where filter column is in filter_list
    df_filtered = df[df['ID_no_rep'].isin(filter_list)]
    print(f'{df_filtered.shape} -> After filtering original dataframe for samples with 1 inserted nucleotide and value >= 80% based on filter_list.')
    # Save df_filtered as tsv -> 04_filtered_efficacy_and_1ntIN80_samples_{date}.tsv in table_output_dir
    df_filtered.to_csv(Path(str(table_output_dir) + '/' + f'04_2_filtered_1ntIN_above80_samples_{date}.tsv'), sep='\t', index=False)
    # If shape of df_filtered is 0, skip the rest of the code and continue to the next run
    if df_filtered.shape[0] == 0:
        print(f'No samples with 1 inserted nucleotide and value >= 80% found in run {run+1}/{no_of_runs}, skipping the rest of the code and continuing to the next run.')
        continue
    list_4th_nt = default_list + ['4th_nt']
    # Filter df_IN for 1nd insertion only
    df_1IN = df_filtered[df_filtered['n_inserted'] == 1]
    # Get data
    pivot = pd.pivot_table(df_1IN, index=list_4th_nt, columns='4th_nt_inserted', values='%Reads', aggfunc='sum', fill_value=0)
    # reset pivot and remove index
    df_pivot = pivot.reset_index()
    # remove index name
    df_pivot.columns.name = None
    # Check if there are columns True and False, if not add them with 0 values
    if True not in df_pivot.columns:
        df_pivot[True] = 0
    if False not in df_pivot.columns:
        df_pivot[False] = 0
    # Normalise rows to 100%
    df_pivot.iloc[:, len(list_4th_nt):] = df_pivot.iloc[:, len(list_4th_nt):].div(df_pivot.iloc[:, len(list_4th_nt):].sum(axis=1), axis=0) * 100
    # Melt dataframe
    df_pivot = pd.melt(df_pivot, id_vars=list_4th_nt, value_vars=[True, False])
    # remove from list_4th_nt the 'replicate' element
    list_4th_nt.remove('replicate')
    # Make pivot table with all columns without column replicate and make average of values
    pivot_table_avg = create_pivot_table(df_pivot, list_4th_nt, columns='variable', values='value', aggfunc='mean')
    # Met pivot table_avg reset index
    pivot_table_avg = pd.melt(pivot_table_avg.reset_index(), id_vars=list_4th_nt, value_vars=[True, False])
    # Filter out False values of variable
    pivot_table_avg = pivot_table_avg[pivot_table_avg['variable'] == True]
    # Rename column variable to '4th_nt_inserted'
    pivot_table_avg = pivot_table_avg.rename(columns={'variable': '4th_nt_inserted'})
    len(pivot_table_avg['ID_no_rep'].unique())

    # make hue_order according to present 4th nucleotides
    nucleotide_order_list = {'A':1, 'C':2, 'G':3, 'T':4}
    nucleotide_list = pivot_table_avg['4th_nt'].unique()
    sort_list = []
    for nucleotide in nucleotide_list:
        sort_list.append(nucleotide_order_list[nucleotide])
    hue_order = [x for _, x in sorted(zip(sort_list, nucleotide_list))]     

    title='Is the 4th nucleotide the inserted nucleotide? Mean of replicates'
    no_of_samples = len(pivot_table_avg)
    ylabel='Normalised percentage\n of inserted 4th nucleotides'
    text_bellow_chart = f'Only 1nt IN around the cut are plotted.\
                        \nSequences without DEL/SUB around the cut.\
                        \nIn individual sample True + False = 100%.\
                        \nHorizontal lines show medians for individual 4th nucleotide.\
                        \nOnly data where 1nt IN is >= 80% of all IN events in the sample are plotted.\
                        \nNo of samples: {no_of_samples}'
    name_to_save = f'002_2_IN_4th_nt_inserted_normalised_violinplot_nts_sorted'
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
        , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
        , show_hue_medians=False\
        , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
        , clip_right=True\
        , title=title, xtickt_rotation=0, ylabel=ylabel, xlabel='4th nucleotide inserted'\
        , order=None, plot_legend=True \
        , integrated_properties_dic=True, text_bellow_chart=text_bellow_chart\
        , y_legend_coordinate=-0.05, yaxis_formater=100\
        , x_max=1\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir \
        , data_source=data_source, output_format=output_format\
        , dpi=200\
        , y_min=0\
        , y_max=102, show_plots=show_plots
            )

    # Without legend etc as png
    name_to_save = f'002_2_IN_4th_nt_inserted_normalised_violinplot_nts_sorted_small'
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
        , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
        , show_hue_medians=False\
        , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
        , clip_right=True\
        , title=None, xtickt_rotation=0, ylabel=ylabel, xlabel='4th nucleotide inserted'\
        , order=None, plot_legend=True \
        , integrated_properties_dic=True, text_bellow_chart=None\
        , y_legend_coordinate=-0.05, yaxis_formater=100\
        , x_max=1\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir \
        , data_source=data_source, output_format='png'\
        , dpi=200\
        , y_min=0\
        , y_max=102, show_plots=show_plots
            )

    # Without legend etc as svg
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
        , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
        , show_hue_medians=False\
        , hue='4th_nt', hue_order=hue_order, hue_only_for_data_points=True \
        , clip_right=True\
        , title=None, xtickt_rotation=0, ylabel=ylabel, xlabel='4th nucleotide inserted'\
        , order=None, plot_legend=True \
        , integrated_properties_dic=True, text_bellow_chart=None\
        , y_legend_coordinate=-0.05, yaxis_formater=100\
        , x_max=1\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir \
        , data_source=data_source, output_format='svg'\
        , dpi=200\
        , y_min=0\
        , y_max=102, show_plots=show_plots
            )
    # Save pivot table as tsv
    pivot_table_avg.to_csv(Path(str(charts_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
    # Save statistics
    stat_info = pivot_table_avg.groupby('4th_nt_inserted').describe().unstack()
    save_stat(stat_info, name_to_save, charts_output_dir)

    # Full plot as png
    name_to_save = f'002_3_IN_4th_nt_inserted_normalised_violinplot'
    # palette = [True:'#AAAAAA']
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
        , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
        , show_hue_medians=False\
        , hue=None, hue_order=None, hue_only_for_data_points=False \
        , color_dic={'True':"#D3D3D3"}\
        , clip_right=True\
        , title=title, xtickt_rotation=0, ylabel=ylabel, xlabel='4th nucleotide inserted'\
        , order=None, plot_legend=False \
        , integrated_properties_dic=False, text_bellow_chart=text_bellow_chart\
        , y_legend_coordinate=-0.05, yaxis_formater=100\
        , x_max=1\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir \
        , data_source=data_source, output_format=output_format\
        , dpi=200\
        , y_min=0\
        , y_max=102, show_plots=show_plots
            )

    # Without legend etc as png
    name_to_save = f'002_3_IN_4th_nt_inserted_normalised_violinplot_small'
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
        , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
        , show_hue_medians=False\
        , hue=None, hue_order=None, hue_only_for_data_points=False \
        , color_dic={'True':"#D3D3D3"}\
        , clip_right=True\
        , title=None, xtickt_rotation=0, ylabel=ylabel, xlabel='4th nucleotide inserted'\
        , order=None, plot_legend=False \
        , integrated_properties_dic=False, text_bellow_chart=None\
        , y_legend_coordinate=-0.05, yaxis_formater=100\
        , x_max=1\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir \
        , data_source=data_source, output_format='png'\
        , dpi=200\
        , y_min=0\
        , y_max=102, show_plots=show_plots
            )

    # Without legend etc as svg
    plot = plot_sns_violin_and_print(data=pivot_table_avg, x='4th_nt_inserted', y='value', number_of_columns=1.2\
        , show_datapoints=True, jitter=0.3, stripplot_offset=25, dodge=True\
        , show_hue_medians=False\
        , hue=None, hue_order=None, hue_only_for_data_points=False \
        , color_dic={'True':"#D3D3D3"}\
        , clip_right=True\
        , title=None, xtickt_rotation=0, ylabel=ylabel, xlabel='4th nucleotide inserted'\
        , order=None, plot_legend=False \
        , integrated_properties_dic=False, text_bellow_chart=None\
        , y_legend_coordinate=-0.05, yaxis_formater=100\
        , x_max=1\
        , name_to_save=name_to_save\
        , project_output_dir=charts_output_dir \
        , data_source=data_source, output_format='svg'\
        , dpi=200\
        , y_min=0\
        , y_max=102, show_plots=show_plots
            )
    # Save pivot table as tsv
    pivot_table_avg.to_csv(Path(str(charts_output_dir) + '/' + name_to_save + '_' + date + '.tsv'), sep='\t')
    # Save statistics
    stat_info = pivot_table_avg.groupby('4th_nt_inserted').describe().unstack()
    save_stat(stat_info, name_to_save, charts_output_dir)

print(f'Run {run+1}/{no_of_runs} sample {filter_value} analysed, charts and tables saved in {charts_output_dir} and {table_output_dir} respectively.')
# %%
