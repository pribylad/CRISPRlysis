# This script splits a merged CRISPR analysis results table into multiple tables based on a specified parameter (e.g., 'cell_cycle').
# \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\03_merged_tables\03_filtered_efficacy_samples_limit10efficacy_260105.tsv
# Each resulting table is saved as a separate TSV file in the output directory.

# %%
import os
import pandas as pd
from datetime import datetime
##########################################################################################

# input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
# output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/'

input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/03_filtered_efficacy_samples_limit5efficacy_All.tsv'
output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/02_filtered_efficacy_95/with_Ponce_samples/'

# input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_mutants_260413.tsv'
# output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_no_filter_all_samples/mutants/'

# If output directory does not exist, create it
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
# load the filtered dataframe
df = pd.read_csv(input_file, sep='\t')
df['hours_post_inoculation'] = df['hours_post_inoculation'].astype(str)
df['hours_post_inoculation'].unique()
# %%
df.columns
df['source'].unique()
df['mutant'].unique()
# %%
##########################################################################################
##########################################################################################
# FILTER BY PARAMETER WHICH DOES NOT NEED ADDITIONAL PROCESSING
##########################################################################################
##########################################################################################

filter_criterion = 'kingdom'  # parameter to split the table by
split_parameter_dic = {'kingdom': ['plant', 'animal']}
##########################################################################################

# df['index'] = df['organism'] + ' - ' + df['cell_line'] + ' - ' + df['hours_post_inoculation']
# df['index'].unique()
# filter_criterion = 'index'
# # Make dictionary where key is 'index' and value is list of unique valuses o index
# split_parameter_dic = {'index': df['index'].unique().tolist()}

# # For Cucuy et al., 2026
# filter_criterion = 'ProExoInt'
# split_parameter_dic = {'ProExoInt': ['Pro', 'Exo', 'Int']}
# # Keep only the first 3 characters of the values in the column 'ProExoInt' to match the split_parameter_dic
# df['ProExoInt'] = df['ProExoInt'].str[:3]
# # Remove rows where 'ProExoInt' is not in split_parameter_dic
# df = df[df['ProExoInt'].isin(split_parameter_dic['ProExoInt'])]

# # For mutants
# filter_criterion = 'mutant'
# split_parameter_dic = {'mutant': ['WT', 'ku70', 'ku70teb5', 'teb5']}

if filter_criterion not in split_parameter_dic:
    raise ValueError(f"Filter criterion {filter_criterion} not recognized")
else:
    split_parameter = split_parameter_dic[filter_criterion]
##########################################################################################
# get unique values of the split parameter
unique_values = df[filter_criterion].unique()
# Rise error if different values than expected
for value in unique_values:
    if value not in split_parameter:
        raise ValueError(f"Unexpected value '{value}' found in column '{filter_criterion}'")
##########################################################################################
# split the dataframe and save each part to a separate TSV file
for value in split_parameter:
    split_df = df[df[filter_criterion] == value]
    output_file = os.path.join(output_directory, f'04_split_{filter_criterion}_{value}.tsv')
    split_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved {output_file} with shape {split_df.shape}")
    print(f"Sources in this split: {split_df['source'].unique().tolist()}")
    print("--------------------------------------------------")


# %%
