# This script splits a merged CRISPR analysis results table into multiple tables based on a specified parameter (e.g., 'cell_cycle').
# \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\03_merged_tables\03_filtered_efficacy_samples_limit10efficacy_260105.tsv
# Each resulting table is saved as a separate TSV file in the output directory.

# %%
import os
import pandas as pd
from datetime import datetime
##########################################################################################
input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_efficacy_samples_limit5efficacy_260109.tsv'
output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/'

# load the filtered dataframe
df = pd.read_csv(input_file, sep='\t')

# %%
##########################################################################################
##########################################################################################
# FILTER BY PARAMETER WHICH DOES NOT NEED ADDITIONAL PROCESSING
##########################################################################################
##########################################################################################

filter_criterion = 'cell_cycle'  # parameter to split the table by

##########################################################################################
split_parameter_dic = {'cell_cycle': ['dividing', 'nondividing']}

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
# %%
# split the dataframe and save each part to a separate TSV file
for value in split_parameter:
    split_df = df[df[filter_criterion] == value]
    output_file = os.path.join(output_directory, f'04_split_{filter_criterion}_{value}_{datetime.now().strftime("%y%m%d")}.tsv')
    split_df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved {output_file} with shape {split_df.shape}")
    print(f"Sources in this split: {split_df['source'].unique().tolist()}")
    print("--------------------------------------------------")

# # %%
# ##########################################################################################
# ##########################################################################################
# # FILTER BY PARAMETER WHICH DOES NEED ADDITIONAL PROCESSING
# ##########################################################################################
# ##########################################################################################
# filter_criterion = 'Event_around_cut'  # parameter to split the table by
# df[filter_criterion].unique()
# split_parameter_dic = {'Event_around_cut': [df[filter_criterion].unique()]}
# wanted_parameter = 'IN'
# treshold = 24.193858
# ##########################################################################################
# # Make pivot table to get sum of %Reads for wanted paremeter from the column 'Event_around_cut' where index is 'sample'
# pivot_table = df.pivot_table(index='sample', columns=filter_criterion, values='%Reads', aggfunc='sum', fill_value=0)
# # Get samples where wanted parameter is above treshold
# samples_above_treshold = pivot_table[pivot_table[wanted_parameter] > treshold].index.tolist()
# samples_below_treshold = pivot_table[pivot_table[wanted_parameter] <= treshold].index.tolist()
# # Split the dataframe based on the samples
# df_above_treshold = df[df['sample'].isin(samples_above_treshold)]
# df_below_treshold = df[df['sample'].isin(samples_below_treshold)]
# # Save the split dataframes to separate TSV files
# output_file_above = os.path.join(output_directory, f'04_split_{filter_criterion}_{wanted_parameter}_above_treshold_{datetime.now().strftime("%y%m%d")}.tsv')
# df_above_treshold.to_csv(output_file_above, sep='\t', index=False)
# print(f"Saved {output_file_above} with shape {df_above_treshold.shape}")
# print(f"No of unique samples above treshold: {df_above_treshold['sample'].nunique()}")
# print("--------------------------------------------------")
# output_file_below = os.path.join(output_directory, f'04_split_{filter_criterion}_{wanted_parameter}_below_treshold_{datetime.now().strftime("%y%m%d")}.tsv')
# df_below_treshold.to_csv(output_file_below, sep='\t', index=False)
# print(f"Saved {output_file_below} with shape {df_below_treshold.shape}")
# print(f"No of unique samples below treshold: {df_below_treshold['sample'].nunique()}")
# print("--------------------------------------------------")
# # %%