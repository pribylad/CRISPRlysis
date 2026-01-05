# This script splits a merged CRISPR analysis results table into multiple tables based on a specified parameter (e.g., 'cell_cycle').
# \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\03_merged_tables\03_filtered_efficacy_samples_limit10efficacy_260105.tsv
# Each resulting table is saved as a separate TSV file in the output directory.

# %%
import os
import pandas as pd
from datetime import datetime
##########################################################################################
input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/03_filtered_efficacy_samples_limit10efficacy_260105.tsv'
output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/'
filter_criterion = 'cell_cycle'  # parameter to split the table by
##########################################################################################
split_parameter_dic = {'cell_cycle': ['dividing', 'nondividing']}

if filter_criterion not in split_parameter_dic:
    raise ValueError(f"Filter criterion {filter_criterion} not recognized")
else:
    split_parameter = split_parameter_dic[filter_criterion]
# %%
##########################################################################################
# load the filtered dataframe
df = pd.read_csv(input_file, sep='\t')
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