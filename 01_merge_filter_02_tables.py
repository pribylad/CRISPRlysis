# Merge all tables in folder \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\02_tables\ and load them into a single dataframe
# %%
import os
import glob
import pandas as pd
from datetime import datetime
##########################################################################################
source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/part_A'
output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/'
output_file_name = '01_merged_02tables_NoFilter_' + datetime.now().strftime('%y%m%d') + '.tsv'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
##########################################################################################
def load_tsv_files_from_directory(directory_path):
    tsv_files = glob.glob(os.path.join(directory_path, "*.tsv"))
    dataframes = []
    for file in tsv_files:
        df = pd.read_csv(file, sep='\t')
        dataframes.append(df)
    return dataframes
def merge_dataframes(dataframes):
    merged_df = pd.concat(dataframes, ignore_index=True)
    return merged_df
##########################################################################################
dataframes = load_tsv_files_from_directory(source_directory)
merged_df = merge_dataframes(dataframes)
# Save merged dataframe to TSV
output_file = os.path.join(output_directory, output_file_name)
merged_df.to_csv(output_file, sep='\t', index=False)
##########################################################################################
# %%
# print size of the merged dataframe
print(f"Merged dataframe shape: {merged_df.shape}")
print(merged_df['source'].unique())
print(f"Merged dataframe organisms: {merged_df['organism'].unique()}")
# %%
