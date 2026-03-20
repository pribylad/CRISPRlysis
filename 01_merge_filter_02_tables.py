# Merge all tables in folder \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\02_tables\ and load them into a single dataframe
# %%
import os
import glob
import pandas as pd
from datetime import datetime
##########################################################################################
source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/All_260320'
output_file_name = '01_merged_02tables_NoFilter_All_' + datetime.now().strftime('%y%m%d') + '.tsv'

# source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/All_260305'
# output_file_name = '01_merged_02tables_NoFilter_All_' + datetime.now().strftime('%y%m%d') + '.tsv'
# source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/part_A/Animals'
# output_file_name = '01_merged_02tables_NoFilter_Animals_' + datetime.now().strftime('%y%m%d') + '.tsv'
# source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/part_A/Plants'
# output_file_name = '01_merged_02tables_NoFilter_Plants_' + datetime.now().strftime('%y%m%d') + '.tsv'

# source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/Kamoen_2024/WT'
# output_file_name = '01_merged_02tables_NoFilter_Kamoen_WT_' + datetime.now().strftime('%y%m%d') + '.tsv'

# source_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/Weiss_0_8_7_4'
# output_file_name = '01_merged_02tables_NoFilter_Weiss_0874_' + datetime.now().strftime('%y%m%d') + '.tsv'
##########################################################################################

output_directory = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Make dictionary to assign plant/animal/yeast origin and colour for scatter plot
organism_dict = {
    'A. thaliana Col-0': 'plant',
    'A. thaliana': 'plant',
    'Arabidopsis thaliana': 'plant',
    'HCT116': 'animal',
    'HEK293': 'animal',
    'Hela cell lines': 'animal',
    'K562': 'animal',
    'Nicotiana benthamiana': 'plant',
    'Petunia hybrida': 'plant',
    'Silene latifoliaa': 'plant',
    'Solanum lycopersicum': 'plant',
    'Homo sapiens': 'animal',
    'Xenopus tropicalis ': 'animal'
}
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

# Check if all unique values in the 'organism' column of the merged dataframe are present in the organism_dict
unique_organisms = merged_df['organism'].unique()
for organism in unique_organisms:
    if organism not in organism_dict:
        # Rise an error if an organism is not present in the organism_dict
        raise ValueError(f"Organism '{organism}' is not present in the organism_dict. Please add it to the dictionary.")
# Add column 'source' to merged_df to indicate the source of the data (plant/animal/yeast)
merged_df['kingdom'] = merged_df['organism'].map(organism_dict)
# Repair cell line column
merged_df['cell_line'] = merged_df['cell_line'].fillna('xxx')
# Change all hours_post_inoculation values to string
merged_df['hours_post_inoculation'] = merged_df['hours_post_inoculation'].astype(str)
# Change A. thaliana Col-0 to Arabidopsis thaliana
merged_df['organism'] = merged_df['organism'].replace('A. thaliana Col-0', 'Arabidopsis thaliana')
# %%
# Make a list of unique 'name' values and asign a number to each unique name in a new column 'sample_ID'
list_for_sample_ID = ['source', 'organism', 'cell_line', 'mutant', 'gene', 'gRNA_seq', 'hours_post_inoculation']
merged_df['ID_no_rep'] = merged_df[list_for_sample_ID].apply(lambda x: '_'.join(x.astype(str)), axis=1)
print(merged_df['ID_no_rep'].unique())

#  Remove column Unnamed: 0 if it exists
if 'Unnamed: 0' in merged_df.columns:
    merged_df = merged_df.drop(columns=['Unnamed: 0']) 

# Save merged dataframe to TSV
output_file = os.path.join(output_directory, output_file_name)
merged_df.to_csv(output_file, sep='\t', index=False)
##########################################################################################
# %%
# print size of the merged dataframe
print(f"Merged dataframe shape: {merged_df.shape}")
print(merged_df['kingdom'].unique())
print(f"Merged dataframe organisms: {merged_df['organism'].unique()}")
# %%
print(merged_df.columns)