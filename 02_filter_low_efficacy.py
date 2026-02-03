# Description:
# This script processes a merge tables of 02_tables CRISPR analysis results, identifies and filter out samples with low editing efficacy.
# Merged files are: 02_pivot_table_no_filter_<date>.tsv
# Filtered files are: 03_filtered_efficacy_samples_limit<efficacy>%efficacy_<date>.tsv
# Low efficacy samples list are: 03_low_efficacy_samples_limit<efficacy>%efficacy_<date>.txt

# %%
import os
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
##########################################################################################
x=95  # efficacy threshold percentage
output_directory = f'/home/lefnui/CRISPRlysis_outputs/03_merged_tables/filtered_{x}/'

# input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_merged_02tables_NoFilter_Plants_260126.tsv'
# output_nofilter_tsv = f'02_pivot_table_no_filter_{datetime.now().strftime("%y%m%d")}.tsv'
# output_tsv = f'03_filtered_efficacy_samples_limit{100-x}efficacy_Plants_{datetime.now().strftime("%y%m%d")}.tsv'
# output_low_efficacy_txt = f'03_low_efficacy_samples_limit{100-x}efficacy_Plants_{datetime.now().strftime("%y%m%d")}.txt'

# cell_line = 'HCT116'
# hpt = '48'
# input_file = f'/home/lefnui/CRISPRlysis_outputs/van_Overbeek_2016/02_mastertable_{cell_line}_RNP_{hpt}_260108.tsv'
# output_nofilter_tsv = f'02_pivot_table_no_filter_{cell_line}_RNP_{hpt}_{datetime.now().strftime("%y%m%d")}.tsv'
# output_tsv = f'03_filtered_efficacy_samples_limit{100-x}efficacy_{cell_line}_RNP_{hpt}_{datetime.now().strftime("%y%m%d")}.tsv'
# output_low_efficacy_txt = f'03_low_efficacy_samples_limit{100-x}efficacy_{cell_line}_RNP_{hpt}_{datetime.now().strftime("%y%m%d")}.txt'

# input_file = '/home/lefnui/CRISPRlysis_outputs/03_merged_tables/01_merged_02tables_NoFilter_Kamoen_WT_260130.tsv'
# output_nofilter_tsv = f'02_pivot_table_no_filter_Kamoen_WT_{datetime.now().strftime("%y%m%d")}.tsv'
# output_tsv = f'03_filtered_efficacy_samples_limit{100-x}efficacy_Kamoen_WT_{datetime.now().strftime("%y%m%d")}.tsv'
# output_low_efficacy_txt = f'03_low_efficacy_samples_limit{100-x}efficacy_Kamoen_WT_{datetime.now().strftime("%y%m%d")}.txt'

filter_value = '37'
input_file = f'/home/lefnui/CRISPRlysis_outputs/Illa_Berenguer_et_al_2023_{filter_value}/02_mastertable_PDS_WT_{filter_value}_260203.tsv'
output_nofilter_tsv = f'02_pivot_table_no_filter_PDS_WT_{filter_value}_{datetime.now().strftime("%y%m%d")}.tsv'
output_tsv = f'03_filtered_efficacy_samples_limit{100-x}efficacy_PDS_WT_{filter_value}_{datetime.now().strftime("%y%m%d")}.tsv'
output_low_efficacy_txt = f'03_low_efficacy_samples_limit{100-x}efficacy_PDS_WT_{filter_value}_{datetime.now().strftime("%y%m%d")}.txt'
##########################################################################################
# If output_directory does not exists, create on
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
# load the merged dataframe
df = pd.read_csv(input_file, sep='\t')
##########################################################################################
# %%
# create pivot table
pivot_table = pd.pivot_table(df, index=['sample', 'gRNA_seq', 'replicate'], columns=['Event'], values=['%Reads'], aggfunc='sum', fill_value=0)
print(pivot_table)
# save pivot table to TSV
output_file = os.path.join(output_directory, output_nofilter_tsv)
pivot_table.to_csv(output_file, sep='\t')
##########################################################################################
# if no_mutation column contains any value higher than 100, print a warning and the rows where it happens
if (pivot_table[('%Reads', 'no_mutation')] > 100).any():
    print("Warning: %Reads for no_mutation exceeds 100 in the following rows:")
    print(pivot_table[pivot_table[('%Reads', 'no_mutation')] > 100])
# %%
# plot column 'no_mutation' as histogram
plt.figure(figsize=(10, 6))
plt.hist(pivot_table[('%Reads', 'no_mutation')], bins=50, color='blue', alpha=0.7)
plt.title('Histogram of No Mutation %Reads')
plt.xlabel('%Reads of No Mutation')
plt.ylabel('Frequency')
plt.grid(axis='y', alpha=0.75)
plt.show()
# %%
# Make list of samples where sum of %Reads for no_mutation is greater than x
low_efficacy_samples = pivot_table[pivot_table[('%Reads', 'no_mutation')] >= x].index.tolist()
print( f"Samples with low efficacy (no_mutation %Reads >= {x}):")
for sample in low_efficacy_samples:
    print(sample)
# Save low efficacy samples to a text file
output_file = os.path.join(output_directory, output_low_efficacy_txt)
with open(output_file, 'w') as f:
    for sample in low_efficacy_samples:
        f.write(f"{sample}\n")
# %%
# remove low efficacy samples from the original dataframe
filtered_df = df[~df.set_index(['sample', 'gRNA_seq', 'replicate']).index.isin(low_efficacy_samples)]
# save filtered dataframe to TSV
output_file = os.path.join(output_directory, output_tsv)
filtered_df.to_csv(output_file, sep='\t', index=False)
# %%
pivot_table = pd.pivot_table(filtered_df, index=['sample', 'gRNA_seq', 'replicate'], columns=['Event'], values=['%Reads'], aggfunc='sum', fill_value=0)
plt.hist(pivot_table[('%Reads', 'no_mutation')], bins=50, color='blue', alpha=0.7)