# pylint: disable=line-too-long
# pylint: disable=invalid-name
#  %%
import pandas as pd
import os

# %%
path_to_df = '/home/lefnui/CRISPRlysis_outputs/Ellison_et_al_2020/inf/02_mastertable_PDS_WT_infected_260203.tsv'
df = pd.read_csv(path_to_df, sep='\t')
df['hours_post_inoculation'].unique()
df['hours_post_inoculation'] = '336'

df['hours_post_inoculation'].unique()
# save the modified dataframe and replace old file
df.to_csv(path_to_df, sep='\t', index=False)


# %%
df_noFilter = pd.read_csv('/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/01_no_filter_all_samples/01_merged_02tables_NoFilter_All_260402.tsv', sep='\t')
df_noFilter.head()

# %%
print(f'Number of unique guide sequences in no filter data: {df_noFilter["gRNA_seq"].nunique()}')
print(f'Number of source publications in no filter data: {df_noFilter["source"].nunique()}')
print(f'sources: {df_noFilter["source"].unique()}')

# %%
df_filter = pd.read_csv('/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/02_filtered_efficacy_95/03_filtered_efficacy_samples_limit5efficacy_All_260402.tsv', sep='\t')
# %%
print(f'Number of unique guide sequences in filtered data: {df_filter["gRNA_seq"].nunique()}')
print(f'Number of source publications in filtered data: {df_filter["source"].nunique()}')
print(f'sources: {df_filter["source"].unique()}')

# %%
# For all .tsv files in the directory /home/lefnui/CRISPRlysis_outputs/02_source_tables/All_260330/, read them as pandas dataframe and print how many unique gRNAs they have
directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/All_260330/'
for filename in os.listdir(directory):
    if filename.endswith('.tsv'):
        df = pd.read_csv(os.path.join(directory, filename), sep='\t')
        print(f'{filename}: {df["gRNA_seq"].nunique()} unique gRNAs')