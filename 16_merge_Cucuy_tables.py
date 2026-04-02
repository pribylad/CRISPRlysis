# pylint: disable=line-too-long
# pylint: disable=invalid-name

# This code merges the Suppplementary data from the Cucuy et apl., 2026 (data from "within_atac_peak", "leaf_exp", "t0_exp", and "description" columns) with individual tables
# \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\02_source_tables\All_260330\02_mastertable_WT_Solanum lycopersicum_Exo_260328.tsv
# \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\02_source_tables\All_260330\02_mastertable_WT_Solanum lycopersicum_Int_260328.tsv
# \\wsl$\Ubuntu\home\lefnui\CRISPRlysis_outputs\02_source_tables\All_260330\02_mastertable_WT_Solanum lycopersicum_Pro_260328.tsv

# %%
import os
import pandas as pd
from datetime import datetime
##########################################################################################
input_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/All_260330/'
output_directory = '/home/lefnui/CRISPRlysis_outputs/02_source_tables/All_260330/'

# Load the Cucuy supplementary data
cucuy_data_file = '/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/03_maybe_to use/Cucuy_et_al_2026/tpj70809-sup-0002-datas1.xlsx'
cucuy_df = pd.read_excel(cucuy_data_file, sheet_name='main_supp_table')
cucuy_df
# %%
print(f'Number of unique guide sequences in Cucuy data: {cucuy_df["Guide sequence"].nunique()}')
print(f'Number of duplicated guide sequences in Cucuy data: {cucuy_df["Guide sequence"].duplicated().sum()}')
# %%
# Load the individual tables and merge them with the Cucuy data
# df_exo = pd.read_csv(os.path.join('/home/lefnui/CRISPRlysis_outputs/Cucuy_et_al_2026_Exo', '02_mastertable_WT_Solanum lycopersicum_Exo_260328.tsv'), sep='\t')
# df_int = pd.read_csv(os.path.join('/home/lefnui/CRISPRlysis_outputs/Cucuy_et_al_2026_Int', '02_mastertable_WT_Solanum lycopersicum_Int_260328.tsv'), sep='\t')
# df_pro = pd.read_csv(os.path.join('/home/lefnui/CRISPRlysis_outputs/Cucuy_et_al_2026_Pro', '02_mastertable_WT_Solanum lycopersicum_Pro_260328.tsv'), sep='\t')
df_05_pivot = pd.read_csv(os.path.join('/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/02_filtered_efficacy_95/04_Solanum_lycopersicum_M82_48_260402/05_pivot_DEL_orientation_normalised_avg_260402.tsv'), sep='\t')
# To each row in df_exo, df_int and df_pro, add the columns "within_atac_peak", "leaf_exp", "t0_exp", and "description" from cucuy_df based on matching "Guide sequence" column in cucuy_df with "gRNA_seq" column in df_exo, df_int and df_pro
def merge_with_cucuy(df, cucuy_df):
    merged_df = pd.merge(df, cucuy_df[['Guide sequence', 'within_atac_peak', 'leaf_exp', 't0_exp', 'description']], left_on='gRNA_seq', right_on='Guide sequence', how='left')
    merged_df.drop(columns=['Guide sequence'], inplace=True)
    return merged_df
# df_exo_merged = merge_with_cucuy(df_exo, cucuy_df)
# df_int_merged = merge_with_cucuy(df_int, cucuy_df)
# df_pro_merged = merge_with_cucuy(df_pro, cucuy_df)
df_05_pivot_merged = merge_with_cucuy(df_05_pivot, cucuy_df)
# Save the merged tables as new TSV files
# %%
output_file_exo = os.path.join(output_directory, '02_mastertable_WT_Solanum_lycopersicum_Exo_merged_' + datetime.now().strftime('%y%m%d') + '.tsv')
output_file_int = os.path.join(output_directory, '02_mastertable_WT_Solanum_lycopersicum_Int_merged_' + datetime.now().strftime('%y%m%d') + '.tsv')
output_file_pro = os.path.join(output_directory, '02_mastertable_WT_Solanum_lycopersicum_Pro_merged_' + datetime.now().strftime('%y%m%d') + '.tsv')
output_file_05_pivot = os.path.join('/mnt/d/Lab210/DATA/04_1_bp_MH_DNA_repair/data/Obsidian_1bp_CMEJ_vault/01_data/00_final_analysis/03_merged_tables/02_filtered_efficacy_95/04_Solanum_lycopersicum_M82_48_260402', '05_pivot_DEL_orientation_normalised_avg_merged_' + datetime.now().strftime('%y%m%d') + '.tsv')
# df_exo_merged.to_csv(output_file_exo, sep='\t', index=False)
# df_int_merged.to_csv(output_file_int, sep='\t', index=False)
# df_pro_merged.to_csv(output_file_pro, sep='\t', index=False)
df_05_pivot_merged.to_csv(output_file_05_pivot, sep='\t', index=False)
# print(f'Saved merged Exo table to {output_file_exo} with shape {df_exo_merged.shape}')
# print(f'Saved merged Int table to {output_file_int} with shape {df_int_merged.shape}')
# print(f'Saved merged Pro table to {output_file_pro} with shape {df_pro_merged.shape}')
print(f'Saved merged 05_pivot table to {output_file_05_pivot} with shape {df_05_pivot_merged.shape}')
# %%
df_05_pivot_merged.head()
# %%
