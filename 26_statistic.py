# %%
# Load table from where the statistic should be done
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels
from statsmodels.stats.multitest import multipletests
print("funguje:", multipletests.__name__)
# %%
# Decide if you want to analyse only one table or compare two tables (e.g. plant vs. animal)
One_table = False

# # 1a - IN/DEL/SUB/... categories
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/001_1_INDELSUB_summary_violinplot_260518.tsv"
# column_to_analyze = "Event_around_cut"
# add_no_mutation = True
# make_index_from = ["ID_no_rep"]

# # 1b - number of IN nucleotides
# # df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# column_to_analyze = "variable"
# add_no_mutation = False
# make_index_from = ["ID_no_rep"]

# # 1c - 9 del categories
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# column_to_analyze = "variable"
# add_no_mutation = False
# make_index_from = ["organism_cell_line_hpi", "gene", "gRNA_seq"]

# 1d - DEL end participating
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# column_to_analyze = "variable"
# add_no_mutation = False
# make_index_from = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]

# TODO this can not be done as the two tables have different number of samples and they are independent, different test should be used
# #2
# # Paths to the two tables to compare, e.g. plant vs. animal
# df_path_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
# df_path_2 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/001_1_INDELSUB_summary_violinplot_260518.tsv"
# # Make output path for statistic results, same as input but with suffix "_statistic"
# output_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/00_statistic/001_1_INDELSUB_summary_violinplot_260518_statistic_plantVSanimal.tsv"
# column_to_analyze = "Event_around_cut"
# add_no_mutation = True
# make_index_from = ["ID_no_rep"]


if One_table == False:
    df1 = pd.read_csv(df_path_1, sep="\t")
    df2 = pd.read_csv(df_path_2, sep="\t")

    # Make index from columns in make_index_from
    if len(make_index_from) > 1:
        df1['index'] = df1[make_index_from].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
        df2['index'] = df2[make_index_from].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    else:
        df1['index'] = df1[make_index_from[0]]
        df2['index'] = df2[make_index_from[0]]
    # Check unique values in index
    print("Unique values in index for DataFrame 1:", df1['index'].nunique())
    print("Unique values in index for DataFrame 2:", df2['index'].nunique())

    pivot_table1 = df1.pivot_table(index='index', columns=column_to_analyze, values='value')
    pivot_table2 = df2.pivot_table(index='index', columns=column_to_analyze, values='value')

    if add_no_mutation:
        pivot_table1['no_mutation'] = 100 - pivot_table1.sum(axis=1)
        pivot_table2['no_mutation'] = 100 - pivot_table2.sum(axis=1)

    df_stat1 = pivot_table1.reset_index()
    df_stat2 = pivot_table2.reset_index()

    #  rename columns of df_stat1 and add suffix "_plant" to each column except ID_no_rep
    df_stat1 = df_stat1.rename(columns={col: col + "_plant" for col in df_stat1.columns if col != 'index'})
    # rename columns of df_stat2 and add suffix "_animal" to each column except ID_no_rep
    df_stat2 = df_stat2.rename(columns={col: col + "_animal" for col in df_stat2.columns if col != 'index'})

    # merge df_stat1 and df_stat2 by ID_no_rep and fill Nan with 0
    df_stat = pd.merge(df_stat1, df_stat2, on='index', how='outer')
    df_stat = df_stat.fillna(0)
else:
    output_path = df_path.replace(".tsv", "_statistic.tsv")
    df = pd.read_csv(df_path, sep="\t")
    # Make index from columns in make_index_from
    if len(make_index_from) > 1:
        df['index'] = df[make_index_from].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    else:
        df['index'] = df[make_index_from[0]]
    # Check unique values in index
    print("Unique values in index:", df['index'].nunique())
    # Make pivot table for the statistic, rows ID_no_rep, columns Event_around_cut, values value
    pivot_table = df.pivot_table(index='index', columns=column_to_analyze, values='value')
    # add column "no_mutation" with value 100-sum of the other columns
    if add_no_mutation:
        pivot_table['no_mutation'] = 100 - pivot_table.sum(axis=1)

    # make dataframe from pivot_table and use index as index
    df_stat = pivot_table.reset_index()
df_stat.index = df_stat['index']
# remove column ID_no_rep
df_stat = df_stat.drop(columns=['index'])
# Print names of columns
print(f"{len(df_stat.columns.tolist())} columns in the data for statistic:")
print(df_stat.columns.tolist())

print("Data prepared for statistic:")
print(df_stat.head())

# get list of columns
columns = df_stat.columns.tolist()

# CLR transformace
def clr(x):
    x = np.array(x)
    log_x = np.log(x + 0.001) # 0.001 is a pseudocount to avoid log(0), smaller is not  recommended to reduce distortion near zero
    return log_x - np.mean(log_x, axis=1, keepdims=True)

X = df_stat[columns].values
X_clr = clr(X)

print("Data after CLR transformation:")
print(X_clr[:5])
# ANOVA‑like test (all categories vs. zero in CLR)
p_list = []
for i, cat in enumerate(X_clr.T):
    t, p = stats.ttest_1samp(cat, 0)
    p_list.append(p)

# Corection FDR
_, p_corr, _, _ = multipletests(p_list, method="fdr_bh")

print("CLR‑ANOVA‑like test (all categories vs. zero):")
for k, cat in enumerate(columns):
    print(f"{cat:12} | p = {p_list[k]:.4f} | p_fdr = {p_corr[k]:.4f}")

# Pair comparison of all pairs
categories = columns
tests = []
pvals = []

for i in range(len(categories)):
    for j in range(i + 1, len(categories)):
        catA = categories[i]
        catB = categories[j]
        A = df_stat[catA].values
        B = df_stat[catB].values
        diff = A - B
        _, pval = stats.ttest_1samp(diff, 0)   # nebo stats.wilcoxon(diff) pro Wilcoxon
        tests.append(f"{catA} vs. {catB}")
        pvals.append(pval)

# Korekce na mnohonásobná porovnání
_, pvals_corrected, _, _ = multipletests(pvals, method="fdr_bh")

print("\All pair-wise comparisons (FDR corrected):")
for t, p, pc in zip(tests, pvals, pvals_corrected):
    print(f"{t:20} | p = {p:.4f} | p_fdr = {pc:.4f}")

# print pairwise comparison results in a table
results_df = pd.DataFrame({
    "Comparison": tests,
    "p-value": pvals,
    "p-value (FDR)": pvals_corrected
})
print("\nPair-wise comparison results:")
print(results_df)

# save results to tsv file
results_df.to_csv(output_path, sep="\t", index=False)
# %%
