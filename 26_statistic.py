# %%
# Load table from where the statistic should be done
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels
from statsmodels.stats.multitest import multipletests
print("funguje:", multipletests.__name__)
# %%

# # 1a - IN/DEL/SUB/... categories
# # df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
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

# 1c - 9 del categories
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_4_DEL_orientation_summary_violinplot_small.tsv"
df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
column_to_analyze = "variable"
add_no_mutation = False
make_index_from = ["organism_cell_line_hpi", "gene", "gRNA_seq"]
dict_rename = {'upstream_ambiguous':5, 'bidirectional_ambiguous':1, 'downstream_ambiguous':3,
               'upstream_NO_ambi':6, 'bidirectional_NO_ambi':2, 'downstream_NO_ambi':4,
               'cannot_be_determined_ambiguous':7, 'DEL_exceeds_analysed_window':8,
               'gap_outside_of_predicted_cut':9}

# 1d - DEL end participating
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# df_path = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# column_to_analyze = "variable"
# add_no_mutation = False
# make_index_from = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]



output_path = df_path.replace(".tsv", "_statistic.txt")
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
columns
# %%
# rename columns according to dict_rename
df_stat = df_stat.rename(columns=dict_rename)
print("Columns after renaming:")
print(df_stat.columns.tolist())
columns = df_stat.columns.tolist()
columns
# %%
# rename values in columns according to dict_rename
for col in columns:
    df_stat[col] = df_stat[col].replace(dict_rename)
print("Values in columns after renaming:")
for col in columns:
    print(f"{col}: {df_stat[col].unique()}")

# %%
# CLR transformace
def clr(x):
    x = np.array(x)
    log_x = np.log(x + 0.001) # 0.001 is a pseudocount to avoid log(0), smaller is not  recommended to reduce distortion near zero
    return log_x - np.mean(log_x, axis=1, keepdims=True)

X = df_stat[columns].values
X_clr = clr(X)
X_clr
# %%
print("Data after CLR transformation:")
print(X_clr[:5])
# Test normality of each category after CLR transformation with Shapiro-Wilk test
print("Shapiro-Wilk test for normality of each category after CLR transformation:")
for i, cat in enumerate(columns):
    stat, p = stats.shapiro(X_clr[:, i])
    print(f"{cat:12} | Shapiro-Wilk test statistic: {stat:.4f} | p-value: {p:.4f}")
    # Print if the category is normally distributed or not based on p-value
    if p > 0.05:
        print(f"{cat:12} is normally distributed")
    else:
        print(f"{cat:12} is NOT normally distributed")
# %%
# If data is not normally distributed, use non-parametric test (Kruskal-Wallis) to compare groups
print("Kruskal-Wallis test for comparing groups (non-parametric):")
stat, p = stats.kruskal(*[X_clr[:, i] for i in range(X_clr.shape[1])])
print(f"Kruskal-Wallis test statistic: {stat:.4f} | p-value: {p:.4f}")
if p < 0.05:
    print("There is a significant difference between groups (p < 0.05)")
else:
    print("There is no significant difference between groups (p >= 0.05)")
# %%
# Perform Kruskal-Wallis test for pairwise comparisons between categories
print("Pairwise comparisons between categories (Kruskal-Wallis test):")
for i in range(len(columns)):
    for j in range(i + 1, len(columns)):
        catA = columns[i]
        catB = columns[j]
        stat, p = stats.kruskal(X_clr[:, i], X_clr[:, j])
        print(f"{catA:12} vs. {catB:12} | Kruskal-Wallis test statistic: {stat:.4f} | p-value: {p:.4f}")
# %%
# make Dunn posthoc test for pairwise comparisons between categories with Holm-Bonferroni correction
print("Dunn posthoc test for pairwise comparisons between categories (with Holm-Bonferroni correction):")
p_values = []
comparisons = []
for i in range(len(columns)):
    for j in range(i + 1, len(columns)):
        catA = columns[i]
        catB = columns[j]
        stat, p = stats.kruskal(X_clr[:, i], X_clr[:, j])
        p_values.append(p)
        comparisons.append(f"{catA} vs. {catB}")
# Correct p-values for multiple comparisons with Holm-Bonferroni method
_, p_values_corrected, _, _ = multipletests(p_values, method="holm")
i=0
for comp, p, p_corr in zip(comparisons, p_values, p_values_corrected):
    print(f"{comp:20} | p-value: {p:.4f} | p-value (Holm-Bonferroni corrected): {p_corr:.4f}")

    # save results to tsv file
    # if txt file does not exist, create
    if i == 0:
        with open(output_path, "w") as f:
            f.write(f"{comp:20} \t p-value: {p:.4f} \t p-value (Holm-Bonferroni corrected): {p_corr:.4f}")
            f.write("\n\n")
    else:
        with open(output_path, "a") as f:
            f.write(f"{comp:20} \t p-value: {p:.4f} \t p-value (Holm-Bonferroni corrected): {p_corr:.4f}")
            f.write("\n\n")
    i += 1
# convert the output file to a dataframe and print it
df_output = pd.read_csv(output_path, sep="\t", header=None, names=["Comparison", "p-value", "p-value (Holm-Bonferroni corrected)"])
# remove "p-value: " and "p-value (Holm-Bonferroni corrected): " from the columns "p-value" and "p-value (Holm-Bonferroni corrected)"
df_output["p-value"] = df_output["p-value"].str.replace("p-value: ", "").astype(float)
df_output["p-value (Holm-Bonferroni corrected)"] = df_output["p-value (Holm-Bonferroni corrected)"].str.replace("p-value (Holm-Bonferroni corrected): ", "").astype(float)
# replace the output file with the new dataframe
df_output.to_csv(output_path.replace(".txt", ".tsv"), sep="\t", index=False)
print("\nPairwise comparison results:")
print(df_output)
# %%
# print only the comparisons that are NOT significant after Holm-Bonferroni correction
print("\nNon-significant comparisons after Holm-Bonferroni correction (p >= 0.05):")
print(df_output[df_output["p-value (Holm-Bonferroni corrected)"] >= 0.05])
# %%
# For Non-significant comparisons after Holm-Bonferroni correction, plot barplot with error bars for the two categories being compared
import matplotlib.pyplot as plt
import seaborn as sns
for comp, p_corr in zip(df_output["Comparison"], df_output["p-value (Holm-Bonferroni corrected)"]):
    if p_corr >= 0.05:
        catA, catB = comp.split(" vs. ")
        print(catA, catB)
        # check if catA and catB are in columns
        if catA not in columns or catB not in columns:
            print(f"One of the categories in comparison {catA} or {catB} is not in the columns, skipping plot.")
            # convert them to int and check if they are in columns
            try:                
                catA_int = int(catA)
                catB_int = int(catB)
                if catA_int not in columns or catB_int not in columns:
                    print(f"One of the categories in comparison {catA} or {catB} is not in the columns, skipping plot.")
                    continue
                else:
                    catA = catA_int
                    catB = catB_int
            except ValueError:
                print(f"Categories {catA} or {catB} are not valid integers, skipping plot.")
                continue
        data_to_plot = df_stat[[catA, catB]]
        data_to_plot = data_to_plot.melt(var_name="Category", value_name="Value")
        plt.figure(figsize=(6, 4))
        sns.barplot(x="Category", y="Value", data=data_to_plot, ci="sd")
        plt.title(f"{comp} (p-value (Holm-Bonferroni corrected): {p_corr:.4f})")
        # plt.savefig(output_path.replace(".txt", f"_{comp.replace(' ', '_')}_barplot.png"))
        plt.show()
        plt.close()
# %%
# # ANOVA‑like test (all categories vs. zero in CLR)
# p_list = []
# for i, cat in enumerate(X_clr.T):
#     t, p = stats.ttest_1samp(cat, 0)
#     p_list.append(p)

# # Corection FDR
# _, p_corr, _, _ = multipletests(p_list, method="fdr_bh")

# print("CLR‑ANOVA‑like test (all categories vs. zero):")
# for k, cat in enumerate(columns):
#     print(f"{cat:12} | p = {p_list[k]:.4f} | p_fdr = {p_corr[k]:.4f}")

# # Pair comparison of all pairs
# categories = columns
# tests = []
# pvals = []

# for i in range(len(categories)):
#     for j in range(i + 1, len(categories)):
#         catA = categories[i]
#         catB = categories[j]
#         A = df_stat[catA].values
#         B = df_stat[catB].values
#         diff = A - B
#         _, pval = stats.ttest_1samp(diff, 0)   # nebo stats.wilcoxon(diff) pro Wilcoxon
#         tests.append(f"{catA} vs. {catB}")
#         pvals.append(pval)

# # Korekce na mnohonásobná porovnání
# _, pvals_corrected, _, _ = multipletests(pvals, method="fdr_bh")

# print("\All pair-wise comparisons (FDR corrected):")
# for t, p, pc in zip(tests, pvals, pvals_corrected):
#     print(f"{t:20} | p = {p:.4f} | p_fdr = {pc:.4f}")

# # print pairwise comparison results in a table
# results_df = pd.DataFrame({
#     "Comparison": tests,
#     "p-value": pvals,
#     "p-value (FDR)": pvals_corrected
# })
# print("\nPair-wise comparison results:")
# print(results_df)


# %%
