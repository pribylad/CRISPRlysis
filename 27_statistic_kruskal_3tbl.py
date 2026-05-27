# %%
# Load table from where the statistic should be done
import pandas as pd
import scipy
from scipy import stats
import seaborn as sns
import os
import matplotlib.pyplot as plt
# import statsmodels
import scikit_posthocs as sp
# print version of scikit_posthocs
print("scikit_posthocs version:", sp.__version__)
# print version of scipy
print("scipy version:", scipy.__version__)
# %%
output_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/00_statistic/"
# INFO cluters are numbered based on the manuscript!, 0=1, 1=3, 2=2
# df_cluster_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_cluster_2 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_cluster_3 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "Clusters_9DELcategories_posthoc_results"
# dict_rename = {'upstream_ambiguous':5, 'bidirectional_ambiguous':1, 'downstream_ambiguous':3,
#                'upstream_NO_ambi':6, 'bidirectional_NO_ambi':2, 'downstream_NO_ambi':4,
#                'cannot_be_determined_ambiguous':7, 'DEL_exceeds_analysed_window':8,
#                'gap_outside_of_predicted_cut':9}

df_cluster_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
df_cluster_2 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
df_cluster_3 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
output_name = "Clusters_participating_end_posthoc_results"
dict_rename = {"3' NOT blunt": 30, "5' possible blunt": 51, "5' NOT blunt": 50,
               "eider 3' or 5' possible blunt": 351, "eider 3' or 5' NOT blunt": 350}

# load the three tables into pandas dataframes
df_cluster_1 = pd.read_csv(df_cluster_1, sep="\t")
df_cluster_2 = pd.read_csv(df_cluster_2, sep="\t")
df_cluster_3 = pd.read_csv(df_cluster_3, sep="\t")
# %%
df_cluster_1.head()
# %%
# print shape of the three dataframes
print("Shape of Cluster 1:", df_cluster_1.shape)
print("No of unique samples in Cluster 1:", len(df_cluster_1["variable"])/df_cluster_1["variable"].nunique())
print('')
print("Shape of Cluster 2:", df_cluster_2.shape)
print("No of unique samples in Cluster 2:", len(df_cluster_2["variable"])/df_cluster_2["variable"].nunique())
print('')
print("Shape of Cluster 3:", df_cluster_3.shape)
print("No of unique samples in Cluster 3:", len(df_cluster_3["variable"])/df_cluster_3["variable"].nunique())
print( df_cluster_1.columns)
# %%
# rename values in column "variable"
print(df_cluster_1["variable"].unique())

# rename values in column "variable" in all three dataframes
df_cluster_1["variable"] = df_cluster_1["variable"].replace(dict_rename)
df_cluster_2["variable"] = df_cluster_2["variable"].replace(dict_rename)
df_cluster_3["variable"] = df_cluster_3["variable"].replace(dict_rename)
print(df_cluster_1["variable"].unique())
# %%
# add column "cluster" to each dataframe with the name of the cluster
df_cluster_1["cluster"] = 1
df_cluster_2["cluster"] = 2
df_cluster_3["cluster"] = 3
# values in "variable" join with values in "cluster" to create new column "variable_cluster" and separeate by '_'
df_cluster_1["variable_cluster"] = df_cluster_1["variable"].astype(str) + "_" + df_cluster_1["cluster"].astype(str)
df_cluster_2["variable_cluster"] = df_cluster_2["variable"].astype(str) + "_" + df_cluster_2["cluster"].astype(str)
df_cluster_3["variable_cluster"] = df_cluster_3["variable"].astype(str) + "_" + df_cluster_3["cluster"].astype(str)
print(df_cluster_1["variable_cluster"].unique())
# %%
# merge the three dataframes into one dataframe
df_all = pd.concat([df_cluster_1, df_cluster_2, df_cluster_3], ignore_index=True)
print(df_all.shape)
print(df_all["variable_cluster"].unique())
# %%
# test normality of individual "variable_cluster" groups with Shapiro-Wilk test
for group in df_all["variable_cluster"].unique():
    stat, p = stats.shapiro(df_all[df_all["variable_cluster"] == group]["value"])
    print("Group:", group, "Shapiro-Wilk test statistic:", stat, "p-value:", p)
    # Print if the group is normally distributed or not based on p-value
    if p > 0.05:
        print("Group:", group, "is normally distributed")
    else:
        print("Group:", group, "is NOT normally distributed")
# %%
# Do Kruskal-Wallis test to compare groups in "variable_cluster" (data did not pass normality test)
stat, p = stats.kruskal(*[df_all[df_all["variable_cluster"] == group]["value"] for group in df_all["variable_cluster"].unique()])
print("Kruskal-Wallis test statistic:", stat, "p-value:", p)
if p < 0.05:
    print("There is a significant difference between groups in variable_cluster, medians were compared.")
else:
    print("There is NO significant difference between groups in variable_cluster, medians were compared.")
# %%
# Do posthoc test to compare groups in "variable_cluster" with Dunn's test and adjust p-values with Holm-Bonferroni method
# compare only groups that share the first part of the name, e.g. 1_1, 1_2, 1_3, etc. and do not compare groups that do not share the first part of the name, e.g. 1_1 and 2_1
# create list of groups to compare based on the first part of the name
groups_to_compare = {}
for group in df_all["variable_cluster"].unique():
    group_name = group.split("_")[0]
    if group_name not in groups_to_compare:
        groups_to_compare[group_name] = []
    groups_to_compare[group_name].append(group)

# sort groups_to_compare by group name
groups_to_compare = dict(sorted(groups_to_compare.items()))
groups_to_compare
# %%
i=0
for group_name, groups in groups_to_compare.items():
    print("Comparing:", group_name, groups)

    df_filtered = df_all[df_all["variable_cluster"].isin(groups)]

    p_values = sp.posthoc_dunn(
        df_filtered,
        val_col="value",
        group_col="variable_cluster",
        p_adjust="holm"
    )
    # make another table where if p<0.0001, value is "****", if p<0.001, value is "***", if p<0.01, value is "**", if p<0.05, value is "*", else value is "ns"
    p_values_significance = p_values.copy()
    p_values_significance[(p_values < 0.001)] = "***"
    p_values_significance[(p_values >= 0.001) & (p_values < 0.01)] = "**"
    p_values_significance[(p_values >= 0.01) & (p_values < 0.05)] = "*"
    p_values_significance[p_values >= 0.05] = "ns"

    print(p_values)
    print(p_values_significance)

    # save p_values  and p_values_significance to txt file with name "Clusters_9DELcategories_posthoc_results_group_name.txt" in output_dir and separate by tab and each iteration add to the same file and separate by two empty rows
    # if txt file does not exist, create
    if i == 0:
        with open(output_dir + output_name + ".txt", "w") as f:
            p_values.to_csv(f, sep="\t")
            p_values_significance.to_csv(f, sep="\t")
            f.write("\n\n")
    else:
        with open(output_dir + output_name + ".txt", "a") as f:
            p_values.to_csv(f, sep="\t")
            p_values_significance.to_csv(f, sep="\t")
            f.write("\n\n")
    # make heatmap of p_values with seaborn and save it as "Clusters_9DELcategories_posthoc_heatmap_group_name.png" in output_dir

    plt.figure(figsize=(1, 1))
    # 
    sns.heatmap(p_values, annot=False, cmap="coolwarm", cbar=True, vmin=0, vmax=1)
    plt.title("Normality not passed Shapiro-Wilk test p < 0.05\n Kruskal-Wallis test p < 0.05\nDunn's posthoc test results for variable_cluster groups\np-values adjusted with Holm-Bonferroni method")
    plt.savefig(output_dir + output_name + "_" + group_name + ".png", dpi=300)
    plt.show()
    i+=1

# %%

