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
# INFO comparing plants vs animals
# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_4_DEL_orientation_summary_violinplot_small.tsv"
# df_2 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "Kingdoms_9DELcategories_posthoc_results"
# dict_rename = {'upstream_ambiguous':5, 'bidirectional_ambiguous':1, 'downstream_ambiguous':3,
#                'upstream_NO_ambi':6, 'bidirectional_NO_ambi':2, 'downstream_NO_ambi':4,
#                'cannot_be_determined_ambiguous':7, 'DEL_exceeds_analysed_window':8,
#                'gap_outside_of_predicted_cut':9}

# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
# df_2 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/001_1_INDELSUB_summary_violinplot_260518.tsv"
# output_name = "Kingdoms_INDELSUB_posthoc_results"
# dict_rename = {}

df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
df_2 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
output_name = "Kingdoms_IN_posthoc_results"
dict_rename = {}



# load the three tables into pandas dataframes
df_1 = pd.read_csv(df_1, sep="\t")
df_2 = pd.read_csv(df_2, sep="\t")

if "Event_around_cut" in df_1.columns:
    df_1 = df_1.rename(columns={"Event_around_cut": "variable"})
if "Event_around_cut" in df_2.columns:
    df_2 = df_2.rename(columns={"Event_around_cut": "variable"})
# %%
df_1.head()
# %%
# print shape of the three dataframes
print("Shape of df 1:", df_1.shape)
print("No of unique samples in df 1:", len(df_1["variable"])/df_1["variable"].nunique())
print('')
print("Shape of df 2:", df_2.shape)
print("No of unique samples in df 2:", len(df_2["variable"])/df_2["variable"].nunique())
print('')
print( df_1.columns)
# %%
# rename values in column "variable"
print(df_1["variable"].unique())

# rename values in column "variable" in all three dataframes
df_1["variable"] = df_1["variable"].replace(dict_rename)
df_2["variable"] = df_2["variable"].replace(dict_rename)
print(df_1["variable"].unique())
# %%
# add column "kingdom" to each dataframe with the name of the kingdom
df_1["kingdom"] = 1
df_2["kingdom"] = 2
# values in "variable" join with values in "kingdom" to create new column "variable_kingdom" and separeate by '_'
df_1["variable_kingdom"] = df_1["variable"].astype(str) + "_" + df_1["kingdom"].astype(str)
df_2["variable_kingdom"] = df_2["variable"].astype(str) + "_" + df_2["kingdom"].astype(str)
print(df_1["variable_kingdom"].unique())
# %%
# merge the three dataframes into one dataframe
df_all = pd.concat([df_1, df_2], ignore_index=True)
print(df_all.shape)
print(df_all["variable_kingdom"].unique())
# %%
# test normality of individual "variable_kingdom" groups with Shapiro-Wilk test
for group in df_all["variable_kingdom"].unique():
    stat, p = stats.shapiro(df_all[df_all["variable_kingdom"] == group]["value"])
    print("Group:", group, "Shapiro-Wilk test statistic:", stat, "p-value:", p)
    # Print if the group is normally distributed or not based on p-value
    if p > 0.05:
        print("Group:", group, "is normally distributed")
    else:
        print("Group:", group, "is NOT normally distributed")
# %%
# Do Kruskal-Wallis test to compare groups in "variable_kingdom" (data did not pass normality test)
stat, p = stats.kruskal(*[df_all[df_all["variable_kingdom"] == group]["value"] for group in df_all["variable_kingdom"].unique()])
print("Kruskal-Wallis test statistic:", stat, "p-value:", p)
if p < 0.05:
    print("There is a significant difference between groups in variable_kingdom, medians were compared.")
else:
    print("There is NO significant difference between groups in variable_kingdom, medians were compared.")
# %%
# Do posthoc test to compare groups in "variable_kingdom" with Dunn's test and adjust p-values with Holm-Bonferroni method
# compare only groups that share the first part of the name, e.g. 1_1, 1_2, 1_3, etc. and do not compare groups that do not share the first part of the name, e.g. 1_1 and 2_1
# create list of groups to compare based on the first part of the name
groups_to_compare = {}
for group in df_all["variable_kingdom"].unique():
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

    df_filtered = df_all[df_all["variable_kingdom"].isin(groups)]

    p_values = sp.posthoc_dunn(
        df_filtered,
        val_col="value",
        group_col="variable_kingdom",
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

    # save p_values  and p_values_significance to txt file with name "kingdoms_9DELcategories_posthoc_results_group_name.txt" in output_dir and separate by tab and each iteration add to the same file and separate by two empty rows
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
    # make heatmap of p_values with seaborn and save it as "kingdoms_9DELcategories_posthoc_heatmap_group_name.png" in output_dir

    plt.figure(figsize=(1, 1))
    # 
    sns.heatmap(p_values, annot=False, cmap="coolwarm", cbar=True, vmin=0, vmax=1)
    plt.title("Normality not passed Shapiro-Wilk test p < 0.05\n Kruskal-Wallis test p < 0.05\nDunn's posthoc test results for variable_kingdom groups\np-values adjusted with Holm-Bonferroni method")
    plt.savefig(output_dir + output_name + "_" + group_name + ".png", dpi=300)
    plt.show()
    i+=1

# %%

