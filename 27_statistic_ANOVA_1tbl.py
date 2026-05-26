# %%
# Load table from where the statistic should be done
import pandas as pd
import scipy
from scipy import stats
import seaborn as sns
import os
import numpy as np
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
df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_4_DEL_orientation_summary_violinplot_small.tsv"
output_name = "plants_9DELcategories_posthoc_results"
# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "Kingdoms_9DELcategories_posthoc_results"
dict_rename = {'upstream_ambiguous':5, 'bidirectional_ambiguous':1, 'downstream_ambiguous':3,
               'upstream_NO_ambi':6, 'bidirectional_NO_ambi':2, 'downstream_NO_ambi':4,
               'cannot_be_determined_ambiguous':7, 'DEL_exceeds_analysed_window':8,
               'gap_outside_of_predicted_cut':9}

# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/001_1_INDELSUB_summary_violinplot_260518.tsv"
# output_name = "Kingdoms_INDELSUB_posthoc_results"
# dict_rename = {}

# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# df_1 = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# output_name = "Kingdoms_IN_posthoc_results"
# dict_rename = {}



# load the table into pandas dataframes
df_1 = pd.read_csv(df_1, sep="\t")

if "Event_around_cut" in df_1.columns:
    df_1 = df_1.rename(columns={"Event_around_cut": "variable"})
# %%
df_1.head()
# %%
# print shape of the dataframe
print("Shape of df 1:", df_1.shape)
print("No of unique samples in df 1:", len(df_1["variable"])/df_1["variable"].nunique())
print('')
print(df_1.columns)
# %%
# rename values in column "variable"
print(df_1["variable"].unique())

# rename values in column "variable"
df_1["variable"] = df_1["variable"].replace(dict_rename)
# sort table by column "variable"
df_1 = df_1.sort_values(by="variable")
print(df_1["variable"].unique())
# %%
# make CLR transformation of "value" column for each group in "variable" column
# first, add 1 to all values in "value" column to avoid log of zero
df_1["value"] = df_1["value"] + 1
# then, for each group in "variable" column, divide "value" column by the sum of "value" column for that group and take log of the result
df_1["value"] = df_1.groupby("variable")["value"].transform(lambda x: x / x.sum())
df_1["value"] = df_1["value"].apply(lambda x: np.log(x))
# %%
df_1
# %%
# test normality of individual "variable" groups with Shapiro-Wilk test
for group in df_1["variable"].unique():
    stat, p = stats.shapiro(df_1[df_1["variable"] == group]["value"])
    print("Group:", group, "Shapiro-Wilk test statistic:", stat, "p-value:", p)
    # Print if the group is normally distributed or not based on p-value
    if p > 0.05:
        print("Group:", group, "is normally distributed")
    else:
        print("Group:", group, "is NOT normally distributed")
# %%
# Do Kruskal-Wallis test to compare groups in "variable" (data did not pass normality test)
stat, p = stats.kruskal(*[df_1[df_1["variable"] == group]["value"] for group in df_1["variable"].unique()])
print("Kruskal-Wallis test statistic:", stat, "p-value:", p)
if p < 0.05:
    print("There is a significant difference between groups in variable, medians were compared.")
else:
    print("There is NO significant difference between groups in variable, medians were compared.")
# %%
# Do posthoc test to compare groups in "variable" with Dunn's test and adjust p-values with Holm-Bonferroni method
# compare all groups in "variable" column with each other and save the results in a table
p_values = sp.posthoc_dunn(df_1, val_col="value", group_col="variable", p_adjust="holm")
# make another table where if p<0.0001, value is "****", if p<0.001, value is "***", if p<0.01, value is "**", if p<0.05, value is "*", else value is "ns"
p_values_significance = p_values.copy()
p_values_significance[(p_values < 0.001)] = "***"
p_values_significance[(p_values >= 0.001) & (p_values < 0.01)] = "**"
p_values_significance[(p_values >= 0.01) & (p_values < 0.05)] = "*"
p_values_significance[p_values >= 0.05] = "ns"
print(p_values)
print(p_values_significance)
# %%


# save p_values  and p_values_significance to txt file with name "kingdoms_9DELcategories_posthoc_results_group_name.txt" in output_dir and separate by tab and each iteration add to the same file and separate by two empty rows
# if txt file does not exist, create
with open(output_dir + output_name + ".txt", "w") as f:
    p_values.to_csv(f, sep="\t")
    p_values_significance.to_csv(f, sep="\t")
    f.write("\n\n")

# make heatmap of p_values with seaborn and save it as "kingdoms_9DELcategories_posthoc_heatmap_group_name.png" in output_dir

plt.figure(figsize=(5, 4))
sns.heatmap(p_values, annot=False, cmap="coolwarm", cbar=True, vmin=0, vmax=1)
plt.title("Normality not passed Shapiro-Wilk test p < 0.05\n Kruskal-Wallis test p < 0.05\nDunn's posthoc test results for variable groups\np-values adjusted with Holm-Bonferroni method")


# set tight layout
plt.tight_layout()
# save figure
plt.savefig(output_dir + output_name + ".png", dpi=300)
plt.show()


# %%

