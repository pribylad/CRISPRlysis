
# %%
"""
Compositional data testing pipeline (revised):
  1. Pivot long -> wide, unique ID from organism/cell_line/mutant/gene/gRNA_seq/hpi
  2. Half-minimum zero imputation
  3. CLR (Centered Log-Ratio) transform
  4. Shapiro-Wilk normality test per component (diagnostic)
  5. Skewness test per pair → paired t-test (symmetric) or Wilcoxon (skewed)
  6. FDR correction (Benjamini-Hochberg)
  7. Compact Letter Display (CLD) violin plot

Requires: pandas, numpy, scipy, statsmodels, matplotlib, seaborn
  pip install statsmodels
"""

import sys
import string
import numpy as np
import pandas as pd
import scipy
from scipy.stats import shapiro, ttest_1samp, wilcoxon, skewtest
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore")

# ── 0. Load long-format data ──────────────────────────────────────────────────
debugging = True

# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "plants_9DELcategories_posthoc_results"
# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "animals_9DELcategories_posthoc_results"
df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
output_name = "Cluster_0_A_9DELcategories_posthoc_results" # it is cluster A in the paper
# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "Cluster_2_B_9DELcategories_posthoc_results" # it is cluster B in the paper
# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "Cluster_1_C_9DELcategories_posthoc_results" # it is cluster C in the paper
ID_COLS = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]
dict_rename = {
    'upstream_ambiguous': 5,
    'bidirectional_ambiguous': 1,
    'downstream_ambiguous': 3,
    'upstream_NO_ambi': 6,
    'bidirectional_NO_ambi': 2,
    'downstream_NO_ambi': 4,
    'cannot_be_determined_ambiguous': 7,
    'DEL_exceeds_analysed_window': 8,
    'gap_outside_of_predicted_cut': 9,
}

# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
# output_name = "plants_INDELSUB_posthoc_results"
# # df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/001_1_INDELSUB_summary_violinplot_260518.tsv"
# # output_name = "animals_INDELSUB_posthoc_results"
# ID_COLS = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]
# dict_rename = {}

# # df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# # output_name = "plants_IN_posthoc_results"
# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# output_name = "animals_IN_posthoc_results"
# ID_COLS = ["ID_no_rep"]
# dict_rename = {}


# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_0_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# output_name = "Cluster_0_A_DEL_end_posthoc_results" # it is cluster A in the paper
# # df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_2_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# # output_name = "Cluster_2_B_DEL_end_posthoc_results" # it is cluster B in the paper
# # df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/07_cluster_1_DEL_03_filtered_efficacy_samples_limit5efficacy_All_260518/003_5_DEL_end_summary_norm_violinplot_small_260522.tsv"
# # output_name = "Cluster_1_C_DEL_end_posthoc_results" # it is cluster C in the paper
# ID_COLS = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]
# dict_rename = {}

output_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/00_statistic/"


class Tee:
    """Write to both stdout and a file simultaneously."""
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()


_log_file = open(output_dir + output_name + ".txt", "w", encoding="utf-8")
sys.stdout = Tee(sys.__stdout__, _log_file)

df_long = pd.read_csv(df_dir, sep="\t")
print(df_long.columns)

if "Event_around_cut" in df_long.columns:
    df_long = df_long.rename(columns={"Event_around_cut": "variable"})

# ── 1. Build unique sample ID and pivot to wide ───────────────────────────────
df_long["sample_id"] = df_long[ID_COLS].astype(str).agg(" - ".join, axis=1)

if df_long["sample_id"].nunique()*df_long['variable'].nunique() != len(df_long):
    raise ValueError("Duplicate sample_id values. Ensure ID_COLS produce unique identifiers.")

df_wide = df_long.pivot_table(
    index="sample_id",
    columns="variable",
    values="value",
)
df_wide.columns.name = None
df_wide.index.name = "sample_id"

if debugging:
    df_wide.to_csv(output_dir + output_name + "_01_wide.tsv", sep="\t")

components = list(df_wide.columns)
X_raw = df_wide.values.astype(float)


print("=" * 65)
print("COMPOSITIONAL DATA ANALYSIS — ADAPTIVE PAIRED TEST ON CLR DIFFERENCES")
print("=" * 65)
print(f"\nSamples    : {len(df_wide)}")
print(f"Variables  : {len(components)}")
for c in components:
    print(f"  - {c} - {dict_rename.get(c, c)}")

# ── 2. Half-minimum zero imputation ──────────────────────────────────────────
X = X_raw.copy()

if np.isnan(X).any():
    raise ValueError("NaN values in data. Handle before running this analysis.")
if (X < 0).any():
    raise ValueError("Negative values in data. Compositional data cannot be negative.")

nonzero_vals = X[X > 0]
half_min = nonzero_vals.min() / 2 if len(nonzero_vals) > 0 else 1e-6
n_zeros = (X == 0).sum()
X[X == 0] = half_min

print(f"\nZero imputation : {n_zeros} zero cells → half-minimum = {half_min:.6g}")

# ── 3. CLR transform ──────────────────────────────────────────────────────────
# Re-normalise rows to sum 1 (half-min imputation shifts the sums slightly)
X = X / X.sum(axis=1, keepdims=True)
log_X = np.log(X)
X_clr = log_X - log_X.mean(axis=1, keepdims=True)

clr_df = pd.DataFrame(X_clr, columns=components, index=df_wide.index)

if debugging:
    clr_df.to_csv(output_dir + output_name + "_02_clr.tsv", sep="\t")

print(f"CLR transform   : row sums ≈ 0 (max abs = {np.abs(X_clr.sum(axis=1)).max():.2e})  ✓")

# Sort components by dict_rename values for display
components_sorted = sorted(components, key=lambda c: dict_rename.get(c, float("inf")))


# ── 4. Shapiro-Wilk normality test (diagnostic) ───────────────────────────────
print("\n" + "=" * 65)
print("STEP 1: Shapiro-Wilk Normality Test (CLR values, per variable)")
print("-" * 65)
print(f"  {'Variable':<36}  {'W':>7}  {'p-value':>8}  {'Normal?':>8}")
print(f"  {'-'*36}  {'-'*7}  {'-'*8}  {'-'*8}")
print(f"  scipy.stats.shapiro  — scipy {scipy.__version__}")
print("-" * 65)

any_nonnormal = False
for comp in components_sorted:
    vals = clr_df[comp].values
    if len(vals) >= 3:
        W, p = shapiro(vals)
        normal = "yes" if p >= 0.05 else "NO *"
        if p < 0.05:
            any_nonnormal = True
        print(f"  {comp:<36}  {W:>7.4f}  {p:>8.1e}  {normal:>8}")
    else:
        print(f"  {comp:<36}  {'n<3':>7}  {'n<3':>8}  {'n/a':>8}")

print(f"\n  * p < 0.05 = non-normal distribution of CLR values.")
if any_nonnormal:
    print(f"  → Non-normality detected. With n = {len(df_wide)}, the CLT justifies")
    print(f"    proceeding with the paired t-test (means are asymptotically normal).")
    print(f"    Skewed difference distributions will additionally trigger Wilcoxon for ALL pairs.")
else:
    print("  → All CLR distributions consistent with normality — t-test valid.")


# CLR distribution histograms (visual inspection)
print("\n" + "=" * 65)
print("STEP 1b: CLR distribution histograms")
print("-" * 65)

n_cols = 3
n_rows = (len(components) + n_cols - 1) // n_cols
plt.figure(figsize=(n_cols * 5, n_rows * 4))
for i, comp in enumerate(components_sorted):
    plt.subplot(n_rows, n_cols, i + 1)
    sns.histplot(clr_df[comp], kde=True)
    rename_val = dict_rename.get(comp, None)
    title = f"[{rename_val}] - {comp}" if rename_val is not None else comp
    plt.title(title)
    plt.xlabel("CLR value")
    plt.ylabel("Frequency")
plt.tight_layout()
if debugging:
    plt.savefig(output_dir + output_name + "_03_clr_histograms.png", dpi=200)
plt.show()


# ── 5. Global skewness screening → single test method for all pairs ───────────
#
# Pre-pass: compute CLR differences for every pair and run D'Agostino skewtest.
# If ANY pair's differences are significantly skewed (p < 0.05), the entire
# analysis uses Wilcoxon signed-rank for all pairs (conservative, consistent).
# Otherwise all pairs use the paired t-test.

print("\n" + "=" * 65)
print("STEP 2a: Global skewness screening across all pairs")
print("-" * 65)
print(f"  scipy.stats.skewtest — scipy {scipy.__version__}")
print("-" * 65)

pair_list = list(combinations(components_sorted, 2))

skew_report = []
any_skewed = False
for catA, catB in pair_list:
    diff = clr_df[catA].values - clr_df[catB].values
    if len(diff) > 7:
        _, p_skew = skewtest(diff)
        skew_val = float(pd.Series(diff).skew())
        is_skewed = p_skew < 0.05
    else:
        p_skew = np.nan
        skew_val = np.nan
        is_skewed = False
    if is_skewed:
        any_skewed = True
    skew_report.append({
        "catA": catA, "catB": catB,
        "skewness": skew_val, "p_skewtest": p_skew, "is_skewed": is_skewed,
    })

n_skewed_pairs = sum(r["is_skewed"] for r in skew_report)
global_test = "Wilcoxon" if any_skewed else "t-test"

print(f"\n  Pairs screened         : {len(pair_list)}")
print(f"  Pairs with skewed diffs: {n_skewed_pairs}")
print(f"\n  {'→ At least one skewed pair detected.' if any_skewed else '→ No skewed pairs detected.'}")
print(f"  Global test chosen     : {global_test} (applied to ALL {len(pair_list)} pairs)")

# Print per-pair skewness detail
print(f"\n  {'catA':<36}  {'catB':<36}  {'skewness':>8}  {'p_skew':>8}  {'skewed?':>8}")
print(f"  {'-'*36}  {'-'*36}  {'-'*8}  {'-'*8}  {'-'*8}")
for r in skew_report:
    sk_str  = f"{r['skewness']:>8.3f}" if not np.isnan(r["skewness"]) else "     n/a"
    skp_str = f"{r['p_skewtest']:>8.2e}" if not np.isnan(r["p_skewtest"]) else "     n/a"
    flag    = "YES *" if r["is_skewed"] else "no"
    print(f"  {r['catA']:<36}  {r['catB']:<36}  {sk_str}  {skp_str}  {flag:>8}")


# ── 5b. Run chosen test on all pairs ─────────────────────────────────────────
print("\n" + "=" * 65)
print(f"STEP 2b: {global_test} on all {len(pair_list)} pairs")
print("-" * 65)
print(f"  {'Wilcoxon signed-rank (scipy.stats.wilcoxon)' if global_test == 'Wilcoxon' else 'One-sample t-test (scipy.stats.ttest_1samp)'}")
print(f"  scipy {scipy.__version__}  |  {len(df_wide)} paired observations per comparison")
print("-" * 65)

t_results = []
for sr in skew_report:
    catA, catB = sr["catA"], sr["catB"]
    diff = clr_df[catA].values - clr_df[catB].values

    if global_test == "Wilcoxon":
        stat, p_raw = wilcoxon(diff, alternative="two-sided", method="auto")
    else:
        stat, p_raw = ttest_1samp(diff, popmean=0)

    t_results.append({
        "catA": catA,
        "catB": catB,
        "n": len(diff),
        "mean_diff": diff.mean(),
        "median_diff": float(np.median(diff)),
        "skewness": sr["skewness"],
        "p_skewtest": sr["p_skewtest"],
        "is_skewed": sr["is_skewed"],
        "test_used": global_test,
        "statistic": stat,
        "p_raw": p_raw,
    })

# ── 6. Benjamini-Hochberg FDR correction ──────────────────────────────────────
p_vals = [r["p_raw"] for r in t_results]
reject, p_adj, _, _ = multipletests(p_vals, method="fdr_bh")
for i, r in enumerate(t_results):
    r["p_adj"] = p_adj[i]
    r["reject"] = bool(reject[i])

results_df = pd.DataFrame(t_results)
if debugging:
    results_df.to_csv(output_dir + output_name + "_04_test_results.tsv", sep="\t", index=False)

# ── Print full pairwise results table ────────────────────────────────────────
print(
    f"\n  {'catA':<36}  {'catB':<36}  {'mean diff':>9}  "
    f"{'stat':>7}  {'p_raw':>8}  {'p_BH':>8}  {'sig':>5}"
)
print(
    f"  {'-'*36}  {'-'*36}  {'-'*9}  "
    f"{'-'*7}  {'-'*8}  {'-'*8}  {'-'*5}"
)
for r in t_results:
    stars = (
        "***" if r["p_adj"] < 0.001
        else ("**" if r["p_adj"] < 0.01
              else ("*" if r["p_adj"] < 0.05 else "ns"))
    )
    print(
        f"  {r['catA']:<36}  {r['catB']:<36}  "
        f"{r['mean_diff']:>9.4f}  {r['statistic']:>7.3f}  "
        f"{r['p_raw']:>8.2e}  {r['p_adj']:>8.2e}  {stars:>5}"
    )

# Significant pairs
print(f"\n  Significant pairs (BH-FDR p < 0.05)  [{global_test}]:")
sig_pairs_found = False
for r in t_results:
    if r["reject"]:
        stars = "***" if r["p_adj"] < 0.001 else ("**" if r["p_adj"] < 0.01 else "*")
        rA = dict_rename.get(r["catA"], r["catA"])
        rB = dict_rename.get(r["catB"], r["catB"])
        print(f"    [{rA}] {r['catA']}  vs  [{rB}] {r['catB']}  →  p_BH = {r['p_adj']:.2e} {stars}")
        sig_pairs_found = True
if not sig_pairs_found:
    print("    None after BH-FDR correction.")

# Non-significant pairs
print(f"\n  Non-significant pairs (BH-FDR p ≥ 0.05)  [{global_test}]:")
ns_found = False
for r in t_results:
    if not r["reject"]:
        rA = dict_rename.get(r["catA"], r["catA"])
        rB = dict_rename.get(r["catB"], r["catB"])
        print(f"    [{rA}] {r['catA']}  vs  [{rB}] {r['catB']}  →  p_BH = {r['p_adj']:.2e}")
        ns_found = True
if not ns_found:
    print("    All pairs are significantly different after BH-FDR correction.")


# ── 7. CLD violin plot ────────────────────────────────────────────────────────

def compute_cld(categories, p_matrix, alpha=0.05):
    """
    Compact Letter Display — Piepho (2004) algorithm.

    Two categories share a letter iff BH-adjusted p >= alpha (not sig diff).

    Parameters
    ----------
    categories : list of str  — category names in display order
    p_matrix   : pd.DataFrame — symmetric adjusted p-value matrix
    alpha      : float        — significance threshold

    Returns
    -------
    dict : {category_name: letter_string}
    """
    cats = list(categories)

    def sig(a, b):
        return p_matrix.loc[a, b] < alpha

    groups = {t: {0} for t in cats}
    letter_groups = {0: set(cats)}
    next_group = [1]

    def split_group(letter, keep_anchor, split_anchor):
        members = set(letter_groups[letter])
        new_members = {t for t in members if not sig(t, split_anchor) or t == split_anchor}
        new_members.discard(keep_anchor)
        old_members = members - new_members
        if not new_members or old_members == members:
            return
        letter_groups[letter] = old_members
        for t in members - old_members:
            groups[t].discard(letter)
        nl = next_group[0]; next_group[0] += 1
        letter_groups[nl] = new_members
        for t in new_members:
            groups[t].add(nl)
        for t in members - old_members - new_members:
            placed = False
            for el, eg in letter_groups.items():
                if all(not sig(t, x) for x in eg):
                    eg.add(t); groups[t].add(el); placed = True; break
            if not placed:
                sl = next_group[0]; next_group[0] += 1
                letter_groups[sl] = {t}; groups[t].add(sl)

    changed = True
    while changed:
        changed = False
        for i, catA in enumerate(cats):
            for catB in cats[i + 1:]:
                if not sig(catA, catB):
                    continue
                for letter in list(groups[catA] & groups[catB]):
                    before = {k: frozenset(v) for k, v in letter_groups.items()}
                    split_group(letter, catA, catB)
                    if {k: frozenset(v) for k, v in letter_groups.items()} != before:
                        changed = True

    # Remove redundant letters
    for t in cats:
        for l in list(groups[t]):
            if len(groups[t]) <= 1:
                break
            can_remove = all(
                (groups[t] - {l}) & groups[t2]
                for t2 in cats if t2 != t and not sig(t, t2)
            )
            if can_remove:
                groups[t].discard(l)
                letter_groups[l].discard(t)

    all_used = sorted({l for ls in groups.values() for l in ls})
    i2l = {idx: string.ascii_lowercase[i] for i, idx in enumerate(all_used) if i < 26}
    return {t: "".join(i2l[l] for l in sorted(groups[t]) if l in i2l) for t in cats}


# Build symmetric adjusted p-value matrix
p_matrix = pd.DataFrame(
    np.ones((len(components_sorted), len(components_sorted))),
    index=components_sorted,
    columns=components_sorted,
)
for r in t_results:
    p_matrix.loc[r["catA"], r["catB"]] = r["p_adj"]
    p_matrix.loc[r["catB"], r["catA"]] = r["p_adj"]

cld_letters = compute_cld(components_sorted, p_matrix, alpha=0.05)

print("\n" + "=" * 65)
print("STEP 3: Compact Letter Display (CLD)")
print("-" * 65)
print("  Categories sharing a letter are NOT significantly different (BH p ≥ 0.05)\n")
for comp in components_sorted:
    rv = dict_rename.get(comp, "")
    label = f"[{rv}]" if rv != "" else "    "
    print(f"  {label}  {comp:<36}  →  '{cld_letters[comp]}'")

# ── Violin plot with CLD letters ──────────────────────────────────────────────
# set figure size with no. of columns based on the number of 
sns.set(style="whitegrid")
fig, ax = plt.subplots(figsize=(len(components_sorted), 7))

x_labels = [
    f"[{dict_rename[c]}]\n{c}" if c in dict_rename else c
    for c in components_sorted
]

plot_data = clr_df[components_sorted].copy()
plot_data.columns = x_labels
plot_long = plot_data.melt(var_name="Category", value_name="CLR value")

sns.violinplot(
    data=plot_long,
    x="Category",
    y="CLR value",
    order=x_labels,
    inner="box",
    ax=ax,
    palette="Set2",
)

y_top = clr_df[components_sorted].max().max()
letter_y = y_top + 0.20 * abs(y_top)
for i, (comp, label) in enumerate(zip(components_sorted, x_labels)):
    ax.text(
        i, letter_y, cld_letters[comp],
        ha="center", va="bottom", fontsize=13, fontweight="bold",
    )

ax.axhline(0, color="grey", lw=0.8, ls="--", alpha=0.6)
ax.set_title(
    f"CLR-transformed deletion orientation categories\n"
    f"(n = {len(df_wide)} samples; {global_test} + BH-FDR CLD, α = 0.05)",
    fontsize=13,
)
ax.set_xlabel("Deletion orientation category", fontsize=11)
ax.set_ylabel("CLR value", fontsize=11)
plt.xticks(rotation=15, ha="right", fontsize=9)
plt.tight_layout()

if debugging:
    plt.savefig(output_dir + output_name + "_05_violin_CLD.png", dpi=200)
plt.show()


print("\n" + "=" * 65)
print("INTERPRETATION")
print("=" * 65)
print(f"""
  Data:
    {len(df_wide)} biological samples × {len(components)} deletion orientation categories.
    Categories are compositional (sum to 100% per sample).

  Zero replacement (half-minimum imputation):
    {n_zeros} zero cells replaced with half the global minimum non-zero
    value ({half_min:.6g}). Scales with data magnitude (Martín-Fernández).

  CLR transformation:
    log(xᵢ / geometric_mean(x)) maps compositions into unconstrained
    real space. Row sums = 0 confirms correct transformation.

  Shapiro-Wilk (diagnostic):
    Checks per-category CLR normality. With n = {len(df_wide)}, CLT justifies
    the t-test even if some distributions deviate from normality.

  Global skewness screening:
    D'Agostino skewtest applied to the CLR difference distribution of every
    pair ({len(pair_list)} pairs total). If ANY pair is significantly skewed
    (p_skew < 0.05), Wilcoxon signed-rank is used for ALL pairs to ensure
    a consistent, conservative analysis.
    Skewed pairs found: {n_skewed_pairs} / {len(pair_list)}
    → Test used for all pairs: {global_test}

  {global_test}:
    {'Wilcoxon signed-rank: tests whether the median CLR difference equals zero. Non-parametric; robust to asymmetric difference distributions. Preserves within-sample pairing.' if global_test == 'Wilcoxon' else 'One-sample t-test: tests whether the mean CLR difference equals zero. More powerful under near-normality. Preserves within-sample pairing.'}

  BH-FDR correction:
    Benjamini-Hochberg applied to all {len(pair_list)} p-values.
    More powerful than Holm/FWER; standard in biological studies.

  CLD violin plot:
    Violins show the full CLR distribution per category.
    Letters above each violin: shared letter = NOT significantly
    different (BH-adjusted p ≥ 0.05).
""")

sys.stdout = sys.__stdout__
_log_file.close()
print("Done. Results saved to:", output_dir)