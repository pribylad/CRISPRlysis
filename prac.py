# %%
"""
Compositional data testing pipeline:
  1. Pivot long -> wide, unique ID from organism/cell_line/mutant/gene/gRNA_seq/hpi
  2. Half-minimum zero imputation
  3. CLR (Centered Log-Ratio) transform
  4. Shapiro-Wilk normality test per component per group
  5. Kruskal-Wallis test per component
  6. Dunn post-hoc test with Holm correction per component

Requires: pandas, numpy, scipy, scikit_posthocs
  pip install scikit-posthocs
"""

import sys

import numpy as np
import pandas as pd
from scipy.stats import shapiro, kruskal
import scikit_posthocs as sp
import warnings
warnings.filterwarnings("ignore")

# ── 0. Load long-format data ──────────────────────────────────────────────────
# # df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/003_4_DEL_orientation_summary_violinplot_small.tsv"
# # output_name = "plants_9DELcategories_posthoc_results"
# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/003_4_DEL_orientation_summary_violinplot_small.tsv"
# output_name = "animals_9DELcategories_posthoc_results"
# ID_COLS = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]
# dict_rename = {'upstream_ambiguous':5, 'bidirectional_ambiguous':1, 'downstream_ambiguous':3,
#                'upstream_NO_ambi':6, 'bidirectional_NO_ambi':2, 'downstream_NO_ambi':4,
#                'cannot_be_determined_ambiguous':7, 'DEL_exceeds_analysed_window':8,
#                'gap_outside_of_predicted_cut':9}

# # df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/001_1_INDELSUB_summary_violinplot_260518.tsv"
# # output_name = "plants_INDELSUB_posthoc_results"
# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/001_1_INDELSUB_summary_violinplot_260518.tsv"
# output_name = "animals_INDELSUB_posthoc_results"
# ID_COLS = ["organism", "cell_line", "mutant", "gene", "gRNA_seq", "hours_post_inoculation"]
# dict_rename = {}

# df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_plants/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
# output_name = "plants_IN_posthoc_results"
df_dir = "/home/lefnui/CRISPRlysis_outputs/04_charts/02_filtered_efficacy_95/with_Ponce_samples/02_animals/002_1_IN_no_of_IN_nucleotides_summary_violinplot_small_260514.tsv"
output_name = "animals_IN_posthoc_results"
ID_COLS = ["ID_no_rep"]
dict_rename = {}



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

if "Event_around_cut" in df_long.columns:
    df_long = df_long.rename(columns={"Event_around_cut": "variable"})

# ── 1. Build unique sample ID and pivot to wide ───────────────────────────────


df_long["sample_id"] = df_long[ID_COLS].astype(str).agg(" - ".join, axis=1)

df_wide = df_long.pivot_table(
    index="sample_id",
    columns="variable",
    values="value",
    aggfunc="first"
)
df_wide.columns.name = None
df_wide.index.name = "sample_id"

components = list(df_wide.columns)
X_raw = df_wide.values.astype(float)

print("=" * 65)
print("COMPOSITIONAL DATA ANALYSIS — COMPARING VARIABLES")
print("=" * 65)
print(f"\nSamples    : {len(df_wide)}")
print(f"Variables  : {len(components)}")
for c in components:
    print(f"  - {c} - {dict_rename.get(c, c)}")

# ── 2. Half-minimum zero imputation ──────────────────────────────────────────
X = X_raw.copy()
X = np.nan_to_num(X, nan=0.0)
X[X < 0] = 0.0

nonzero_vals = X[X > 0]
half_min = nonzero_vals.min() / 2 if len(nonzero_vals) > 0 else 1e-6
n_zeros = (X == 0).sum()
X[X == 0] = half_min

print(f"\nZero imputation : {n_zeros} zero cells → half-minimum = {half_min:.6g}")

# ── 3. CLR transform ──────────────────────────────────────────────────────────
X = X / X.sum(axis=1, keepdims=True)
log_X = np.log(X)
X_clr = log_X - log_X.mean(axis=1, keepdims=True)

clr_df = pd.DataFrame(X_clr, columns=components, index=df_wide.index)

print(f"CLR transform   : row sums ≈ 0 (max abs = {np.abs(X_clr.sum(axis=1)).max():.2e})  ✓")
# ── 4. Shapiro-Wilk per variable ──────────────────────────────────────────────
print("\n" + "=" * 65)
print("STEP 1: Shapiro-Wilk Normality Test (CLR values, per variable)")
print("-" * 65)
print(f"  {'Variable':<36}  {'W':>7}  {'p-value':>8}  {'Normal?':>8}")
print(f"  {'-'*36}  {'-'*7}  {'-'*8}  {'-'*8}")

any_nonnormal = False
for comp in components:
    vals = clr_df[comp].values
    if len(vals) >= 3:
        W, p = shapiro(vals)
        normal = "yes" if p >= 0.05 else "NO *"
        if p < 0.05:
            any_nonnormal = True
        print(f"  {comp:<36}  {W:>7.4f}  {p:>8.4f}  {normal:>8}")
    else:
        print(f"  {comp:<36}  {'n<3':>7}  {'n<3':>8}  {'n/a':>8}")

print(f"\n  * p < 0.05 = non-normal distribution of CLR values.")
print(f"  → {'Non-normality detected: Kruskal-Wallis is appropriate.' if any_nonnormal else 'All normal: Kruskal-Wallis still valid as robust non-parametric test.'}")

# ── 5. Kruskal-Wallis across all variables ────────────────────────────────────
#
# Here each VARIABLE is treated as a group.
# The CLR-transformed values across all 377 samples are compared between
# the 9 variables to test: do these categories have different distributions?

print("\n" + "=" * 65)
print("STEP 2: Kruskal-Wallis Test — all variables compared")
print("-" * 65)

group_vals = [clr_df[comp].values for comp in components]
H, p_kw = kruskal(*group_vals)
sig = "***" if p_kw < 0.001 else ("**" if p_kw < 0.01 else ("*" if p_kw < 0.05 else "ns"))

print(f"\n  H-statistic : {H:.4f}")
print(f"  p-value     : {p_kw:.4e}")
print(f"  df          : {len(components) - 1}")
print(f"  Result      : {sig}  —  {'Variables differ significantly' if p_kw < 0.05 else 'No significant difference detected'}")

# ── 6. Dunn post-hoc with Holm correction ─────────────────────────────────────
print("\n" + "=" * 65)
print("STEP 3: Dunn Post-hoc Test — pairwise variable comparisons (Holm)")
print("=" * 65)

if p_kw >= 0.05:
    print("\n  KW not significant — post-hoc not warranted.")
else:
    # Build long format: one row per (sample, variable) with CLR value
    long_clr = clr_df.reset_index().melt(
        id_vars="sample_id", var_name="variable", value_name="clr_value"
    )

    dunn = sp.posthoc_dunn(
        long_clr, val_col="clr_value", group_col="variable", p_adjust="holm"
    )

    # Reorder to match original component order
    dunn = dunn.loc[components, components]

    # ── Print full matrix ──
    col_w = 14
    short = [c[:col_w] for c in components]

    print(f"\n  Holm-adjusted p-values (* p<0.05  ** p<0.01  *** p<0.001  ns)\n")
    # header
    print(f"  {'':36}", end="")
    for s in short:
        print(f"  {s:>{col_w}}", end="")
    print()
    print("  " + "-" * (36 + (col_w + 2) * len(components)))

    for row_c in components:
        print(f"  {row_c[:36]:<36}", end="")
        for col_c in components:
            if row_c == col_c:
                cell = "—"
            else:
                p_val = dunn.loc[row_c, col_c]
                stars = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))
                cell = f"{p_val:.3f}{stars}"
            print(f"  {cell:>{col_w}}", end="")
        print()

    # ── Summary of significant pairs ──
    print(f"\n  Significant pairs (Holm p < 0.05):")
    found = False
    for i, r in enumerate(components):
        for j, c in enumerate(components):
            if j <= i:
                continue
            p_val = dunn.loc[r, c]
            if p_val < 0.05:
                stars = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else "*")
                print(f"    {r}  vs  {c}  →  p = {p_val:.4f} {stars}")
                found = True
    if not found:
        print("    None after Holm correction.")
    # ── Summary of NON-significant pairs ──
    print(f"\n  Non-significant pairs (Holm p >= 0.05):")
    found = False
    for i, r in enumerate(components):
        for j, c in enumerate(components):
            if j <= i:
                continue
            p_val = dunn.loc[r, c]
            if p_val >= 0.05:
                print(f"    {r}  vs  {c}  →  p = {p_val:.4f}")
                found = True
    if not found:
        print("    None after Holm correction.")
    # ── Summary of NON-significant pairs and replace component names using dict_rename──
    print(f"\n  Non-significant pairs (Holm p >= 0.05) with renamed components:\n")
    found = False
    for i, r in enumerate(components):
        for j, c in enumerate(components):
            if j <= i:
                continue
            p_val = dunn.loc[r, c]
            if p_val >= 0.05:
                r_renamed = dict_rename.get(r, r)
                c_renamed = dict_rename.get(c, c)
                print(f"    {r_renamed}  vs  {c_renamed}  →  p = {p_val:.4f}")
                found = True
    if not found:
        print("    None after Holm correction.")

print("\n" + "=" * 65)
print("INTERPRETATION")
print("=" * 65)
print("""
What is being tested:
  Each of the 9 composition categories (variables) is treated as a
  group. The CLR-transformed values across all samples are compared
  to ask: do these categories occupy different positions in the
  composition, i.e. do they differ in their relative abundance?

  Kruskal-Wallis: omnibus test — do ANY variables differ?
  Dunn + Holm:    which specific pairs of variables differ?

Why CLR first:
  Raw percentages sum to 100% and are not independent. CLR(xᵢ) =
  log(xᵢ / geometric_mean(x)) removes this constraint so comparisons
  between variables are valid.
""")

