# CRISPRlysis

**Annotation of CRISPR–Cas repair outcome — estimating cut type (blunt vs staggered), deletion orientation, and the DNA end participating in repair.**

This code is part of the publication Přibylová et al., 2026 - 5′ Complementarity-Mediated End Joining (5′CMEJ) DNA repair; author's manuscript [doi.org/10.64898/2026.07.23.740330](https://doi.org/10.64898/2026.07.23.740330).

**CRISPRlysis** takes the allele frequency tables produced by CRISPResso2* from amplicon Next Generation Sequencing (NGS) data and turns them into a fully annotated table of repair events, in which every read is classified by the position, orientation and complementarity context of its deletion relative to the predicted cleavage site. (Any table which has the same column system as CRISPResso2 can be used.) It is optimised for *Sp*Cas9, but also should work for Cas12a and any other endonuclease, just the parameters must be adjusted. It then produces a standard set of publication-ready charts and accompanying descriptive statistics for each analysis.

*Clement, K., H. Rees, M.C. Canver, J.M. Gehrke, R. Farouni, J.Y. Hsu, M.A. Cole, D.R. Liu, J.K. Joung, D.E. Bauer, and L. Pinello. 2019. CRISPResso2 provides accurate and rapid genome editing sequence analysis. *Nat. Biotechnol.* 37:224–226. doi:[10.1038/s41587-019-0032-3](https://doi.org/10.1038/s41587-019-0032-3).

---

## What the script analyses

In brief, CRISPRlysis performs a complex analysis of mutation events, covering:

- **Overall event composition** — frequencies and ratios of insertion, deletion and substitution events around the cut site.
- **Insertions** — the number of inserted nucleotides and, in the case of single-nucleotide insertions, whether they are templated by the 4th nucleotide upstream of the PAM.
- **Deletion orientation** — the ratio of deletions that lie upstream or downstream of the predicted cut site, or that span it.
- **Cut geometry** — the ratio of possible blunt and staggered cut formation, how each contributes to deletions upstream and downstream of the cut, and the effect of the 4th nucleotide upstream of the PAM sequence.
- **End participation** — which 5′/3′ DNA end is likely to participate in the repair.

Examples of every output are shown [below](#graphical-outputs) and in the
[example analysis](CRISPRlysis_outputs/Pribylova_et_al_2022/WT) included in this repository.

## Scope

CRISPRlysis can analyse any CRISPResso2 (version 2.2.11) `Alleles_frequency_table_*.txt` output file, or any file with the same structure. It is not restricted to a single nuclease: any type of cut (*Sp*Cas9, *Lb*Cas12a, …) can be analysed, and only the nuclease and cut positions need to be modified in `parameters.txt`. The default graphical outputs were optimised for *Sp*Cas9, but can be
adapted to other endonucleases by changing the script parameters — see the [Parameters](#parameters) section below.

---

## Quick start

```bash
git clone https://github.com/pribylad/CRISPRlysis.git
cd CRISPRlysis
python CRISPRlysis.py
```

`CRISPRlysis.py` reads **`parameters.txt`** from the same folder. That file is a commented Python dictionary; edit it and re-run with your datasets. **Run the script from inside the repository folder**: `parameters.txt` is opened as `./parameters.txt`, i.e. relative to the current working directory. No command-line arguments are used.

The repository ships with a working example: `parameters.txt` as committed points at `./testing_dataset` and reproduces the outputs under `CRISPRlysis_outputs/Pribylova_et_al_2022/WT/` (the date in the file and folder names is the date of the run).

### Repository layout

```
CRISPRlysis.py            entry point — reads parameters.txt, builds tables, draws all charts
parameters.txt            the only file you normally edit
scripts/
  make_mastertables.py    alignment parsing, gap resolution, MH detection, cut-type logic
  plot_functions.py       violin, stacked horizontal/vertical barchart and stats helpers
  dict_functions.py       colour and hatch palettes for cut types and orientations
  general_functions.py    reporting, path handling, file discovery
testing_dataset/          32 CRISPResso2 outputs (8 sgRNAs × 4 replicates), SpCas9, N. benthamiana
CRISPRlysis_outputs/      worked example produced from testing_dataset
```

---

## Parameters

The dictionary in `parameters.txt` is documented inline; the entries that matter most:

| Key                                            | Meaning                                                                                                     |
| ---------------------------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `source_dir`                                 | Folder searched recursively for`Alleles_frequency_table_around*.txt`                                      |
| `project_output_dir`                         | Where tables, charts and`report.txt` are written                                                          |
| `cut_position`                               | Optional. The cut is detected for every sample from the gRNA and the PAM; this value (the`offset_around_cut_to_plot` of the CRISPResso2 run) is used only as a fallback for nucleases without a PAM rule |
| `cut_maximum`                                | Most 5′ cleavage position relative to the optimum (e.g.`-2` for *Sp*Cas9)                              |
| `cut_minimum`                                | Most 3′ cleavage position relative to the optimum (e.g.`0` for *Sp*Cas9)                               |
| `nuclease`                                   | `SpCas9` or `LbCas12a`or yours                                                                          |
| `minimum_Reads`                              | Read-count threshold applied when building`01_mastertable` (typically 5)                                  |
| `data_source_original`                       | Provenance label written into every row and onto every chart                                                |
| `add_columns`                                | Constant metadata columns — organism, delivery, cell-cycle status, cell line, hpi, ...                     |
| `filter_out_values`                          | `{column: [values to drop]}`; `{}` disables filtering                                                   |
| `subtables_param`                            | Columns to split on - separate analysis will be done on each unique combination, e.g.`['gene', 'mutant']` |
| `plot_INDELs`, `plot_INs`, `plot_DELs`, `plot_DEL_individual_gRNAs`, `plot_DEL_dinochart`, `plot_participating_end`, `plot_participating_end_dinochart` | Switches (`True`/`False`) for the chart families described below; all are required in `parameters.txt`. `plot_MH` and `plot_MH_dinochart` may still be present in older files but are **not read** by the code |
| `show_plots`                                 | `True` shows each chart on screen while the script runs (keep `False` for batch runs)                       |
| `sort_alphabetically`                        | `True` sorts the samples of `01_mastertable` by name before plotting                                        |
| `IN_greater_than`                            | Insertion sizes above this are pooled into a single category (max 7)                                        |
| `output_format`                              | `png`, `svg`, or any matplotlib-supported format                                                        |
| `make_mastertable_00/01`, `make_subtables` | Set to `False` to reuse an existing table via the matching `path_to_existing_mastertable_00`, `path_to_existing_mastertable_01` or `path_to_existing_subtable_02` key (keep the key, `''` if unused) |

### Optional parameters (added 2026-11)

Every key below is optional. **With none of them set, the tables are those of the previous release** (plus a few extra columns written by the automatic cut detection, e.g. `Cut_position`, `Cut_position_in_alignment`, and a cut-detection line per sample in `report.txt`). The audited analyses of Pallaseni et al. 2024 and of the bioRxiv meta-analysis used the values in the right-hand column.

| Key | Default | Used in the audited runs | Meaning |
| --- | --- | --- | --- |
| `use_per_row_cut` | `False` | `True` | Analyse every allele at its own cut column. Needed when a read carries an insertion, because CRISPResso2's fixed-width window then shifts the reference. Changes results for those rows. |
| `insertion_side` | `'left'` | `'left'` | On which side of an aligned insertion the cut column is counted: `'left'` (the inserted bases read as downstream of the cut) or `'right'`. |
| `substitution_window` | `3` | `'allowed_cut'` | Which columns count as "at the cut" for substitutions: an integer *n* (cut-*n*+1 … cut+*n*), or `'allowed_cut'` (the nucleotides whose deletion counts as at the cut, 4 columns for −2,0). |
| `revert_substitutions_outside_window` | `False` | `'all'` | Substitutions outside the window are treated as sequencing noise and reverted to the reference, so a read whose only difference is a distant substitution counts as unedited. `'all'` applies it to every allele, `'pure_sub'` only to alleles without an indel. |
| `protect_substitutions_adjacent_to_indel` | `1` | `3` | Substitutions within this many columns of a gap run are kept (they may be part of the indel). |
| `skip_multi_gap_run_alleles` | `True` | `True` | Alleles with two or more separate gap runs (an aligner-split event) are not corrected. |
| `collapse_identical_alleles` | `False` | `True` | After the correction, alleles that became identical are merged (reads summed) before the `minimum_Reads` filter. |
| `plot_INDEL_horizontal_barcharts`, `plot_IN_…`, `plot_DEL_…`, `plot_participating_end_horizontal_barcharts` | `True` | `False` | One bar per sample in files of 40; switch off for pooled panels with thousands of targets. |
| `random_gen_data` | `False` | `False` | Only for the randomised validation datasets (one row = one unique random target): the original sequence is taken from `Reference_Sequence` of each row. Leave unset for CRISPResso2 data. |
| `list_of_columns` | not set | Cucuy 2026, Weiss 2022 | Column names for folder names that do not follow the naming convention below. **It must be carried over every time such a dataset is re-run**: without it the folder names are split with the default scheme, the columns that `subtables_param` refers to do not exist and the run stops with a `KeyError` (e.g. `'generation'`). |

The cut detection (`scripts/cut_site.py`) tolerates reference variants that overhang the amplicon (CRISPResso2 can shift the quantification window by a few bp between alleles): the amplicon is assembled from overlapping variants (at least 20 bp overlap) before the consistency check.

For *Sp*Cas9 the code checks that positions 2–3 of the inferred PAM are `GG` and raises a `ValueError` if not — this catches a misspecified `cut_position` before it silently propagates into every downstream column. We recommend to check 03_07_*.png files if the sequence is correct.

### Input naming convention

Two names are read from the CRISPResso2 output:

- the **folder name** of each run, which must start with a 14-character prefix (`CRISPResso_on_`, as written by CRISPResso2); the code drops these first 14 characters and splits the rest at `_`;
- the **file name** `Alleles_frequency_table_around_sgRNA_<gRNA sequence>.txt` (the gRNA sequence is the last element; CRISPResso2 writes this name when no gRNA name is given).

The sample name is built as `sgRNA_` + folder name without the prefix. The default scheme expects eight `_`-separated fields after the prefix:

| Position | Field in the folder name | Column in the tables | Example |
| --- | --- | --- | --- |
| 1 | Order | `order` | `1` |
| 2 | Nuclease (Cas variant) | `Cas_variant` | `SpCas9` |
| 3 | Targeted gene | `gene` | `P35S` |
| 4 | gRNA / target number | `gRNA_num` | `1` |
| 5 | Mutant (genotype) | `mutant` | `WT` |
| 6 | Number of targets | `no_of_targets` | `1` |
| 7 | Overhangs | `overhang` | `0` |
| 8 | Replicate | `replicate` | `1` |

Example: the folder `CRISPResso_on_1_SpCas9_P35S_1_WT_1_0_1` gives the sample `sgRNA_1_SpCas9_P35S_1_WT_1_0_1`. These column names are the ones to use in `subtables_param` (e.g. `['gene', 'mutant']`) and in `filter_out_values`. Any other scheme can be supplied through the optional `list_of_columns` parameter (the names are assigned to the fields after the prefix in the given order).

---

## Tabular outputs

| File                                   | Content                                                                                     |
| -------------------------------------- | ------------------------------------------------------------------------------------------- |
| `00_mastertable_<date>.tsv`          | Every allele from every input file, annotated — 62 columns (64 with the optional substitution correction) |
| `01_mastertable_<date>.tsv`          | Read-depth filtered, sample metadata split into columns (the fields of the folder name and the `add_columns` entries), sorted — 76 columns (78 with the optional substitution correction) |
| `02_mastertable_<subset>_<date>.tsv` | One table per unique combination of`subtables_param`                                      |
| `report.txt`                         | Full run log: files read, filters applied, table shapes, every chart written, total runtime |

Per-amplicon annotated tables (`00_mastertable_<sample>.tsv`) are written to the subfolder `00_per_amplicon/` of the output folder (nothing is written next to the CRISPResso2 input), so any anomaly in the mastertable can be traced back to a single input file.

### Key annotation columns

| Column                                                                           | Meaning                                                                                                                                     |
| -------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `Gap_min_max_ranges`                                                           | Leftmost and rightmost equivalent placement of the deletion                                                                                 |
| `Gap_around_cut`                                                               | Whether the deletion falls within the permissible cleavage window                                                                           |
| `n_deleted_around_cut`                                                         | Size of the deletion around the cut (independent of CRISPResso2's`n_deleted`)                                                             |
| `MH_to_proximal_len`, `MH_to_distal_len`, `Maximal_MH_len`, `Maximal_MH` | Complementarity on each side of the junction, its maximum, and its sequence                                                                 |
| `Gap_category_1`, `Gap_category_2`                                           | Deletion category (1–12) for each extreme placement; `97` = no deletion (the coordinates are not numbers), `98` = the range of the deletion could not be determined |
| `Cut_type`                                                                     | `0` blunt, negative = staggered upstream, positive = staggered downstream,`C` undetermined, `O` outside the window, `X` no deletion |
| `DEL_orientation`                                                              | `U` upstream, `D` downstream, `B` bidirectional, `C` undetermined, `X` no deletion                                                  |
| `End_participating_on_MH`                                                      | Which end carries the complementarity:`U`, `D`, `C`, `X`                                                                            |
| `End_3_5`                                                                      | DNA end participating in repair:`3`, `5`, or `0` (either or cannot be determined)                                                     |
| `Likely_cut_type`                                                              | `0` if the geometry is unambiguous, `1` if inferred from complementarity length                                                         |
| `Event`, `Event_around_cut`                                                  | `DEL`/`IN`/`SUB`/`INDEL`/`DELSUB`/`INSUB`/… or `no_mutation`, before and after discarding substitutions distant from the cut        |

The `Likely_cut_type` flag is deliberate: it lets any downstream analysis be repeated on the confidently assigned subset alone as a sensitivity check.

---

## Graphical outputs

The charts shown below are the output of the shipped example (code of commit `6f6edce`, `parameters.txt` as committed, run of 2026-11-21; the date in the chart titles is the date of the run). They are stored in `image/example_outputs/`.

Charts are written to a dated subfolder named after the subtable, e.g. `CRISPRlysis_outputs/Pribylova_et_al_2022/WT/P35S_WT_261121/`. **Every chart is accompanied by a `<name>_stat.txt` file** containing count, mean, standard deviation, minimum, quartiles, median and maximum for each plotted category, so the numbers behind a figure are always available without re-running the analysis. Charts carry the data source, the number of gRNAs and replicates, and a footnote stating exactly what was filtered and how it was normalised — the figure is self-documenting.

Next to each chart (`.png`) the plotted numbers are saved as `<name>_data.tsv` (the violin chart `01_01` also as `<name>_<date>.tsv`). Files are numbered `NN_MM_…` so that they sort into the four analysis families below. Where more
than 40 samples are present, per-sample charts are automatically split across several files
(`_1`, `_2`, …).

### 01 — Overall INDEL/SUB composition

| File                                  | Chart                                                                                                                              |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| `01_01_INDELSUB_summary_violinplot` | Violin plot of IN / DEL / SUB / OTHERS as a percentage of reads, one point per sample overlaid. Exported alongside its own`.tsv` |
| `01_02_i_…_not_normalised`         | Stacked horizontal barchart, one bar per sample, absolute % of reads                                                               |
| `01_03_i_…_normalised`             | The same, rescaled so IN + DEL + SUB + OTHERS = 100%                                                                               |

`OTHERS` pools the compound events (INDEL, INSUB, DELSUB, INDELSUB). All categories are classified only for events falling around the cut.

![INDEL/SUB summary violin plot](image/example_outputs/01_01_INDELSUB_summary_violinplot.png)

![INDELSUB horizontal barchart, not normalised](image/example_outputs/01_02_1_INDELSUB_horizontal_barchart_all_samples_not_normalised.png)

![INDELSUB horizontal barchart, normalised](image/example_outputs/01_03_1_INDELSUB_horizontal_barchart_all_samples_normalised.png)

### 02 — Insertions

| File                                                    | Chart                                                                                                          |
| ------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `02_01_IN_no_of_IN_nucleotides_summary_violinplot`    | Violin plot of the normalised distribution of insertion sizes (1 nt, 2 nt, … pooled above`IN_greater_than`) |
| `02_02_IN_no_of_IN_nucleotides_stacked_barchart_name` | Per-sample stacked barchart of the same size distribution                                                      |
| `02_03_i_4th_nt_inserted_not_normalised`              | Is the inserted base a copy of the 4th nucleotide upstream of the PAM? Per-sample stacked barchart, True/False |
| `02_04_i_4th_nt_inserted_normalised`                  | The same as 02_03*, normalised to True + False = 100%                                                          |
| `02_05_IN_4th_nt_inserted_normalised_violinplot`      | Violin plot of the normalised True/False split, points coloured by the identity of the 4th nucleotide          |
| `02_06_IN_4th_nt_inserted_normalised_violinplot`      | The same with medians drawn separately for each 4th nucleotide                                                 |

The 02_03–02_06 series is restricted to 1 nt insertions around the cut in sequences with no deletion or substitution around the cut, which is what makes the templated-insertion signal interpretable.

![Insertion size, violin plot](image/example_outputs/02_01_IN_no_of_IN_nucleotides_summary_violinplot.png)

![Insertion size, stacked barchart per sample](image/example_outputs/02_02_IN_no_of_IN_nucleotides_stacked_barchart_name.png)

![4th nucleotide inserted, not normalised](image/example_outputs/02_03_1_4th_nt_inserted_not_normalised.png)

![4th nucleotide inserted, normalised](image/example_outputs/02_04_1_4th_nt_inserted_normalised.png)

![4th nucleotide inserted, violin plot](image/example_outputs/02_05_IN_4th_nt_inserted_normalised_violinplot.png)

![4th nucleotide inserted, violin plot with medians per nucleotide](image/example_outputs/02_06_IN_4th_nt_inserted_normalised_violinplot.png)

### 03 — Deletions, cut type and orientation

| File                                                              | Chart                                                                                                                                                                                                                                                                             |
| ----------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `03_01` / `03_02_DEL_orientation_summary_violinplot`          | Violin plots of deletion orientation (upstream / bidirectional / downstream / cannot be determined), not normalised and normalised                                                                                                                                                |
| `03_03_i_Stagg_blunt_ratio_all_samples`                         | Per-sample stacked horizontal barchart of the full cut-type spectrum: blunt, −1 nt and −2 nt staggered (for*Sp*Cas9), each with and without complementarity, plus undetermined and window-spanning deletions                                                                  |
| `03_04_DEL_Stagg_blunt_summary_violinplot`                      | Violin plot of the same categories across samples as in 03_03*                                                                                                                                                                                                                    |
| `03_05` / `03_06_DEL_Stagg_blunt_oriented_summary_violinplot` | The cut-type spectrum split by orientation relative to the cut — upstream blunt, downstream −1 nt staggered with complementarity, and so on — not normalised and normalised. This is the chart that separates upstream- from downstream-biased repair (optimised fo*Sp*Cas9) |
| `03_07_nnn_…_dino_30nt_NotNormalised_<sample>`                 | Per-sample deletion-spectrum chart (see below)                                                                                                                                                                                                                                    |
| `03_08_DEL_Stagg_blunt_summary_barchart_4th_nt`                 | Cut-type composition grouped by the 4th nucleotide upstream of the PAM                                                                                                                                                                                                            |
| `03_09` / `03_10_DEL_Stagg_blunt_summary_violinplot_4th_nt`   | The same relationship as in 03_08* as violin plots, keyed by nucleotide and by cut type respectively                                                                                                                                                                              |

![Deletion orientation, not normalised](image/example_outputs/03_01_DEL_orientation_summary_violinplot.png)

![Deletion orientation, normalised](image/example_outputs/03_02_DEL_orientation_summary_violinplot.png)

![Cut-type spectrum per sample](image/example_outputs/03_03_1_Stagg_blunt_ratio_all_samples.png)

![Cut-type spectrum, violin plot](image/example_outputs/03_04_DEL_Stagg_blunt_summary_violinplot.png)

![Cut-type spectrum by orientation, not normalised](image/example_outputs/03_05_DEL_Stagg_blunt_oriented_summary_violinplot.png)

![Cut-type spectrum by orientation, normalised](image/example_outputs/03_06_DEL_Stagg_blunt_oriented_summary_violinplot_normalised.png)

![Cut type by 4th nucleotide, barchart](image/example_outputs/03_08_DEL_Stagg_blunt_summary_barchart_4th_nt.png)

![Cut type by 4th nucleotide, violin plot (nucleotide)](image/example_outputs/03_09_DEL_Stagg_blunt_summary_violinplot_4th_nt.png)

![Cut type by 4th nucleotide, violin plot (cut type)](image/example_outputs/03_10_DEL_Stagg_blunt_summary_violinplot_4th_nt.png)

### 04 — DNA end participating in repair

| File                                                             | Chart                                                                                                        |
| ---------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| `04_01_i` / `04_02_i_DEL_end_summary_…_horizontal_barchart` | Per-sample stacked barchart of the participating end — 3′, 5′, or either — normalised and not normalised |
| `04_03` / `04_04_DEL_end_summary_…_violinplot`              | The same distributions as violin plots across samples                                                        |
| `04_05_nnn_end_participating_…_<sample>`                      | Per-sample deletion-spectrum chart coloured by participating end                                             |

Restricted to deletions that carry complementarity, since the participating end is undefined without it.

![Participating end, normalised barchart](image/example_outputs/04_01_1_DEL_end_summary_norm_horizontal_barchart.png)

![Participating end, not normalised barchart](image/example_outputs/04_02_1_DEL_end_summary_NOTnorm_horizontal_barchart.png)

![Participating end, normalised violin plot](image/example_outputs/04_03_DEL_end_summary_norm_violinplot.png)

![Participating end, not normalised violin plot](image/example_outputs/04_04_DEL_end_summary_NOTnorm_violinplot.png)

### The per-sample deletion-spectrum chart

The `03_07` and `04_05` series are the most information-dense output. Each is a vertical stacked barchart in which:

- the **x axis** is the number of deleted nucleotides, signed by direction — negative to the upstream of the probable cut, positive to the downstream, spanning 30 nt either side;
- the **y axis** is the summed percentage of reads;
- a **vertical line** marks the cut site, and the **target sequence around the cut** is printed in the title with `|` at the '0' (blunt) cleavage position;
- **stacking and hatching** encode cut type, orientation and the presence of complementarity (`03_07`) or the participating 3′/5′ end (`04_05`).

Bidirectional deletions are excluded and orientation-undetermined deletions are split between upstream and downstream, so the upstream/downstream asymmetry of the plot reads directly as the asymmetry of end processing at that target. One file is produced per sample, numbered sequentially.

![Deletion-spectrum chart of one sample (03_07)](image/example_outputs/03_07_example_sample_1.png)

![Deletion-spectrum chart of one sample coloured by the participating end (04_05)](image/example_outputs/04_05_example_sample_1.png)

---

## Worked example

`CRISPRlysis_outputs/Pribylova_et_al_2022/WT/` contains a complete run over the shipped `testing_dataset`: 8 sgRNAs × 4 replicates targeting *P35S* in *Nicotiana benthamiana*, *Sp*Cas9 delivered by agroinfiltration, sampled 120 h post-inoculation in non-dividing tissue, with `cut_position = 120`, `cut_maximum = -2`, `cut_minimum = 0`. It is the fastest way to see the shape of every output before running your own data.

*Note on the stored example:* the example was produced on 2026-11-21 with the code of the commit `6f6edce` and `parameters.txt` as committed (default settings of the optional parameters). The analyses of the manuscript used the optional parameters listed above (right-hand column of the table); with them the absolute percentages of reads are about 10 % higher than in the example (e.g. mean insertions 12.1 % → 13.2 % of reads, deletions 8.8 % → 9.9 %), the shape of the charts is the same.

---

## Requirements

Python 3 with `pandas`, `numpy`, `matplotlib`, `seaborn`, and `more-itertools`. Tested with:

| Package        | Version |
| -------------- | ------- |
| Python         | 3.11.5  |
| pandas         | 2.1.4   |
| matplotlib     | 3.8.0   |
| seaborn        | 0.13.2  |
| numpy          | 1.26.4  |
| more-itertools | 10.1.0  |

---

## Licence

Creative Commons Attribution 4.0 International (CC-BY-4.0). See `LICENSE.txt`.

---

## Citation

If you use CRISPRlysis, please cite the accompanying manuscript.
*[doi.org/10.64898/2026.07.23.740330](https://doi.org/10.64898/2026.07.23.740330).*

---

## Usage of AI

This code was developed since 2021, without AI. Since late 2025 AI was used to polish the code, generates descriptions of functions and those minor tasks. All logic behind DEL classification etc. is human based.

---

## Author

Adéla Přibylová, Faculty of Science, Charles University, Prague, Czech Republic.
