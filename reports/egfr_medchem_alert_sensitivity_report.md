# EGFR Medchem-Alert Sensitivity Report

PAINS, Brenk, and external unwanted-substructure SMARTS alerts were used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.

This sensitivity run uses Morgan fingerprints and the same split logic/random state as the primary Random Forest baseline, with 30 trees for runtime. The official 100-tree primary benchmark files are preserved unchanged.

## Subset Composition

| subset | row_count | fraction_retained | molecule_count | scaffold_count | pIC50_mean | pIC50_median | pIC50_std | pIC50_mean_shift_vs_full |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| main_all_model_ready | 10593 | 1.000 | 10593 | 3685 | 6.925 | 7.026 | 1.338 | 0.000 |
| pains_excluded | 9746 | 0.920 | 9746 | 3369 | 6.971 | 7.071 | 1.320 | 0.046 |
| pains_brenk_excluded | 4248 | 0.401 | 4248 | 1572 | 6.871 | 6.921 | 1.372 | -0.054 |
| strict_medchem_clean | 3669 | 0.346 | 3669 | 1401 | 6.677 | 6.680 | 1.349 | -0.248 |

## Morgan RF Sensitivity Metrics

| subset | split | train_size | test_size | MAE | RMSE | R2 | Pearson | Spearman |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| main_all_model_ready | random_split | 8474 | 2119 | 0.522 | 0.719 | 0.713 | 0.845 | 0.849 |
| main_all_model_ready | scaffold_split | 8447 | 2146 | 0.688 | 0.891 | 0.595 | 0.782 | 0.770 |
| pains_excluded | random_split | 7796 | 1950 | 0.527 | 0.715 | 0.712 | 0.844 | 0.840 |
| pains_excluded | scaffold_split | 7780 | 1966 | 0.687 | 0.887 | 0.574 | 0.770 | 0.760 |
| pains_brenk_excluded | random_split | 3398 | 850 | 0.533 | 0.726 | 0.725 | 0.852 | 0.856 |
| pains_brenk_excluded | scaffold_split | 3389 | 859 | 0.781 | 0.999 | 0.596 | 0.783 | 0.793 |
| strict_medchem_clean | random_split | 2935 | 734 | 0.581 | 0.775 | 0.696 | 0.837 | 0.841 |
| strict_medchem_clean | scaffold_split | 2914 | 755 | 0.817 | 1.009 | 0.579 | 0.790 | 0.794 |

## Top-Ranked Alert Composition

- Top-20 clean count: 20
- Top-20 alert-flagged count: 0
- Diverse top-20 clean count: 20
- Diverse top-20 alert-flagged count: 0

## Interpretation

The strict medchem-clean subset changed scaffold-split RMSE by 0.118 pIC50 units relative to the full sensitivity baseline.
This is a sensitivity analysis, not a replacement for the primary full-dataset benchmark.
