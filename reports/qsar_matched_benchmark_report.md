# QSAR Matched Benchmark Report

This report assembles the completed EGFR baseline outputs into a matched benchmark summary and trains a lightweight full-data primary model for downstream CLI/demo use.

## Split Benchmarks

| split | model | feature_set | train_size | test_size | MAE | RMSE | R2 | Pearson | Spearman |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| random_split | Random Forest | rdkit_descriptors | 8474 | 2119 | 0.745 | 0.979 | 0.468 | None | None |
| random_split | Random Forest | morgan_fingerprints | 8474 | 2119 | 0.516 | 0.712 | 0.719 | None | None |
| random_split | Random Forest | combined | 8474 | 2119 | 0.526 | 0.719 | 0.713 | None | None |
| scaffold_split | Random Forest | morgan_fingerprints | 8471 | 2122 | 0.667 | 0.871 | 0.550 | None | None |
| random_split | Mean baseline | none | 8474 | 2119 | 1.104 | 1.342 | -0.001 | None | None |

## Scaffold Split Integrity

- Train molecules: 8,471
- Test molecules: 2,122
- Train scaffolds: 3,117
- Test scaffolds: 570
- Scaffold overlap: 0
- Standardized molecule overlap: 0

## Cross-Validation Summary

| validation_scheme | MAE | RMSE | R2 |
| --- | --- | --- | --- |
| random_kfold | 0.511 | 0.703 | 0.724 |
| scaffold_groupkfold | 0.618 | 0.823 | 0.616 |

## Primary Model Selection

- Best model by scaffold performance: Random Forest
- Scaffold MAE: 0.667
- Scaffold RMSE: 0.871
- Scaffold R2: 0.550
- Primary full-data model saved for CLI/demo: Ridge on Morgan fingerprints
