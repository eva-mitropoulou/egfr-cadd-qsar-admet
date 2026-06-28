# Applicability Domain Report

Applicability domain was estimated from max Tanimoto similarity to the training set during scaffold GroupKFold validation.

| similarity_bin | count | MAE | RMSE | mean_similarity |
| --- | --- | --- | --- | --- |
| low | 149 | 0.957 | 1.199 | 0.256 |
| medium | 4185 | 0.762 | 0.971 | 0.582 |
| high | 6259 | 0.513 | 0.697 | 0.804 |

- Low-similarity MAE: 0.957
- High-similarity MAE: 0.513
- Out-of-domain molecules: 149

Predictions with max Tanimoto < 0.3 are flagged as out-of-domain.
