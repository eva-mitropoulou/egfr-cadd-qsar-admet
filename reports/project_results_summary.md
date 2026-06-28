# EGFR QSAR/CADD Project Results Summary

## Dataset

- Model-ready molecules: 10,593
- Activity endpoint: ChEMBL EGFR IC50 converted to median pIC50 per molecule.
- Broad sanity filters removed extreme activity/property artifacts while retaining most molecules.

## Feature-Set Comparison

| feature_set | model | MAE | RMSE | R2 |
| --- | --- | --- | --- | --- |
| RDKit descriptors | Random forest | 0.745 | 0.979 | 0.468 |
| Morgan fingerprints | Random forest | 0.516 | 0.712 | 0.719 |
| Fingerprints + descriptors | Random forest | 0.526 | 0.719 | 0.713 |

Morgan fingerprints captured EGFR-relevant substructure information better than broad descriptors alone.

## Random vs Scaffold Split

| validation | feature_set | MAE | RMSE | R2 |
| --- | --- | --- | --- | --- |
| Random split | Morgan fingerprints | 0.516 | 0.712 | 0.719 |
| Scaffold split | Morgan fingerprints | 0.667 | 0.871 | 0.550 |

Random split was optimistic; scaffold split gave a harder estimate of generalization to new chemotypes.

## Cross-Validation

| validation_scheme | MAE | RMSE | R2 |
| --- | --- | --- | --- |
| Random KFold | 0.511 +/- 0.009 | 0.703 +/- 0.008 | 0.724 +/- 0.006 |
| Scaffold GroupKFold | 0.618 +/- 0.037 | 0.823 +/- 0.042 | 0.616 +/- 0.034 |

Model performance was stable under random CV but dropped under scaffold-aware CV.

## Applicability Domain

| similarity_bin | count | MAE | RMSE |
| --- | --- | --- | --- |
| <0.3 | 149 | 0.957 | 1.199 |
| 0.3-0.5 | 792 | 0.842 | 1.072 |
| 0.5-0.7 | 3372 | 0.745 | 0.947 |
| >0.7 | 6280 | 0.514 | 0.697 |

Prediction error decreased as maximum Tanimoto similarity to the training chemistry increased.

## Candidate Ranking

- Ranked molecules: 10,593
- Top 20 predicted pIC50 range: 8.795-10.029
- Diverse top 20 predicted pIC50 range: 8.689-10.029
- Diverse top 20 unique scaffolds: 20
- Diverse top 20 model risk: 17 low-risk, 3 medium-risk
- Lipinski-clean diverse top 20: 19/20

Top diverse candidates are saved in `reports/top_20_diverse_candidates.md`.

## Limitations

- This is a retrospective ChEMBL portfolio project, not prospective drug discovery.
- IC50 values come from heterogeneous assays and publications.
- The ADMET layer is proxy/drug-likeness triage, not true ADMET prediction.
- Docking was not included.
- Predictions are less reliable outside the model applicability domain.

## Generated Figures

- `figures/model_performance_comparison.png`
- `figures/random_vs_scaffold_validation.png`
- `figures/applicability_domain_bins.png`
