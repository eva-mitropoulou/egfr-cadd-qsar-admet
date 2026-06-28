# Project Completion Audit

This audit lists available artifacts and missing critical outputs without exposing raw molecule records.

## Artifact Counts

- raw_data: 1
- processed_data: 10
- results: 12
- reports: 77
- figures: 20
- report_figures: 16
- models: 3
- notebooks: 4
- scripts: 2
- src: 37

## Missing Critical Artifacts

- None

## CSV Shape Metadata

- `data/raw/egfr_chembl_ic50_raw.csv`: 26600 rows, 9 columns
- `data/processed/egfr_descriptors.csv`: 10834 rows, 14 columns
- `data/processed/egfr_ic50_clean.csv`: 10834 rows, 5 columns
- `data/processed/egfr_model_ready.csv`: 10593 rows, 14 columns
- `data/processed/egfr_standardized_molecules.csv`: 10593 rows, 19 columns
- `data/processed/egfr_structure_candidates.csv`: 3 rows, 4 columns
- `data/processed/features_combined_index.csv`: 10593 rows, 5 columns
- `data/processed/features_morgan_index.csv`: 10593 rows, 5 columns
- `data/processed/features_rdkit_descriptors.csv`: 10593 rows, 13 columns
- `results/applicability_domain_predictions.csv`: 10593 rows, 8 columns
- `results/applicability_domain_summary.csv`: 4 rows, 4 columns
- `results/combined_baseline_metrics.csv`: 4 rows, 6 columns
- `results/cross_validation_metrics.csv`: 20 rows, 9 columns
- `results/descriptor_baseline_metrics.csv`: 4 rows, 6 columns
- `results/fingerprint_baseline_metrics.csv`: 4 rows, 6 columns
- `results/ranked_candidates.csv`: 10593 rows, 22 columns
- `results/ranked_candidates_portfolio.csv`: 10593 rows, 20 columns
- `results/ranked_candidates_with_validation.csv`: 10593 rows, 22 columns
- `results/scaffold_fingerprint_metrics.csv`: 4 rows, 9 columns
- `results/top_20_candidates.csv`: 20 rows, 20 columns
- `results/top_20_diverse_candidates.csv`: 20 rows, 20 columns

## Reproducibility Status

`usable`
