# EGFR CADD and QSAR Decision Workflow

EGFR cheminformatics and CADD workflow built around public ChEMBL IC50
records. The project keeps the full path in one place: data curation,
descriptor and fingerprint features, QSAR baselines, scaffold aware validation,
uncertainty checks, simple drug likeness triage, and a small structure based
redocking check.

The project is retrospective. It works with existing records and known
structures for QSAR benchmarking, model-risk triage, and pose-recovery review.

## Project Contents

- ChEMBL EGFR IC50 curation from 26,600 raw activity rows.
- Molecule-level pIC50 aggregation and a 10,593-row model-ready set.
- RDKit descriptor, Morgan fingerprint, and combined-feature QSAR baselines.
- Random split, scaffold split, cross-validation, assay-aware validation, and
  document-aware validation.
- Applicability-domain analysis with max Tanimoto similarity.
- Conformal-style uncertainty checks using residual quantiles and applicability-domain proxies.
- SAR-support/error analysis, including descriptor importance, fingerprint-bit
  importance, activity-cliff candidates, and scaffold-level error summaries.
- ADMET-style and model-risk-aware ranking over existing molecules.
- An exploratory custom PyTorch GCN baseline retained as comparative benchmark evidence against the
  Morgan Random Forest baseline.
- EGFR co-crystal contact analysis for 1M17, 2ITY, 4HJO, and 5UG9.
- Retrospective Vina redocking pose-recovery audit on 5UG9 with ligand 8AM with a -9.471 kcal/mol score and 0.968 A
  pose-recovery RMSD.

## Project Snapshot

| Check | Result |
|---|---:|
| Raw ChEMBL IC50 rows | 26,600 |
| Clean molecule-level pIC50 rows | 10,834 |
| Model-ready molecules | 10,593 |
| Best random-split Morgan RF | MAE 0.516, RMSE 0.712, R2 0.719 |
| Best scaffold-split Morgan RF | MAE 0.667, RMSE 0.871, R2 0.550 |
| High-similarity applicability-domain MAE | 0.513 |
| Low-similarity applicability-domain MAE | 0.957 |
| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 A |

## Reproducing The Reports

The final reports and metrics are committed. Raw and processed ChEMBL tables are
local regenerable artifacts, following the same pattern as my antibody workflow.

To run the lightweight public report checks from existing artifacts:

```bash
make reproduce-small
make test
```

Full rebuilds require the local Python/RDKit environment and regenerated
ChEMBL-derived tables under `data/raw/` and `data/processed/`.

## Useful Outputs

- `reports/final_egfr_cadd_qsar_report.md`
- `reports/final_egfr_cv_bullets.md`
- `reports/egfr_assay_aware_validation_report.md`
- `reports/egfr_conformal_uncertainty_report.md`
- `reports/egfr_sar_interpretability_report.md`
- `reports/egfr_redocking_audit_report.md`
- `portfolio_assets/egfr_project_card.md`

Machine-readable summaries are under `reports/metrics/`.

## Interpretation Context

- ChEMBL IC50 values come from heterogeneous assays and papers.
- Scaffold and assay/document splits are more conservative than random splits.
- ADMET-style triage uses simple drug-likeness and model-risk proxy rules.
- Redocking is a retrospective co-crystal pose-recovery check.
