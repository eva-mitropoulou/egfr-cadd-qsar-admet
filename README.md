# EGFR CADD/QSAR Decision Workflow

This is my EGFR cheminformatics/CADD workflow built around public ChEMBL IC50
records. I used it to keep the whole path in one place: data curation,
descriptor/fingerprint features, QSAR baselines, scaffold-aware validation,
uncertainty checks, simple drug-likeness triage, and a small structure-based
redocking check.

The project is retrospective. It works with existing records and known
structures; it does not generate molecules or claim that any compound is a drug
candidate.

## What Is In Here

- ChEMBL EGFR IC50 curation from 26,600 raw activity rows.
- Molecule-level pIC50 aggregation and a 10,593-row model-ready set.
- RDKit descriptor, Morgan fingerprint, and combined-feature QSAR baselines.
- Random split, scaffold split, cross-validation, assay-aware validation, and
  document-aware validation.
- Applicability-domain analysis with max Tanimoto similarity.
- Split-conformal prediction intervals for pIC50.
- SAR-support/error analysis, including descriptor importance, fingerprint-bit
  importance, activity-cliff candidates, and scaffold-level error summaries.
- ADMET-style and model-risk-aware ranking over existing molecules.
- A PyTorch GCN benchmark, kept because it is useful that it did not beat the
  Morgan Random Forest baseline.
- EGFR co-crystal contact analysis for 1M17, 2ITY, 4HJO, and 5UG9.
- Vina redocking on 5UG9 / 8AM with a -9.471 kcal/mol score and 0.968 A
  pose-recovery RMSD.

## Current Snapshot

| Check | Result |
|---|---:|
| Raw ChEMBL IC50 rows | 26,600 |
| Clean molecule-level pIC50 rows | 10,834 |
| Model-ready molecules | 10,593 |
| Best random-split Morgan RF | MAE 0.516, RMSE 0.712, R2 0.719 |
| Best scaffold-split Morgan RF | MAE 0.667, RMSE 0.871, R2 0.550 |
| High-similarity applicability-domain MAE | 0.513 |
| Low-similarity applicability-domain MAE | 0.957 |
| Redocking case | 5UG9 / 8AM, RMSD 0.968 A |

## Reproducing The Reports

The final reports and metrics are committed. Raw and processed ChEMBL tables are
local artifacts and are not committed, following the same pattern as my antibody
workflow.

To rerun only the lightweight report/evidence hardening stages from existing
artifacts:

```bash
source .venv/bin/activate 2>/dev/null || true
python scripts/agentic_harden_egfr_evidence.py --harden
```

Or:

```bash
bash scripts/reproduce_egfr_final_reports.sh
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
- `reports/egfr_final_hardening_status.md`
- `portfolio_assets/egfr_project_card.md`

Machine-readable summaries are under `reports/metrics/`.

## Caveats

- ChEMBL IC50 values come from heterogeneous assays and papers.
- Scaffold and assay/document splits are more conservative than random splits.
- ADMET-style triage here means simple drug-likeness/model-risk rules, not true
  ADMET prediction.
- Redocking is a retrospective co-crystal check, not a binding free-energy
  calculation.
- The workflow does not claim prospective discovery, therapeutic efficacy,
  clinical relevance, or production-grade prediction.
