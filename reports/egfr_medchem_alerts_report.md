# EGFR Medicinal-Chemistry Alert Report

PAINS, Brenk, NIH, and external unwanted-substructure SMARTS alerts were used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.

## Input And Output

- Input table: `data/processed/egfr_standardized_molecules.csv`
- Output table: `data/processed/egfr_model_ready_with_medchem_alerts.csv`
- Input rows: 10,593
- Output rows: 10,593
- Invalid molecules during alert parsing: 0

## Alert Counts

- PAINS-flagged molecules: 847 (8.0%)
- Brenk-flagged molecules: 6,074 (57.3%)
- NIH-flagged molecules: 3,824 (36.1%)
- External unwanted-substructure flagged molecules: 0 (0.0%)
- Any medicinal-chemistry alert: 6,924 (65.4%)

## External SMARTS Catalog

- External CSV found: False
- External CSV path: not found
- SMARTS total: 0
- SMARTS valid: 0
- SMARTS invalid: 0

## Family Overlap

- PAINS only: 271
- Brenk only: 5,498
- External unwanted CSV only: 0
- PAINS + Brenk: 576
- Brenk + external unwanted CSV: 0
- PAINS + external unwanted CSV: 0
- All alert families: 0
- No PAINS/Brenk/external unwanted alert: 4,248

## Interpretation

The alert layer is used for risk labeling, triage penalties, and sensitivity checks. It does not prove that a molecule is inactive, toxic, or an assay artifact.
