# EGFR CADD/QSAR Decision Workflow

Retrospective EGFR workflow using ChEMBL activity records, RDKit features,
Morgan fingerprints, scaffold validation, uncertainty checks, structure-contact
analysis, and Vina redocking.

## Snapshot

- 26,600 raw EGFR IC50 activity rows
- 10,593 model-ready molecules
- Best scaffold-split Morgan RF: RMSE 0.871, R2 0.550
- Applicability-domain MAE: 0.513 for high-similarity chemistry vs 0.957 for low-similarity chemistry
- Assay/document-aware validation and split-conformal intervals added
- SAR/error analysis: 607 activity-cliff candidates and 387 scaffold-error rows
- Structure work: 4 EGFR co-crystals parsed, 68 ligand-contact residue rows, 5UG9 / 8AM redocking RMSD 0.968 A

## Notes

This is an existing-record benchmarking and triage workflow. It does not claim
new molecule design, therapeutic efficacy, clinical relevance, or production
readiness.
