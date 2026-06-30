# EGFR CADD and QSAR Decision Workflow

Retrospective EGFR inhibitor-like molecule prioritization using ChEMBL, RDKit, Morgan fingerprints, scaffold validation, uncertainty, applicability-domain analysis, medicinal-chemistry alerts, and ADMET-style/model-risk triage.

## Snapshot

- 26,600 raw EGFR IC50 activity rows
- 10,593 model-ready molecules
- Best scaffold-split QSAR model: Morgan Random Forest, R2 0.550
- Applicability-domain MAE changed from 0.957 at low similarity to 0.513 at high similarity
- PAINS-flagged molecules: 847 (8.0%)
- Brenk-flagged molecules: 6,074 (57.3%)
- Top-20 without medchem alerts: 20/20
- Structure module: four EGFR co-crystals parsed; retrospective 5UG9/8AM Vina redocking pose-recovery RMSD 0.968 A
- Exploratory custom PyTorch dense GCN benchmark did not beat the Morgan RF baseline

- Top-5 docking score sanity check: 5/5 molecules docked; Vina score range -8.991 to -8.386 kcal/mol

## Positioning

A complete, model-risk-aware CADD and QSAR workflow for existing public EGFR records. No molecule generation or efficacy claim.

PAINS, Brenk, and external unwanted-substructure SMARTS alerts were used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.
