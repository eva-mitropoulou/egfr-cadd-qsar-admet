# EGFR CADD and QSAR Decision Workflow

Retrospective EGFR inhibitor-like molecule prioritization using ChEMBL, RDKit, Morgan fingerprints, scaffold validation, uncertainty, applicability-domain analysis, and ADMET-style and model-risk triage.

## Recruiter Signal

- 26,600 raw EGFR IC50 activity rows
- 10,593 model-ready molecules
- Best scaffold-split QSAR model: Random Forest, R2 0.550
- Applicability-domain MAE improved from 0.957 to 0.513 from low to high similarity
- Diverse top-20 triage: 20 unique scaffolds, 18/20 Lipinski-clean
- Structure module: 4 EGFR co-crystals parsed; 68 ligand-contact residue rows; retrospective Vina redocking pose-recovery audit `completed`
- Exploratory custom PyTorch dense GCN baseline on NVIDIA GeForce RTX 4090; scaffold R2 0.198; did not beat Morgan RF

## Positioning

A complete, model-risk-aware CADD and QSAR workflow for existing public EGFR records. No molecule generation or efficacy claim.
