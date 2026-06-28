# EGFR Project CV Bullets

- Built a retrospective EGFR CADD and QSAR decision workflow over 26,600 ChEMBL IC50 records and 10,593 model-ready molecules, including RDKit standardization, Morgan fingerprints, scaffold validation, applicability-domain analysis, uncertainty scoring, and ADMET-style and model-risk triage.
- Benchmarked QSAR models with random and scaffold splits; best scaffold-split model achieved MAE 0.667, RMSE 0.871, R2 0.550 while surfacing a random-to-scaffold performance drop.
- Demonstrated applicability-domain behavior: low-similarity compounds had MAE 0.957 versus 0.513 for high-similarity compounds, then used this signal in candidate triage.
- Produced a diverse top-20 existing-molecule prioritization table with 20 unique scaffolds, 20/20 low-or-medium model risk, and 18/20 Lipinski-clean molecules.
- Added structure-based EGFR co-crystal analysis across 4 parsed PDB structures with 68 ligand-contact residue rows plus a retrospective Vina redocking pose-recovery audit.
- Ran an exploratory custom PyTorch GCN baseline on the EGFR pIC50 task using NVIDIA GeForce RTX 4090; scaffold-split R2 was 0.198, underperforming the Morgan RF baseline.
