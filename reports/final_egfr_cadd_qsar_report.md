# EGFR CADD/QSAR Decision Workflow Final Report

Final project title: EGFR CADD/QSAR Decision Workflow with Molecular Standardization, Scaffold Validation, Uncertainty, ADMET-Style Triage, Structure-Based Analysis, and Active-Learning Simulation

## Scope

This is a retrospective modeling, benchmarking, and triage workflow over existing public/project EGFR inhibitor-like records. It does not create new molecules, claim therapeutic efficacy, or claim production-grade deployment.

## Dataset

- Raw ChEMBL activity rows: 26,600
- Clean molecule-level pIC50 rows: 10,834
- Model-ready molecule rows: 10,593
- Target: CHEMBL203

## Molecular Standardization and Features

Molecules were curated into pIC50 labels, standardized/audited with RDKit where feasible, and represented with RDKit descriptors, Morgan fingerprints, and combined feature matrices.

## QSAR Benchmarks

- Best random-split Morgan RF: MAE 0.516, RMSE 0.712, R2 0.719
- Best scaffold-split Morgan RF: MAE 0.667, RMSE 0.871, R2 0.550
- Scaffold split was used as the primary model-risk estimate because it better tests generalization to new chemotypes.

## Assay/Document-Aware Validation

- Assay-aware validation status: completed
- Assay-group split RMSE/R2: 1.0135992949202917 / 0.4480171002230632
- Assay group overlap count: 0
- Document-aware validation status: completed
- Document-group split RMSE/R2: 1.1425119156890893 / 0.21240743631437198
- Document group overlap count: 0

## Applicability Domain

- Low-similarity MAE: 0.957
- High-similarity MAE: 0.513
- Prediction risk was flagged using max Tanimoto similarity to training chemistry.

## Split-Conformal Uncertainty

- Random split 90% target coverage: empirical coverage 0.9023124115148655
- Scaffold split 90% target coverage: empirical coverage 0.9313771888310459
- Scaffold mean interval width: 3.0943794820870933
- Intervals are retrospective uncertainty summaries, not clinical confidence statements.

## ADMET-Style / Drug-Likeness / Model-Risk Triage

- Ranked existing molecules: 10593
- Diverse top-20 unique scaffolds: 20
- Diverse top-20 Lipinski-clean count: 18/20
- This is proxy drug-likeness/model-risk triage, not true ADMET prediction.

## SAR-Support and Error Analysis

- SAR analysis status: completed
- Activity cliff candidate pairs: 607
- Count-filtered scaffold error rows: 387
- Descriptor and Morgan bit importances are interpreted as model-support evidence only, not causal mechanisms.

## Structure-Based Module

- EGFR co-crystals parsed: 4
- PDB IDs used: 1M17, 2ITY, 4HJO, 5UG9
- Ligand-contact residue rows: 68
- Redocking status: completed_redocking
- 5UG9 / 8AM docking score: -9.471 kcal/mol
- Pose recovery RMSD: 0.968 angstrom
- Added EGFR co-crystal structure analysis and Vina redocking validation on a known ligand.

## GNN Benchmark

- GNN status: completed
- Backend: custom_pytorch_dense_gcn
- Device: NVIDIA GeForce RTX 4090
- GNN scaffold split R2: 0.19800341794433707
- The GPU PyTorch dense GCN was run honestly and did not outperform the Morgan RF baseline.

## Retrospective Active Learning

- Strategies tested: 6
- Best strategy: applicability_domain_aware_high_score
- Active learning was simulated over existing labels only.

## CLI/Demo and Reproducibility

- CLI: `src/app/predict_egfr_cli.py`
- Hardening runner: `python scripts/agentic_harden_egfr_evidence.py --harden`
- Reproduce final evidence layers: `bash scripts/reproduce_egfr_final_reports.sh`

## Limitations

- Retrospective public/project data only.
- ChEMBL IC50 values come from heterogeneous assays.
- No new molecules were generated.
- No therapeutic efficacy is claimed.
- No production-grade deployment is claimed.
- ADMET-style triage is not true ADMET prediction.
- Redocking is retrospective co-crystal validation, not a binding free-energy calculation or prospective docking campaign.

FINAL_STATUS = DONE_WITH_WARNINGS
