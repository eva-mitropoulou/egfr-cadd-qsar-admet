# EGFR CADD and QSAR Decision Workflow Final Report

Final project title: EGFR CADD and QSAR Decision Workflow with Molecular Standardization, Scaffold Validation, Uncertainty, ADMET-Style Triage, Structure-Based Analysis, and Active-Learning Simulation

## Project Role

This is a retrospective modeling, benchmarking, and triage workflow over existing public/project EGFR inhibitor-like records.

## Dataset

- Raw ChEMBL activity rows: 26,600
- Clean molecule-level pIC50 rows: 10,834
- Model-ready molecule rows: 10,593
- Target: CHEMBL203

## Molecular Standardization

- Standardized molecules: 10,593
- Invalid molecules: 0
- Duplicate standardized molecules: 42
- MolStandardize available: True

## Feature Generation

RDKit descriptors, Morgan fingerprints, and combined descriptor and fingerprint matrices were generated and checked for label alignment.

## QSAR Benchmarks

- Best random-split model: Random Forest with MAE 0.516, RMSE 0.712, R2 0.719
- Best scaffold-split model: Random Forest with MAE 0.667, RMSE 0.871, R2 0.550
- Scaffold R2 drop relative to random split: 0.168

## Assay/Document-Aware Validation

- Assay-group split: RMSE 1.014, R2 0.448, group overlap 0
- Document-group split: RMSE 1.143, R2 0.212, group overlap 0

## Applicability Domain

- Low-similarity MAE: 0.957
- High-similarity MAE: 0.513
- Out-of-domain count: 149

## Conformal-Style Uncertainty Check

- Uncertainty score: Applicability-domain reliability proxy: 1 - max Tanimoto similarity
- Uncertainty-error Spearman correlation: 0.273
- 90 percent interval coverage: 0.900

This is a retrospective uncertainty proxy using residual quantiles and applicability-domain context.

## ADMET-Style, Drug-Likeness, And Model-Risk Triage

- Ranked existing molecules: 10,593
- Diverse top-20 unique scaffolds: 20
- Diverse top-20 low/medium risk count: 20/20
- Diverse top-20 Lipinski-clean count: 18/20

This is transparent drug-likeness and model-risk triage over existing molecules.

## Structure-Based Module

- Structure module status: completed_redocking
- Available structures: 4
- Parsed co-crystals with ligand: 4
- PDB IDs used: 1M17, 2ITY, 4HJO, 5UG9
- Redocking audit status: completed
- Pose recovery RMSD: 0.968 angstrom
- Overlay artifact status: overlay_figure_created
- Interaction fingerprint status: completed_heuristic_contacts
- Binding-site contact residue rows: 68

The structure module completed co-crystal retrieval, binding-site interaction analysis, and a retrospective Vina redocking pose-recovery audit for one prepared EGFR co-crystal.

## Exploratory Custom PyTorch GCN Baseline

- GNN status: completed
- Backend: custom_pytorch_dense_gcn
- CUDA available: True
- Device: NVIDIA GeForce RTX 4090
- GNN random split: MAE 0.886, RMSE 1.115, R2 0.310
- GNN scaffold split: MAE 0.909, RMSE 1.149, R2 0.198
- GNN beat Morgan RF on scaffold RMSE: False
The exploratory custom PyTorch dense GCN baseline is retained as comparative graph-model evidence against the Morgan Random Forest baseline in this run.

## Retrospective Active Learning

- Strategies tested: 6
- Best strategy: applicability_domain_aware_high_score
- Best final recovery fraction: 0.3996229971724788

## CLI/Demo

A CLI script is available at `src/app/predict_egfr_cli.py` and writes prediction outputs without raw SMILES.

## Interpretation Context

- Public/project ChEMBL IC50 values come from heterogeneous assays.
- Scaffold, assay, document, and applicability-domain checks carry the main validation interpretation.
- Docking and protein-ligand MD are optional/future structure-based extensions.
- The redocking result is a retrospective pose-recovery sanity check.

FINAL_STATUS = DONE
