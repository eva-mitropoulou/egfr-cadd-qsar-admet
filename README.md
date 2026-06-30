# EGFR CADD and QSAR Decision Workflow

This project builds an EGFR CADD and QSAR workflow using public ChEMBL IC50 records. I curated molecule-level activity data, standardized molecules, generated descriptor and fingerprint features, trained QSAR models, and then used the scoring workflow to rank existing EGFR inhibitor-like records for closer expert review.

The goal is not to claim new EGFR drugs. The goal is to build a transparent retrospective workflow that shows how far public activity records can support activity prediction, uncertainty checks, scaffold-aware validation, drug-likeness review, and a small structure-based sanity check.

The workflow is supported by several validation and benchmarking layers: molecule-level pIC50 aggregation, RDKit descriptor and Morgan fingerprint features, random and scaffold splits, assay-aware and document-aware validation, applicability-domain analysis, conformal-style uncertainty checks, SAR/error analysis, Lipinski/QED/PAINS/Brenk review flags, an exploratory graph neural network benchmark, co-crystal contact analysis, and retrospective Vina redocking.

## Table of Contents

- [At a Glance](#at-a-glance)
- [Project Workflow](#project-workflow)
- [Current Snapshot](#current-snapshot)
- [Selected Model View](#selected-model-view)
- [How To Read This](#how-to-read-this)
- [Scope and Limits](#scope-and-limits)
- [Reproduce](#reproduce)
- [Useful Outputs](#useful-outputs)

## At a Glance

| Part | What it does |
|---|---|
| Data curation | Starts from public ChEMBL EGFR IC50 rows and builds a molecule-level pIC50 table. |
| Molecular standardization | Standardizes structures and tracks duplicates before modeling. |
| Feature generation | Builds RDKit descriptors, Morgan fingerprints, and combined feature matrices. |
| QSAR benchmarks | Compares model families under random, scaffold, assay-aware, and document-aware splits. |
| Applicability domain | Uses max Tanimoto similarity to separate higher- and lower-support predictions. |
| Uncertainty checks | Uses residual quantiles and applicability-domain proxies to inspect interval behavior. |
| Existing-molecule review | Ranks existing molecules with predicted activity, model-risk context, QED/Lipinski checks, and PAINS/Brenk alerts. |
| Structure module | Runs co-crystal contact analysis, one retrospective redocking pose-recovery audit, and a top-5 structure-aware sanity check. |

## Project Workflow

The workflow begins with public ChEMBL records for EGFR target CHEMBL203. IC50 records are cleaned, aggregated to molecule-level pIC50 values, and filtered into a model-ready table that can support both QSAR benchmarking and later review of existing molecules.

The model-ready set contains 10,593 molecules. The project standardizes those molecules, generates RDKit descriptor and Morgan fingerprint features, and checks that feature matrices stay aligned with labels before training.

The central QSAR comparison uses Random Forest models with Morgan fingerprints as the strongest baseline in the final reports. Random splits show how well the model interpolates across a shuffled molecule table. Scaffold splits ask a harder question: whether the model can generalize across chemical scaffolds.

The project then adds more skeptical checks. Assay-aware and document-aware validation hold out assay or publication groups. Applicability-domain analysis compares performance for high-similarity and low-similarity molecules. Conformal-style uncertainty checks use residual quantiles and similarity context to estimate whether uncertainty bands behave sensibly.

The final review layer ranks existing molecules only. It combines predicted activity with model-risk context, QED and Lipinski-style drug-likeness checks, and PAINS/Brenk alert flags from RDKit filter catalogs. The output keeps a diverse top-20 review set by scaffold. The structure module adds context from EGFR co-crystals and a retrospective redocking check, but it does not turn the QSAR workflow into prospective drug discovery.

## Current Snapshot

| Check | Result |
|---|---:|
| Raw ChEMBL IC50 rows | 26,600 |
| Clean molecule-level pIC50 rows | 10,834 |
| Model-ready molecules | 10,593 |
| Best random-split Morgan RF | MAE 0.516, RMSE 0.712, R2 0.719 |
| Best scaffold-split Morgan RF | MAE 0.667, RMSE 0.871, R2 0.550 |
| Assay-group split | RMSE 1.014, R2 0.448, group overlap 0 |
| Document-group split | RMSE 1.143, R2 0.212, group overlap 0 |
| High-similarity applicability-domain MAE | 0.513 |
| Low-similarity applicability-domain MAE | 0.957 |
| 90 percent uncertainty interval coverage | 0.900 |
| Ranked existing molecules | 10,593 |
| Diverse top-20 low/medium risk count | 20/20 |
| Diverse top-20 Lipinski-clean count | 18/20 |
| PAINS/Brenk catalog status | available |
| PAINS alert count | 847 |
| Brenk alert count | 6,074 |
| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 angstrom |
| Top-5 structure sanity check | 5/5 docked; scores -8.991 to -8.386 kcal/mol; 5 warning labels |

## Selected Model View

The strongest QSAR baseline in the final report is a Morgan-fingerprint Random Forest. It performs well on the random split and remains the best scaffold-split model among the reported baselines, but the scaffold split is clearly harder than the random split.

The exploratory custom PyTorch dense GCN is kept as benchmark evidence, not as the selected scorer. In this run, it did not beat the Morgan Random Forest on scaffold RMSE.

The activity scores are used for retrospective ranking and review of existing molecules. They are not prospective potency guarantees.

## How To Read This

The random split is the optimistic benchmark. It answers how well the model works when train and test molecules are mixed by row.

The scaffold split is the more useful QSAR stress test. It asks whether the model can handle new scaffold families rather than close analogs from the same chemical neighborhood.

The assay-aware and document-aware splits are even more skeptical. They test whether performance survives changes in experimental context and publication source, which matters because public IC50 values come from heterogeneous assays.

The applicability-domain result is one of the clearest practical findings. High-similarity molecules have lower MAE than low-similarity molecules, so the ranking table carries model-risk context instead of treating every prediction as equally supported. The PAINS/Brenk flags play a different role: they do not judge potency, but they help mark molecules that need extra chemistry review before anyone over-interprets a high score.

The redocking result is a retrospective pose-recovery check on a known co-crystal case. The top-5 docking step uses that validated 5UG9 pocket to sanity-check existing ranked molecules by Vina score, pocket localization, and shared contacts with the 8AM reference ligand. In this run, all five docked molecules were kept as structure-sanity warnings because their shared 8AM contact fraction was 0.0.

## Scope and Limits

This is a retrospective public-record modeling and review project. It does not claim new EGFR inhibitors, development-stage molecules, or experimentally validated hits.

The drug-likeness and alert layer is a transparent rule-based review aid. It is not a full ADMET prediction system.

Docking is used as a pose-recovery audit on an existing co-crystal setup. Protein-ligand MD and prospective docking remain outside the default workflow.

The top-5 docking layer is a structure-aware sanity check over already-ranked existing molecules, not binding confirmation.

All rankings are for existing-molecule review and model interpretation, not for direct experimental recommendation.

## Reproduce

The final reports and metrics are committed. Raw and processed ChEMBL-derived tables are local regenerable artifacts.

Fast public checks from existing artifacts:

```bash
make reproduce-small
make test
```

Full rebuilds require the local Python/RDKit environment and regenerated ChEMBL-derived tables under `data/raw/` and `data/processed/`.

## Useful Outputs

- `reports/final_egfr_cadd_qsar_report.md`
- `reports/final_egfr_cv_bullets.md`
- `reports/egfr_assay_aware_validation_report.md`
- `reports/egfr_conformal_uncertainty_report.md`
- `reports/egfr_sar_interpretability_report.md`
- `reports/egfr_ranked_existing_molecules.csv`
- `reports/egfr_redocking_audit_report.md`
- `reports/egfr_top5_structure_sanity_report.md`
- `reports/egfr_top5_structure_sanity_table.csv`
- `portfolio_assets/egfr_project_card.md`

Machine-readable summaries are under `reports/metrics/`.
