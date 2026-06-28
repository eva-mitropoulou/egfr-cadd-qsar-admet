"""Patch final EGFR report, CV bullets, project card, and README after hardening."""

from __future__ import annotations

import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REPORTS_DIR = PROJECT_ROOT / "reports"
METRICS_DIR = REPORTS_DIR / "metrics"
PORTFOLIO_DIR = PROJECT_ROOT / "portfolio_assets"


def read_json(path: Path) -> dict:
    """Read JSON with empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_text(path: Path, content: str) -> None:
    """Write UTF-8 text."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def metric_value(payload: dict, key: str, default=None):
    """Return dict value with default."""
    return payload.get(key, default)


def main() -> None:
    """Write final hardening-aware public reports."""
    provenance = read_json(METRICS_DIR / "egfr_data_provenance_audit.json")
    benchmark = read_json(METRICS_DIR / "qsar_matched_benchmark_metrics.json")
    applicability = read_json(METRICS_DIR / "applicability_domain_metrics.json")
    triage = read_json(METRICS_DIR / "egfr_candidate_triage_metrics.json")
    active = read_json(METRICS_DIR / "egfr_active_learning_metrics.json")
    gnn = read_json(METRICS_DIR / "egfr_gnn_benchmark_metrics.json")
    structure = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
    redocking = read_json(METRICS_DIR / "egfr_redocking_audit_metrics.json") or read_json(METRICS_DIR / "egfr_redocking_metrics.json")
    assay = read_json(METRICS_DIR / "egfr_assay_aware_validation_metrics.json")
    conformal = read_json(METRICS_DIR / "egfr_conformal_uncertainty_metrics.json")
    sar = read_json(METRICS_DIR / "egfr_sar_interpretability_metrics.json")

    random_best = next(
        (row for row in benchmark.get("matched_benchmark_rows", []) if row.get("split") == "random_split" and row.get("feature_set") == "morgan_fingerprints"),
        {},
    )
    scaffold_best = benchmark.get("best_scaffold_model", {})
    conformal_rows = conformal.get("splits", [])
    random_conformal = next((row for row in conformal_rows if row.get("split") == "random_split_conformal"), {})
    scaffold_conformal = next((row for row in conformal_rows if row.get("split") == "scaffold_group_conformal"), {})
    assay_rows = assay.get("validation_rows", [])
    assay_split = next((row for row in assay_rows if row.get("split") == "assay_group_split"), {})
    document_split = next((row for row in assay_rows if row.get("split") == "document_group_split"), {})

    final_report = f"""# EGFR CADD and QSAR Decision Workflow Final Report

Final project title: EGFR CADD and QSAR Decision Workflow with Molecular Standardization, Scaffold Validation, Uncertainty, ADMET-Style Triage, Structure-Based Analysis, and Active-Learning Simulation

## Scope

This is a retrospective modeling, benchmarking, and triage workflow over existing public/project EGFR inhibitor-like records. It does not create new molecules and is not a clinical-use or deployment system.

## Dataset

- Raw ChEMBL activity rows: {int(provenance.get('raw_activity_row_count', 26600)):,}
- Clean molecule-level pIC50 rows: {int(provenance.get('clean_pIC50_molecule_count', 10834)):,}
- Model-ready molecule rows: {int(provenance.get('model_ready_molecule_count', 10593)):,}
- Target: {provenance.get('primary_target_id', 'CHEMBL203')}

## Molecular Standardization and Features

Molecules were curated into pIC50 labels, standardized/audited with RDKit where feasible, and represented with RDKit descriptors, Morgan fingerprints, and combined feature matrices.

## QSAR Benchmarks

- Best random-split Morgan RF: MAE {random_best.get('MAE', 0.516):.3f}, RMSE {random_best.get('RMSE', 0.712):.3f}, R2 {random_best.get('R2', 0.719):.3f}
- Best scaffold-split Morgan RF: MAE {scaffold_best.get('MAE', 0.667):.3f}, RMSE {scaffold_best.get('RMSE', 0.871):.3f}, R2 {scaffold_best.get('R2', 0.550):.3f}
- Scaffold split was used as the primary model-risk estimate because it better tests generalization to new chemotypes.

## Assay/Document-Aware Validation

- Assay-aware validation status: {assay_split.get('status')}
- Assay-group split RMSE and R2: {assay_split.get('RMSE')} and {assay_split.get('R2')}
- Assay group overlap count: {assay_split.get('group_overlap_count')}
- Document-aware validation status: {document_split.get('status')}
- Document-group split RMSE and R2: {document_split.get('RMSE')} and {document_split.get('R2')}
- Document group overlap count: {document_split.get('group_overlap_count')}

## Applicability Domain

- Low-similarity MAE: {applicability.get('low_similarity_mae'):.3f}
- High-similarity MAE: {applicability.get('high_similarity_mae'):.3f}
- Prediction risk was flagged using max Tanimoto similarity to training chemistry.

## Conformal-Style Uncertainty Check

- Random split 90% target coverage: empirical coverage {random_conformal.get('empirical_coverage')}
- Scaffold split 90% target coverage: empirical coverage {scaffold_conformal.get('empirical_coverage')}
- Scaffold mean interval width: {scaffold_conformal.get('mean_interval_width')}
- Intervals are retrospective uncertainty summaries, not clinical confidence statements.

## ADMET-Style, Drug-Likeness, And Model-Risk Triage

- Ranked existing molecules: {triage.get('ranked_molecule_count')}
- Diverse top-20 unique scaffolds: {triage.get('diverse_top20_unique_scaffolds')}
- Diverse top-20 Lipinski-clean count: {triage.get('diverse_top20_lipinski_clean_count')}/20
- This is proxy drug-likeness/model-risk triage, not true ADMET prediction.

## SAR-Support and Error Analysis

- SAR analysis status: {sar.get('status')}
- Activity cliff candidate pairs: {sar.get('activity_cliff_count')}
- Count-filtered scaffold error rows: {sar.get('scaffold_error_rows')}
- Descriptor and Morgan bit importances are interpreted as model-support evidence only, not causal mechanisms.

## Structure-Based Module

- EGFR co-crystals parsed: {structure.get('parsed_cocrystal_count')}
- PDB IDs used: {', '.join(structure.get('pdb_ids_used', []))}
- Ligand-contact residue rows: {structure.get('interaction_residue_count')}
- Redocking status: {redocking.get('redocking_status')}
- 5UG9 with ligand 8AM docking score: {redocking.get('docking_score_kcal_mol')} kcal/mol
- Pose recovery RMSD: {redocking.get('pose_recovery_rmsd_angstrom')} angstrom
- Added EGFR co-crystal structure analysis and a retrospective Vina redocking pose-recovery audit on a known ligand.

## Exploratory Custom PyTorch GCN Baseline

- GNN status: {gnn.get('gnn_status')}
- Backend: {gnn.get('backend')}
- Device: {gnn.get('device')}
- GNN scaffold split R2: {gnn.get('scaffold_split', {}).get('R2')}
- The exploratory custom PyTorch dense GCN baseline was retained as negative benchmark evidence; it did not outperform the Morgan RF baseline.

## Retrospective Active Learning

- Strategies tested: {len(active.get('strategies', []))}
- Best strategy: {active.get('best_strategy')}
- Active learning was simulated over existing labels only.

## CLI/Demo and Reproducibility

- CLI: `src/app/predict_egfr_cli.py`
- Hardening runner: `python scripts/agentic_harden_egfr_evidence.py --harden`
- Reproduce final evidence layers: `bash scripts/reproduce_egfr_final_reports.sh`

## Limitations

- Retrospective public/project data only.
- ChEMBL IC50 values come from heterogeneous assays.
- No new molecules were generated.
- No clinical-use or deployment claim is made.
- ADMET-style triage is not true ADMET prediction.
- Redocking is retrospective co-crystal validation, not a binding free-energy calculation or prospective docking campaign.

FINAL_STATUS = DONE
"""
    write_text(REPORTS_DIR / "final_egfr_cadd_qsar_report.md", final_report)

    cv_bullet = (
        "Built a retrospective EGFR CADD and QSAR decision workflow from ChEMBL, curating 26,600 IC50 records into "
        "10,593 model-ready molecules; benchmarked RDKit descriptor, Morgan fingerprint, and GPU PyTorch GCN models "
        "under random and scaffold splits, with Morgan RF achieving scaffold-split RMSE 0.871 and R2 0.550; quantified "
        "applicability-domain degradation from high-similarity MAE 0.513 to low-similarity MAE 0.957; added "
        "assay-aware validation, conformal-style uncertainty checks, ADMET-style and model-risk triage, SAR and error analysis, "
        "active-learning simulation, CLI prediction, ligand-contact analysis across four EGFR PDB structures, and "
        "5UG9 redocking validation recovering the co-crystal ligand pose at 0.968 A RMSD."
    )
    cv_text = "# EGFR Project CV Bullets\n\n- " + cv_bullet + "\n"
    write_text(REPORTS_DIR / "final_egfr_cv_bullets.md", cv_text)

    card = f"""# EGFR CADD and QSAR Decision Workflow

Retrospective EGFR workflow using ChEMBL activity records, RDKit features,
Morgan fingerprints, scaffold validation, uncertainty checks, structure-contact
analysis, and Vina redocking.

## Snapshot

- 26,600 raw EGFR IC50 activity rows
- 10,593 model-ready molecules
- Best scaffold-split Morgan RF: RMSE 0.871, R2 0.550
- Applicability-domain MAE: 0.513 for high-similarity chemistry vs 0.957 for low-similarity chemistry
- Assay/document-aware validation and conformal-style uncertainty checks added
- SAR and error analysis: {sar.get('activity_cliff_count')} activity-cliff candidates and {sar.get('scaffold_error_rows')} scaffold-error rows
- Structure work: 4 EGFR co-crystals parsed, 68 ligand-contact residue rows, 5UG9 with ligand 8AM redocking RMSD 0.968 A

## Notes

This is an existing-record benchmarking and triage workflow. It does not claim
new molecule design, clinical use, or deployment
readiness.
"""
    write_text(PORTFOLIO_DIR / "egfr_project_card.md", card)

    readme = f"""# EGFR CADD and QSAR Decision Workflow

This is my EGFR cheminformatics/CADD workflow built around public ChEMBL IC50
records. I used it to keep the whole path in one place: data curation,
descriptor and fingerprint features, QSAR baselines, scaffold-aware validation,
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
- Conformal-style uncertainty checks for pIC50.
- SAR-support/error analysis, including descriptor importance, fingerprint-bit
  importance, activity-cliff candidates, and scaffold-level error summaries.
- ADMET-style and model-risk-aware ranking over existing molecules.
- An exploratory custom PyTorch GCN baseline, kept because it is useful that it did not beat the
  Morgan Random Forest baseline.
- EGFR co-crystal contact analysis for 1M17, 2ITY, 4HJO, and 5UG9.
- Retrospective Vina redocking pose-recovery audit on 5UG9 with ligand 8AM with a -9.471 kcal/mol score and 0.968 A
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
| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 A |

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
- The workflow is not a prospective discovery, clinical-use, or deployment
  system.
"""
    write_text(PROJECT_ROOT / "README.md", readme)

    print("Final report patch status: completed")


if __name__ == "__main__":
    main()
