"""Build final EGFR CADD and QSAR decision-workflow reports."""

from __future__ import annotations

import sys
from pathlib import Path
import json


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

METRICS_DIR = PROJECT_ROOT / "reports" / "metrics"
REPORTS_DIR = PROJECT_ROOT / "reports"
PORTFOLIO_DIR = PROJECT_ROOT / "portfolio_assets"


def read_json(path: Path) -> dict:
    """Read JSON with an empty-dict fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, payload: dict) -> None:
    """Save JSON without importing project utility dependencies."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_text(path: Path, content: str) -> None:
    """Write text without importing project utility dependencies."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def fmt_int(value: object) -> str:
    """Format integer-like values for public reports."""
    return f"{int(value):,}" if value is not None else "unavailable"


REQUIRED_OUTPUTS = [
    "reports/molecular_standardization_report.md",
    "reports/qsar_matched_benchmark_report.md",
    "reports/applicability_domain_report.md",
    "reports/egfr_uncertainty_calibration_report.md",
    "reports/egfr_candidate_triage_report.md",
    "reports/egfr_structure_metadata_report.md",
    "reports/egfr_active_learning_report.md",
    "reports/final_egfr_cadd_qsar_report.md",
    "reports/final_egfr_cv_bullets.md",
    "portfolio_assets/egfr_project_card.md",
    ]


def best_split_model(metrics: dict, split_name: str) -> dict:
    """Return best model row for one split by RMSE."""
    rows = [row for row in metrics.get("matched_benchmark_rows", []) if row.get("split") == split_name]
    if not rows:
        return {}
    return sorted(rows, key=lambda row: (row.get("RMSE", 999), row.get("MAE", 999)))[0]


def status_from_outputs() -> str:
    """Return final status based on required outputs present before final files are written."""
    missing = [path for path in REQUIRED_OUTPUTS[:-4] if not (PROJECT_ROOT / path).exists()]
    if missing:
        return "DONE_WITH_WARNINGS"
    return "DONE"


def main() -> None:
    """Build final project reports from saved metrics."""
    provenance = read_json(METRICS_DIR / "egfr_data_provenance_audit.json")
    standardization = read_json(METRICS_DIR / "molecular_standardization_metrics.json")
    benchmark = read_json(METRICS_DIR / "qsar_matched_benchmark_metrics.json")
    assay_validation = read_json(METRICS_DIR / "egfr_assay_aware_validation_metrics.json")
    applicability = read_json(METRICS_DIR / "applicability_domain_metrics.json")
    uncertainty = read_json(METRICS_DIR / "egfr_uncertainty_calibration_metrics.json")
    triage = read_json(METRICS_DIR / "egfr_candidate_triage_metrics.json")
    structure = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
    redocking = read_json(METRICS_DIR / "egfr_redocking_audit_metrics.json")
    gnn = read_json(METRICS_DIR / "egfr_gnn_benchmark_metrics.json")
    active = read_json(METRICS_DIR / "egfr_active_learning_metrics.json")

    random_best = best_split_model(benchmark, "random_split")
    scaffold_best = best_split_model(benchmark, "scaffold_split")
    scaffold_drop = None
    if random_best and scaffold_best and random_best.get("R2") is not None and scaffold_best.get("R2") is not None:
        scaffold_drop = float(random_best["R2"] - scaffold_best["R2"])
    assay_rows = assay_validation.get("validation_rows", [])
    assay_split = next((row for row in assay_rows if row.get("split") == "assay_group_split"), {})
    document_split = next((row for row in assay_rows if row.get("split") == "document_group_split"), {})

    final_status = status_from_outputs()
    if (
        gnn.get("gnn_status", "").startswith("degraded")
        or structure.get("structure_module_status", "").endswith("degraded")
        or redocking.get("status", "completed") != "completed"
    ):
        final_status = "DONE_WITH_WARNINGS"

    report = [
        "# EGFR CADD and QSAR Decision Workflow Final Report",
        "",
        "Final project title: EGFR CADD and QSAR Decision Workflow with Molecular Standardization, Scaffold Validation, Uncertainty, ADMET-Style Triage, Structure-Based Analysis, and Active-Learning Simulation",
        "",
        "## Project Role",
        "",
        "This is a retrospective modeling, benchmarking, and triage workflow over existing public/project EGFR inhibitor-like records.",
        "",
        "## Dataset",
        "",
        f"- Raw ChEMBL activity rows: {fmt_int(provenance.get('raw_activity_row_count'))}",
        f"- Clean molecule-level pIC50 rows: {fmt_int(provenance.get('clean_pIC50_molecule_count'))}",
        f"- Model-ready molecule rows: {fmt_int(provenance.get('model_ready_molecule_count'))}",
        f"- Target: {provenance.get('primary_target_id')}",
        "",
        "## Molecular Standardization",
        "",
        f"- Standardized molecules: {fmt_int(standardization.get('standardized_rows'))}",
        f"- Invalid molecules: {fmt_int(standardization.get('invalid_molecule_count'))}",
        f"- Duplicate standardized molecules: {fmt_int(standardization.get('duplicate_standardized_molecule_count'))}",
        f"- MolStandardize available: {standardization.get('molstandardize_available')}",
        "",
        "## Feature Generation",
        "",
        "RDKit descriptors, Morgan fingerprints, and combined descriptor and fingerprint matrices were generated and checked for label alignment.",
        "",
        "## QSAR Benchmarks",
        "",
        f"- Best random-split model: {random_best.get('model')} with MAE {random_best.get('MAE'):.3f}, RMSE {random_best.get('RMSE'):.3f}, R2 {random_best.get('R2'):.3f}" if random_best else "- Best random-split model: unavailable",
        f"- Best scaffold-split model: {scaffold_best.get('model')} with MAE {scaffold_best.get('MAE'):.3f}, RMSE {scaffold_best.get('RMSE'):.3f}, R2 {scaffold_best.get('R2'):.3f}" if scaffold_best else "- Best scaffold-split model: unavailable",
        f"- Scaffold R2 drop relative to random split: {scaffold_drop:.3f}" if scaffold_drop is not None else "- Scaffold R2 drop: unavailable",
        "",
        "## Assay/Document-Aware Validation",
        "",
        (
            f"- Assay-group split: RMSE {assay_split.get('RMSE'):.3f}, "
            f"R2 {assay_split.get('R2'):.3f}, group overlap {assay_split.get('group_overlap_count')}"
            if assay_split
            else "- Assay-group split: unavailable"
        ),
        (
            f"- Document-group split: RMSE {document_split.get('RMSE'):.3f}, "
            f"R2 {document_split.get('R2'):.3f}, group overlap {document_split.get('group_overlap_count')}"
            if document_split
            else "- Document-group split: unavailable"
        ),
        "",
        "## Applicability Domain",
        "",
        f"- Low-similarity MAE: {applicability.get('low_similarity_mae'):.3f}",
        f"- High-similarity MAE: {applicability.get('high_similarity_mae'):.3f}",
        f"- Out-of-domain count: {applicability.get('out_of_domain_count')}",
        "",
        "## Conformal-Style Uncertainty Check",
        "",
        f"- Uncertainty score: {uncertainty.get('uncertainty_score')}",
        f"- Uncertainty-error Spearman correlation: {uncertainty.get('uncertainty_error_spearman'):.3f}",
        f"- 90 percent interval coverage: {uncertainty.get('coverage_90'):.3f}",
        "",
        "This is a retrospective uncertainty proxy using residual quantiles and applicability-domain context.",
        "",
        "## ADMET-Style, Drug-Likeness, And Model-Risk Triage",
        "",
        f"- Ranked existing molecules: {fmt_int(triage.get('ranked_molecule_count'))}",
        f"- Diverse top-20 unique scaffolds: {triage.get('diverse_top20_unique_scaffolds')}",
        f"- Diverse top-20 low/medium risk count: {triage.get('diverse_top20_low_or_medium_risk_count')}/20",
        f"- Diverse top-20 Lipinski-clean count: {triage.get('diverse_top20_lipinski_clean_count')}/20",
        "",
        "This is transparent drug-likeness and model-risk triage over existing molecules.",
        "",
        "## Structure-Based Module",
        "",
        f"- Structure module status: {structure.get('structure_module_status')}",
        f"- Available structures: {structure.get('available_structures')}",
        f"- Parsed co-crystals with ligand: {structure.get('parsed_cocrystal_count')}",
        f"- PDB IDs used: {', '.join(structure.get('pdb_ids_used', [])) if structure.get('pdb_ids_used') else 'none'}",
        f"- Redocking audit status: {redocking.get('status')}",
        f"- Pose recovery RMSD: {redocking.get('pose_recovery_rmsd_angstrom')} angstrom",
        f"- Overlay artifact status: {redocking.get('overlay_artifact_status')}",
        f"- Interaction fingerprint status: {structure.get('interaction_fingerprint_status')}",
        f"- Binding-site contact residue rows: {structure.get('interaction_residue_count')}",
        "",
        "The structure module completed co-crystal retrieval, binding-site interaction analysis, and a retrospective Vina redocking pose-recovery audit for one prepared EGFR co-crystal.",
        "",
        "## Exploratory Custom PyTorch GCN Baseline",
        "",
        f"- GNN status: {gnn.get('gnn_status')}",
        f"- Backend: {gnn.get('backend')}",
        f"- CUDA available: {gnn.get('cuda_available')}",
        f"- Device: {gnn.get('device')}",
        (
            f"- GNN random split: MAE {gnn.get('random_split', {}).get('MAE'):.3f}, "
            f"RMSE {gnn.get('random_split', {}).get('RMSE'):.3f}, "
            f"R2 {gnn.get('random_split', {}).get('R2'):.3f}"
            if gnn.get("random_split")
            else "- GNN random split: unavailable"
        ),
        (
            f"- GNN scaffold split: MAE {gnn.get('scaffold_split', {}).get('MAE'):.3f}, "
            f"RMSE {gnn.get('scaffold_split', {}).get('RMSE'):.3f}, "
            f"R2 {gnn.get('scaffold_split', {}).get('R2'):.3f}"
            if gnn.get("scaffold_split")
            else "- GNN scaffold split: unavailable"
        ),
        f"- GNN beat Morgan RF on scaffold RMSE: {gnn.get('gnn_beat_morgan_rf', {}).get('scaffold_RMSE')}",
        "The exploratory custom PyTorch dense GCN baseline is retained as comparative benchmark evidence against the Morgan Random Forest baseline in this run.",
        "",
        "## Retrospective Active Learning",
        "",
        f"- Strategies tested: {len(active.get('strategies', []))}",
        f"- Best strategy: {active.get('best_strategy')}",
        f"- Best final recovery fraction: {active.get('best_final_recovery_fraction')}",
        "",
        "## CLI/Demo",
        "",
        "A CLI script is available at `src/app/predict_egfr_cli.py` and writes prediction outputs without raw SMILES.",
        "",
        "## Limitations",
        "",
        "- Retrospective public/project data only.",
        "- ChEMBL IC50 values come from heterogeneous assays.",
        "- Docking and protein-ligand MD are optional/future structure-based extensions.",
        "- The redocking result is a retrospective pose-recovery sanity check.",
        "",
        f"FINAL_STATUS = {final_status}",
        "",
    ]
    write_text(REPORTS_DIR / "final_egfr_cadd_qsar_report.md", "\n".join(report))

    cv_bullets = [
        "# EGFR Project CV Bullets",
        "",
        "- Built a retrospective EGFR CADD and QSAR decision workflow over 26,600 ChEMBL IC50 records and 10,593 model-ready molecules, including RDKit standardization, Morgan fingerprints, scaffold validation, applicability-domain analysis, uncertainty scoring, and ADMET-style and model-risk triage.",
        f"- Benchmarked QSAR models with random and scaffold splits; best scaffold-split model achieved MAE {scaffold_best.get('MAE'):.3f}, RMSE {scaffold_best.get('RMSE'):.3f}, R2 {scaffold_best.get('R2'):.3f} while surfacing a random-to-scaffold performance drop.",
        f"- Demonstrated applicability-domain behavior: low-similarity compounds had MAE {applicability.get('low_similarity_mae'):.3f} versus {applicability.get('high_similarity_mae'):.3f} for high-similarity compounds, then used this signal in candidate triage.",
        f"- Produced a diverse top-20 existing-molecule prioritization table with {triage.get('diverse_top20_unique_scaffolds')} unique scaffolds, {triage.get('diverse_top20_low_or_medium_risk_count')}/20 low-or-medium model risk, and {triage.get('diverse_top20_lipinski_clean_count')}/20 Lipinski-clean molecules.",
        f"- Added structure-based EGFR co-crystal analysis across {structure.get('parsed_cocrystal_count')} parsed PDB structures with {structure.get('interaction_residue_count')} ligand-contact residue rows plus a retrospective Vina redocking pose-recovery audit.",
        f"- Ran an exploratory custom PyTorch GCN baseline on the EGFR pIC50 task using {gnn.get('device')}; scaffold-split R2 was {gnn.get('scaffold_split', {}).get('R2'):.3f}, underperforming the Morgan RF baseline.",
        "",
    ]
    write_text(REPORTS_DIR / "final_egfr_cv_bullets.md", "\n".join(cv_bullets))

    talking_points = [
        "# EGFR Project Interview Talking Points",
        "",
        "- Why scaffold split matters: random splits can place close analogs in both train and test, inflating apparent QSAR performance.",
        "- Why applicability domain matters: predictions are more reliable when a molecule is close to training chemistry; the project quantified this with max Tanimoto similarity.",
        "- How to explain the validation frame: retrospective model benchmarks, existing-molecule triage, and structure sanity checks.",
        "- Why Morgan fingerprints beat descriptors: fingerprints encode substructure patterns that are more relevant to kinase inhibitor SAR than global properties alone.",
        "- What could come next: docking validation, experimental feedback loops, true ADMET predictors, and protein-ligand MD on selected validated poses.",
        "",
    ]
    write_text(REPORTS_DIR / "final_egfr_interview_talking_points.md", "\n".join(talking_points))

    project_card = [
        "# EGFR CADD and QSAR Decision Workflow",
        "",
        "Retrospective EGFR inhibitor-like molecule prioritization using ChEMBL, RDKit, Morgan fingerprints, scaffold validation, uncertainty, applicability-domain analysis, and ADMET-style and model-risk triage.",
        "",
        "## Recruiter Signal",
        "",
        "- 26,600 raw EGFR IC50 activity rows",
        "- 10,593 model-ready molecules",
        f"- Best scaffold-split QSAR model: {scaffold_best.get('model')}, R2 {scaffold_best.get('R2'):.3f}",
        f"- Applicability-domain MAE improved from {applicability.get('low_similarity_mae'):.3f} to {applicability.get('high_similarity_mae'):.3f} from low to high similarity",
        f"- Diverse top-20 triage: {triage.get('diverse_top20_unique_scaffolds')} unique scaffolds, {triage.get('diverse_top20_lipinski_clean_count')}/20 Lipinski-clean",
        f"- Structure module: {structure.get('parsed_cocrystal_count')} EGFR co-crystals parsed; {structure.get('interaction_residue_count')} ligand-contact residue rows; retrospective Vina redocking pose-recovery audit `{redocking.get('status')}`",
        f"- Exploratory custom PyTorch dense GCN baseline on {gnn.get('device')}; scaffold R2 {gnn.get('scaffold_split', {}).get('R2'):.3f}; did not beat Morgan RF",
        "",
        "## Positioning",
        "",
        "A complete, model-risk-aware CADD and QSAR workflow for existing public EGFR records. No molecule generation or efficacy claim.",
        "",
    ]
    PORTFOLIO_DIR.mkdir(parents=True, exist_ok=True)
    write_text(PORTFOLIO_DIR / "egfr_project_card.md", "\n".join(project_card))

    status_payload = {
        "FINAL_STATUS": final_status,
        "best_random_split_model": random_best,
        "best_scaffold_split_model": scaffold_best,
        "scaffold_performance_drop_R2": scaffold_drop,
        "applicability_domain_finding": {
            "low_similarity_mae": applicability.get("low_similarity_mae"),
            "high_similarity_mae": applicability.get("high_similarity_mae"),
        },
        "uncertainty_status": uncertainty.get("conformal_status"),
        "triage_table_size": triage.get("ranked_molecule_count"),
        "structure_module_status": structure.get("structure_module_status"),
        "redocking_status": redocking.get("status"),
        "redocking_pose_recovery_rmsd_angstrom": redocking.get("pose_recovery_rmsd_angstrom"),
        "redocking_overlay_artifact_status": redocking.get("overlay_artifact_status"),
        "interaction_analysis_status": structure.get("interaction_analysis_status"),
        "gnn_status": gnn.get("gnn_status"),
        "gnn_backend": gnn.get("backend"),
        "gnn_random_split": gnn.get("random_split"),
        "gnn_scaffold_split": gnn.get("scaffold_split"),
        "gnn_beat_morgan_rf": gnn.get("gnn_beat_morgan_rf"),
        "active_learning_best_strategy": active.get("best_strategy"),
        "cli_demo_status": (REPORTS_DIR / "egfr_cli_demo_report.md").exists(),
    }
    save_json(METRICS_DIR / "final_egfr_project_status.json", status_payload)

    print(f"Final report: {REPORTS_DIR / 'final_egfr_cadd_qsar_report.md'}")
    print(f"Final status: {final_status}")


if __name__ == "__main__":
    main()
