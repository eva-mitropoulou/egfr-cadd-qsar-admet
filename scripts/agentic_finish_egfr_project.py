"""Autonomous controller to finish the EGFR CADD/QSAR decision workflow."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STATE_PATH = PROJECT_ROOT / "reports" / "metrics" / "agentic_egfr_state.json"
RUN_LOGS_DIR = PROJECT_ROOT / "reports" / "run_logs"


@dataclass
class Stage:
    """Controller stage definition."""

    stage_id: str
    name: str
    command: str
    expected_outputs: list[str]
    quality_gates: list[str]
    repair_strategies: list[str]
    fallback_strategies: list[str]
    max_attempts: int
    criticality: str
    fallback_command: str | None = None


@dataclass
class StageState:
    """Serializable stage execution state."""

    stage_id: str
    name: str
    command: str
    status: str = "PENDING"
    attempts: int = 0
    runtime_seconds: float = 0.0
    expected_outputs_found: list[str] = field(default_factory=list)
    expected_outputs_missing: list[str] = field(default_factory=list)
    quality_gate_results: dict[str, bool] = field(default_factory=dict)
    fallback_used: str | None = None
    warnings: list[str] = field(default_factory=list)
    next_recommended_action: str = ""


def discover_python() -> Path:
    """Find the best project Python interpreter."""
    candidates = [
        PROJECT_ROOT / ".venv" / "bin" / "python",
        PROJECT_ROOT / ".micromamba" / "envs" / "egfr-cadd" / "bin" / "python",
    ]
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return candidate
    return Path(sys.executable)


PYTHON = discover_python()


def stage_definitions() -> list[Stage]:
    """Return all EGFR finish-pipeline stages."""
    py = str(PYTHON)
    return [
        Stage(
            "phase_0",
            "Preflight",
            f"{py} scripts/preflight_egfr_project.py",
            ["reports/metrics/preflight_metrics.json", "reports/preflight_report.md"],
            ["preflight_pass"],
            ["install missing essential Python packages with pip", "create required folders"],
            ["record optional packages as unavailable and continue"],
            2,
            "core",
        ),
        Stage(
            "phase_1",
            "Project inventory and evidence audit",
            f"{py} src/data/project_completion_audit.py",
            ["reports/project_completion_audit.md", "reports/metrics/project_completion_audit.json"],
            ["audit_report_exists", "critical_artifacts_identified"],
            ["repair missing output directories", "rerun audit"],
            ["mark missing artifacts and continue if core data exists"],
            2,
            "core",
        ),
        Stage(
            "phase_2",
            "Data provenance and ChEMBL activity audit",
            f"{py} src/data/audit_egfr_chembl_data.py",
            ["reports/egfr_data_provenance_audit.md", "reports/metrics/egfr_data_provenance_audit.json"],
            ["raw_row_count_reported", "clean_molecule_count_reported", "model_ready_count_large_enough"],
            ["run existing fetch/clean scripts if raw or clean files are missing"],
            ["use existing model-ready dataset if raw audit is degraded"],
            2,
            "critical",
        ),
        Stage(
            "phase_3",
            "Molecular standardization and preparation audit",
            f"{py} src/chem/standardize_molecules.py && {py} src/chem/audit_molecular_representations.py",
            [
                "data/processed/egfr_standardized_molecules.csv",
                "reports/molecular_standardization_report.md",
                "reports/metrics/molecular_standardization_metrics.json",
            ],
            ["standardized_count_reported", "invalid_count_reported", "standardization_policy_reported"],
            ["fallback to canonical SMILES standardization if MolStandardize is unavailable"],
            ["mark advanced standardization as degraded and continue"],
            2,
            "core",
        ),
        Stage(
            "phase_4",
            "Feature generation audit",
            f"{py} src/features/build_rdkit_descriptors.py && {py} src/features/build_morgan_fingerprints.py && {py} src/features/build_combined_features.py",
            [
                "data/processed/features_rdkit_descriptors.csv",
                "data/processed/features_morgan_fingerprints.npz",
                "data/processed/features_combined_descriptors_morgan.npz",
                "reports/feature_generation_report.md",
                "reports/metrics/feature_generation_metrics.json",
            ],
            ["descriptor_shape_reported", "fingerprint_shape_reported", "feature_label_alignment"],
            ["rebuild standardized molecules", "drop constant descriptors if needed"],
            ["continue with Morgan fingerprints only if descriptor build degrades"],
            2,
            "core",
        ),
        Stage(
            "phase_5",
            "Matched QSAR benchmark",
            f"{py} src/models/train_qsar_matched_benchmarks.py",
            [
                "reports/qsar_matched_benchmark_report.md",
                "reports/metrics/qsar_matched_benchmark_metrics.json",
                "reports/figures/random_vs_scaffold_performance.png",
                "reports/figures/predicted_vs_observed_scaffold.png",
                "reports/figures/residuals_by_split.png",
                "models/egfr_primary_full_model.joblib",
            ],
            ["random_and_scaffold_reported", "scaffold_counts_reported", "no_scaffold_or_molecule_leakage", "primary_model_available"],
            ["fallback to RDKit Murcko scaffold splitter", "continue if one noncritical model fails"],
            ["use Morgan Random Forest only if matched benchmark partially degrades"],
            2,
            "critical",
        ),
        Stage(
            "phase_6",
            "Applicability domain and nearest-neighbor reliability",
            f"{py} src/validation/applicability_domain.py",
            [
                "reports/applicability_domain_report.md",
                "reports/metrics/applicability_domain_metrics.json",
                "reports/egfr_applicability_domain_predictions.csv",
                "reports/figures/error_vs_similarity.png",
                "reports/figures/similarity_distribution.png",
            ],
            ["low_high_similarity_reported", "ad_prediction_table_exists"],
            ["recompute Morgan fingerprints", "fallback to existing applicability-domain predictions"],
            ["mark applicability domain degraded if bins cannot be computed"],
            2,
            "core",
        ),
        Stage(
            "phase_7",
            "Uncertainty and calibration",
            f"{py} src/models/uncertainty_calibration.py",
            [
                "reports/egfr_uncertainty_calibration_report.md",
                "reports/metrics/egfr_uncertainty_calibration_metrics.json",
                "reports/egfr_uncertainty_predictions.csv",
                "reports/figures/calibration_curve.png",
                "reports/figures/prediction_interval_coverage.png",
                "reports/figures/uncertainty_vs_error.png",
            ],
            ["uncertainty_score_reported", "coverage_reported"],
            ["fallback to RF ensemble variance only", "skip conformal interval if it fails"],
            ["mark conformal calibration degraded and continue"],
            2,
            "core",
        ),
        Stage(
            "phase_8",
            "ADMET-style and synthetic-feasibility triage",
            f"{py} src/triage/rank_egfr_existing_molecules.py",
            [
                "reports/egfr_candidate_triage_report.md",
                "reports/egfr_ranked_existing_molecules.csv",
                "reports/metrics/egfr_candidate_triage_metrics.json",
                "reports/figures/egfr_triage_funnel.png",
                "reports/figures/predicted_activity_vs_risk.png",
                "reports/figures/applicability_domain_by_rank.png",
            ],
            ["ranked_table_exists", "diverse_top20_reported", "model_risk_flags_included"],
            ["disable PAINS/Brenk if RDKit catalogs fail", "mark synthetic accessibility unavailable"],
            ["continue with QED/Lipinski/model-risk triage"],
            2,
            "core",
        ),
        Stage(
            "phase_9",
            "Structure-based EGFR module",
            f"{py} src/structure/retrieve_egfr_structures.py && {py} src/structure/prepare_structure_metadata.py && {py} src/structure/redocking_validation.py && {py} src/structure/interaction_fingerprints.py",
            [
                "reports/egfr_structure_metadata_report.md",
                "reports/egfr_interaction_fingerprint_report.md",
                "reports/metrics/egfr_structure_module_metrics.json",
                "reports/figures/egfr_structure_workflow.png",
            ],
            ["structure_availability_reported", "redocking_status_reported"],
            ["skip docking if executable unavailable", "use metadata-only report if fetch fails"],
            ["mark structure module degraded and continue"],
            2,
            "optional",
        ),
        Stage(
            "phase_10",
            "Small-molecule GNN benchmark",
            f"{py} src/models/train_gnn_benchmark.py",
            [
                "reports/egfr_gnn_benchmark_report.md",
                "reports/metrics/egfr_gnn_benchmark_metrics.json",
                "reports/figures/rf_vs_gnn_scaffold_split.png",
            ],
            ["gnn_status_reported", "gpu_status_reported"],
            ["try available Chemprop/DeepChem/PyG backend if installed"],
            ["mark GNN benchmark degraded and continue"],
            1,
            "optional",
        ),
        Stage(
            "phase_11",
            "Retrospective active-learning simulation",
            f"{py} src/triage/active_learning_simulation.py",
            [
                "reports/egfr_active_learning_report.md",
                "reports/metrics/egfr_active_learning_metrics.json",
                "reports/figures/egfr_active_learning_discovery_curve.png",
                "reports/figures/egfr_active_learning_scaffold_diversity.png",
            ],
            ["random_baseline_included", "informed_strategies_included", "best_strategy_reported"],
            ["reduce active-learning rounds if runtime fails", "continue with fewer strategies if needed"],
            ["mark active learning degraded if simulation cannot complete"],
            1,
            "core",
        ),
        Stage(
            "phase_12",
            "Optional protein-ligand MD bridge",
            f"{py} src/structure/protein_ligand_md_next_steps.py",
            ["reports/egfr_protein_ligand_md_next_steps.md"],
            ["md_next_steps_reported"],
            ["do not run full MD automatically"],
            ["document future work and continue"],
            1,
            "optional",
        ),
        Stage(
            "phase_13",
            "CLI and demo",
            f"{py} src/app/predict_egfr_cli.py --create-example examples/example_smiles.csv --output reports/example_predictions.csv",
            ["src/app/predict_egfr_cli.py", "reports/example_predictions.csv", "reports/egfr_cli_demo_report.md"],
            ["cli_output_exists", "cli_report_exists"],
            ["create example from existing project molecules", "use primary model if present"],
            ["mark CLI demo degraded and continue"],
            2,
            "core",
        ),
        Stage(
            "phase_14",
            "Final EGFR report and portfolio assets",
            f"{py} src/analysis/build_final_egfr_report.py",
            [
                "reports/final_egfr_cadd_qsar_report.md",
                "reports/final_egfr_cv_bullets.md",
                "reports/final_egfr_interview_talking_points.md",
                "portfolio_assets/egfr_project_card.md",
                "reports/metrics/final_egfr_project_status.json",
                "reports/metrics/agentic_egfr_state.json",
            ],
            ["final_report_exists", "cv_bullets_exists", "portfolio_card_exists"],
            ["rebuild missing final metrics from available reports"],
            ["mark optional modules degraded in final report"],
            2,
            "critical",
        ),
    ]


def load_json(path: Path) -> dict[str, Any]:
    """Load JSON safely."""
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def output_status(expected_outputs: list[str]) -> tuple[list[str], list[str]]:
    """Return found and missing expected outputs."""
    found = []
    missing = []
    for relative in expected_outputs:
        path = PROJECT_ROOT / relative
        if path.exists() and (not path.is_file() or path.stat().st_size > 0):
            found.append(relative)
        else:
            missing.append(relative)
    return found, missing


def gate_results(stage: Stage) -> dict[str, bool]:
    """Evaluate known quality gates."""
    results: dict[str, bool] = {}
    preflight = load_json(PROJECT_ROOT / "reports" / "metrics" / "preflight_metrics.json")
    audit = load_json(PROJECT_ROOT / "reports" / "metrics" / "project_completion_audit.json")
    provenance = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_data_provenance_audit.json")
    std = load_json(PROJECT_ROOT / "reports" / "metrics" / "molecular_standardization_metrics.json")
    feat = load_json(PROJECT_ROOT / "reports" / "metrics" / "feature_generation_metrics.json")
    bench = load_json(PROJECT_ROOT / "reports" / "metrics" / "qsar_matched_benchmark_metrics.json")
    ad = load_json(PROJECT_ROOT / "reports" / "metrics" / "applicability_domain_metrics.json")
    unc = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_uncertainty_calibration_metrics.json")
    triage = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_candidate_triage_metrics.json")
    structure = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_structure_module_metrics.json")
    gnn = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_gnn_benchmark_metrics.json")
    active = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_active_learning_metrics.json")

    checks: dict[str, Callable[[], bool]] = {
        "preflight_pass": lambda: bool(preflight.get("preflight_pass")),
        "audit_report_exists": lambda: (PROJECT_ROOT / "reports" / "project_completion_audit.md").exists(),
        "critical_artifacts_identified": lambda: "missing_critical_artifacts" in audit,
        "raw_row_count_reported": lambda: provenance.get("raw_activity_row_count", 0) > 0,
        "clean_molecule_count_reported": lambda: provenance.get("clean_pIC50_molecule_count", 0) > 0,
        "model_ready_count_large_enough": lambda: provenance.get("model_ready_molecule_count", 0) >= 1000,
        "standardized_count_reported": lambda: std.get("standardized_rows", 0) >= 1000,
        "invalid_count_reported": lambda: "invalid_molecule_count" in std,
        "standardization_policy_reported": lambda: "tautomer_policy" in std and "charge_policy" in std,
        "descriptor_shape_reported": lambda: feat.get("rdkit_descriptors", {}).get("descriptor_rows", 0) >= 1000,
        "fingerprint_shape_reported": lambda: feat.get("morgan_fingerprints", {}).get("fingerprint_rows", 0) >= 1000,
        "feature_label_alignment": lambda: bool(feat.get("combined_features", {}).get("feature_label_alignment")),
        "random_and_scaffold_reported": lambda: any(row.get("split") == "random_split" for row in bench.get("matched_benchmark_rows", []))
        and any(row.get("split") == "scaffold_split" for row in bench.get("matched_benchmark_rows", [])),
        "scaffold_counts_reported": lambda: bench.get("scaffold_split", {}).get("train_scaffolds", 0) > 0
        and bench.get("scaffold_split", {}).get("test_scaffolds", 0) > 0,
        "no_scaffold_or_molecule_leakage": lambda: bench.get("scaffold_split", {}).get("scaffold_overlap", 999) == 0
        and bench.get("scaffold_split", {}).get("standardized_molecule_overlap", 999) == 0,
        "primary_model_available": lambda: (PROJECT_ROOT / "models" / "egfr_primary_full_model.joblib").exists(),
        "low_high_similarity_reported": lambda: ad.get("low_similarity_mae") is not None and ad.get("high_similarity_mae") is not None,
        "ad_prediction_table_exists": lambda: (PROJECT_ROOT / "reports" / "egfr_applicability_domain_predictions.csv").exists(),
        "uncertainty_score_reported": lambda: bool(unc.get("uncertainty_score")),
        "coverage_reported": lambda: unc.get("coverage_90") is not None,
        "ranked_table_exists": lambda: (PROJECT_ROOT / "reports" / "egfr_ranked_existing_molecules.csv").exists(),
        "diverse_top20_reported": lambda: triage.get("diverse_top20_unique_scaffolds", 0) >= 20,
        "model_risk_flags_included": lambda: bool(triage.get("risk_counts")),
        "structure_availability_reported": lambda: "available_structures" in structure,
        "redocking_status_reported": lambda: "redocking_status" in structure,
        "gnn_status_reported": lambda: "gnn_status" in gnn,
        "gpu_status_reported": lambda: "cuda_available" in gnn,
        "random_baseline_included": lambda: "random" in active.get("strategies", []),
        "informed_strategies_included": lambda: len([s for s in active.get("strategies", []) if s != "random"]) >= 2,
        "best_strategy_reported": lambda: bool(active.get("best_strategy")),
        "md_next_steps_reported": lambda: (PROJECT_ROOT / "reports" / "egfr_protein_ligand_md_next_steps.md").exists(),
        "cli_output_exists": lambda: (PROJECT_ROOT / "reports" / "example_predictions.csv").exists(),
        "cli_report_exists": lambda: (PROJECT_ROOT / "reports" / "egfr_cli_demo_report.md").exists(),
        "final_report_exists": lambda: (PROJECT_ROOT / "reports" / "final_egfr_cadd_qsar_report.md").exists(),
        "cv_bullets_exists": lambda: (PROJECT_ROOT / "reports" / "final_egfr_cv_bullets.md").exists(),
        "portfolio_card_exists": lambda: (PROJECT_ROOT / "portfolio_assets" / "egfr_project_card.md").exists(),
    }

    for gate in stage.quality_gates:
        func = checks.get(gate)
        results[gate] = bool(func()) if func is not None else False
    return results


def run_command(command: str, stage_id: str, attempt: int) -> tuple[int, float, Path]:
    """Run a stage command and save stdout/stderr to a log file."""
    RUN_LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = RUN_LOGS_DIR / f"{stage_id}_attempt_{attempt}.log"
    start = time.time()
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        shell=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    runtime = time.time() - start
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    return completed.returncode, runtime, log_path


def repair_from_log(log_path: Path) -> list[str]:
    """Repair obvious package/path issues from logs."""
    warnings: list[str] = []
    text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    (PROJECT_ROOT / "reports" / "metrics").mkdir(parents=True, exist_ok=True)
    RUN_LOGS_DIR.mkdir(parents=True, exist_ok=True)

    match = re.search(r"ModuleNotFoundError: No module named '([^']+)'", text)
    if match:
        package = match.group(1).split(".")[0]
        pip_name = {"sklearn": "scikit-learn", "rdkit": "rdkit-pypi"}.get(package, package)
        subprocess.run(
            [str(PYTHON), "-m", "pip", "install", pip_name],
            cwd=PROJECT_ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        warnings.append(f"attempted_pip_install:{pip_name}")

    if "No such file or directory" in text or "FileNotFoundError" in text:
        for folder in [
            "scripts",
            "src",
            "reports",
            "reports/metrics",
            "reports/figures",
            "reports/run_logs",
            "models",
            "data/processed",
        ]:
            (PROJECT_ROOT / folder).mkdir(parents=True, exist_ok=True)
        warnings.append("ensured_standard_directories")

    subprocess.run(
        [str(PYTHON), "-m", "py_compile", *[str(path) for path in PROJECT_ROOT.glob("src/**/*.py")], *[str(path) for path in PROJECT_ROOT.glob("scripts/*.py")]],
        cwd=PROJECT_ROOT,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    warnings.append("ran_py_compile")
    return warnings


def save_state(state: dict[str, Any]) -> None:
    """Save agentic controller state."""
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")


def run_stage(stage: Stage, state: dict[str, Any]) -> StageState:
    """Run one stage with retries and fallback/degraded handling."""
    stage_state = StageState(stage.stage_id, stage.name, stage.command)
    total_runtime = 0.0
    for attempt in range(1, stage.max_attempts + 1):
        stage_state.attempts = attempt
        return_code, runtime, log_path = run_command(stage.command, stage.stage_id, attempt)
        total_runtime += runtime
        stage_state.runtime_seconds = round(total_runtime, 3)
        found, missing = output_status(stage.expected_outputs)
        gates = gate_results(stage)
        stage_state.expected_outputs_found = found
        stage_state.expected_outputs_missing = missing
        stage_state.quality_gate_results = gates

        if return_code == 0 and not missing and all(gates.values()):
            stage_state.status = "PASS"
            stage_state.next_recommended_action = "continue"
            return stage_state

        stage_state.warnings.append(f"attempt_{attempt}_failed_log:{log_path.relative_to(PROJECT_ROOT)}")
        if attempt < stage.max_attempts:
            stage_state.warnings.extend(repair_from_log(log_path))
            state["stages"][stage.stage_id] = asdict(stage_state)
            save_state(state)

    if stage.fallback_command:
        return_code, runtime, log_path = run_command(stage.fallback_command, stage.stage_id + "_fallback", 1)
        total_runtime += runtime
        stage_state.runtime_seconds = round(total_runtime, 3)
        stage_state.fallback_used = stage.fallback_command
        found, missing = output_status(stage.expected_outputs)
        gates = gate_results(stage)
        stage_state.expected_outputs_found = found
        stage_state.expected_outputs_missing = missing
        stage_state.quality_gate_results = gates
        if return_code == 0 and not missing:
            stage_state.status = "FALLBACK_PASS" if all(gates.values()) else "DEGRADED"
        else:
            stage_state.status = "DEGRADED"
        stage_state.warnings.append(f"fallback_log:{log_path.relative_to(PROJECT_ROOT)}")
    else:
        stage_state.status = "DEGRADED" if stage.criticality == "optional" else "FAIL"

    stage_state.next_recommended_action = (
        "inspect logs and repair core stage" if stage_state.status == "FAIL" else "continue with documented degradation"
    )
    return stage_state


def core_blocked(stage: Stage, stage_state: StageState) -> bool:
    """Return whether controller must stop."""
    if stage_state.status in {"PASS", "FALLBACK_PASS", "DEGRADED"} and stage.criticality != "critical":
        return False
    if stage_state.status in {"PASS", "FALLBACK_PASS"}:
        return False
    if stage.stage_id in {"phase_2", "phase_5", "phase_14"}:
        return True
    return False


def final_summary(state: dict[str, Any]) -> dict[str, Any]:
    """Build concise final summary for terminal output."""
    final_status_payload = load_json(PROJECT_ROOT / "reports" / "metrics" / "final_egfr_project_status.json")
    benchmark = load_json(PROJECT_ROOT / "reports" / "metrics" / "qsar_matched_benchmark_metrics.json")
    applicability = load_json(PROJECT_ROOT / "reports" / "metrics" / "applicability_domain_metrics.json")
    uncertainty = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_uncertainty_calibration_metrics.json")
    triage = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_candidate_triage_metrics.json")
    structure = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_structure_module_metrics.json")
    gnn = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_gnn_benchmark_metrics.json")
    active = load_json(PROJECT_ROOT / "reports" / "metrics" / "egfr_active_learning_metrics.json")

    def best(split: str) -> dict[str, Any]:
        rows = [row for row in benchmark.get("matched_benchmark_rows", []) if row.get("split") == split]
        return sorted(rows, key=lambda row: (row.get("RMSE", 999), row.get("MAE", 999)))[0] if rows else {}

    random_best = final_status_payload.get("best_random_split_model") or best("random_split")
    scaffold_best = final_status_payload.get("best_scaffold_split_model") or best("scaffold_split")
    drop = final_status_payload.get("scaffold_performance_drop_R2")
    if drop is None and random_best and scaffold_best:
        drop = random_best.get("R2", 0) - scaffold_best.get("R2", 0)

    required_outputs = [
        "reports/project_completion_audit.md",
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
        "reports/metrics/agentic_egfr_state.json",
    ]
    missing_required = [path for path in required_outputs if not (PROJECT_ROOT / path).exists()]
    blocked = any(stage_data.get("status") == "FAIL" for stage_data in state["stages"].values())
    if blocked:
        final_status = "BLOCKED"
    elif missing_required or gnn.get("gnn_status", "").startswith("degraded") or structure.get("structure_module_status", "").endswith("degraded"):
        final_status = "DONE_WITH_WARNINGS"
    else:
        final_status = "DONE"

    return {
        "FINAL_STATUS": final_status,
        "status_table": [
            {"stage_id": key, "status": value.get("status"), "attempts": value.get("attempts")}
            for key, value in state["stages"].items()
        ],
        "best_random_split_model": random_best,
        "best_scaffold_split_model": scaffold_best,
        "scaffold_performance_drop": drop,
        "applicability_domain_finding": {
            "low_similarity_mae": applicability.get("low_similarity_mae"),
            "high_similarity_mae": applicability.get("high_similarity_mae"),
        },
        "uncertainty_calibration_status": {
            "method": uncertainty.get("uncertainty_score"),
            "coverage_90": uncertainty.get("coverage_90"),
            "conformal_status": uncertainty.get("conformal_status"),
        },
        "triage_table_size": triage.get("ranked_molecule_count"),
        "structure_module_status": structure.get("structure_module_status"),
        "gnn_status": gnn.get("gnn_status"),
        "active_learning_best_strategy": active.get("best_strategy"),
        "cli_demo_status": (PROJECT_ROOT / "reports" / "egfr_cli_demo_report.md").exists(),
        "first_files_to_inspect": [
            "reports/final_egfr_cadd_qsar_report.md",
            "reports/final_egfr_cv_bullets.md",
            "portfolio_assets/egfr_project_card.md",
            "reports/metrics/agentic_egfr_state.json",
        ],
    }


def print_final(summary: dict[str, Any]) -> None:
    """Print only the final requested summary fields."""
    print(f"FINAL_STATUS: {summary['FINAL_STATUS']}")
    print("status table:")
    for row in summary["status_table"]:
        print(f"{row['stage_id']}\t{row['status']}\tattempts={row['attempts']}")

    def fmt_model(row: dict[str, Any]) -> str:
        if not row:
            return "unavailable"
        return f"{row.get('model')} | MAE={row.get('MAE'):.3f} RMSE={row.get('RMSE'):.3f} R2={row.get('R2'):.3f}"

    print(f"best random-split model: {fmt_model(summary['best_random_split_model'])}")
    print(f"best scaffold-split model: {fmt_model(summary['best_scaffold_split_model'])}")
    drop = summary["scaffold_performance_drop"]
    print(f"scaffold performance drop: {drop:.3f}" if drop is not None else "scaffold performance drop: unavailable")
    ad = summary["applicability_domain_finding"]
    print(
        "applicability-domain finding: "
        f"low-similarity MAE={ad.get('low_similarity_mae'):.3f}, "
        f"high-similarity MAE={ad.get('high_similarity_mae'):.3f}"
    )
    unc = summary["uncertainty_calibration_status"]
    print(
        "uncertainty/calibration status: "
        f"{unc.get('method')} | coverage_90={unc.get('coverage_90')} | {unc.get('conformal_status')}"
    )
    print(f"triage table size: {summary['triage_table_size']}")
    print(f"structure module status: {summary['structure_module_status']}")
    print(f"GNN status: {summary['gnn_status']}")
    print(f"active-learning best strategy: {summary['active_learning_best_strategy']}")
    print(f"CLI/demo status: {summary['cli_demo_status']}")
    print("first files to inspect:")
    for path in summary["first_files_to_inspect"]:
        print(path)


def main() -> None:
    """Run all stages."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--overnight", action="store_true", help="Run full noninteractive finish pipeline.")
    args = parser.parse_args()

    state: dict[str, Any] = {
        "project_root": str(PROJECT_ROOT),
        "python_executable": str(PYTHON),
        "overnight": bool(args.overnight),
        "stages": {},
        "warnings": [],
    }
    save_state(state)

    for stage in stage_definitions():
        stage_state = run_stage(stage, state)
        state["stages"][stage.stage_id] = asdict(stage_state)
        save_state(state)
        if core_blocked(stage, stage_state):
            state["FINAL_STATUS"] = "BLOCKED"
            state["warnings"].append(f"blocked_at:{stage.stage_id}")
            save_state(state)
            break

    summary = final_summary(state)
    state.update(summary)
    save_state(state)
    print_final(summary)


if __name__ == "__main__":
    main()
