"""Agentic hardening runner for final EGFR evidence layers."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PYTHON = str(Path(sys.executable).resolve())

REPORTS_DIR = PROJECT_ROOT / "reports"
RUN_LOGS_DIR = REPORTS_DIR / "run_logs"
METRICS_DIR = REPORTS_DIR / "metrics"
STATE_PATH = METRICS_DIR / "agentic_egfr_hardening_state.json"


@dataclass
class Stage:
    """One hardening runner stage."""

    stage_id: str
    name: str
    command: list[str]
    expected_outputs: list[str]
    quality_gates: list[str]
    repair_strategies: list[str]
    fallback_strategies: list[list[str]] = field(default_factory=list)
    max_attempts: int = 2
    criticality: str = "noncritical"


def ensure_dirs() -> None:
    """Create hardening output directories."""
    for relative in [
        "reports",
        "reports/metrics",
        "reports/figures",
        "reports/run_logs",
        "scripts",
        "tests",
        "src/validation",
        "src/analysis",
        "src/structure",
    ]:
        (PROJECT_ROOT / relative).mkdir(parents=True, exist_ok=True)


def stages() -> list[Stage]:
    """Return hardening stages."""
    py = PYTHON
    return [
        Stage(
            stage_id="stage_0",
            name="Preflight",
            command=[py, "scripts/preflight_egfr_hardening.py"],
            expected_outputs=["reports/metrics/egfr_hardening_preflight.json"],
            quality_gates=["project root resolves", "model-ready data checked", "core package availability checked"],
            repair_strategies=["create missing folders", "compile changed scripts"],
            max_attempts=2,
            criticality="critical",
        ),
        Stage(
            stage_id="stage_1",
            name="Artifact and column inventory",
            command=[py, "src/analysis/egfr_hardening_inventory.py"],
            expected_outputs=["reports/egfr_hardening_inventory.md", "reports/metrics/egfr_hardening_inventory.json"],
            quality_gates=["model-ready table found", "pIC50 and molecule identifier found", "assay/document metadata availability reported"],
            repair_strategies=["repair path/import issues", "compile changed scripts"],
            max_attempts=2,
            criticality="critical",
        ),
        Stage(
            stage_id="stage_2",
            name="Assay-aware and document-aware validation",
            command=[py, "src/validation/assay_aware_validation.py"],
            expected_outputs=[
                "reports/egfr_assay_aware_validation_report.md",
                "reports/metrics/egfr_assay_aware_validation_metrics.json",
                "reports/figures/random_vs_scaffold_vs_assay_document_split.png",
            ],
            quality_gates=["metadata availability reported", "group overlap reported", "assay/document validation attempted"],
            repair_strategies=["repair column-name mappings", "fallback to degraded metadata report"],
            max_attempts=2,
            criticality="noncritical",
        ),
        Stage(
            stage_id="stage_3",
            name="Split-conformal uncertainty intervals",
            command=[py, "src/validation/conformal_uncertainty.py"],
            expected_outputs=[
                "reports/egfr_conformal_uncertainty_report.md",
                "reports/metrics/egfr_conformal_uncertainty_metrics.json",
                "reports/figures/conformal_coverage_by_split.png",
                "reports/figures/conformal_interval_width_vs_similarity.png",
                "reports/figures/conformal_interval_width_distribution.png",
            ],
            quality_gates=["90% target coverage reported", "empirical coverage reported", "mean interval width reported"],
            repair_strategies=["repair split/index alignment", "fallback to existing uncertainty proxy"],
            max_attempts=2,
            criticality="noncritical",
        ),
        Stage(
            stage_id="stage_4",
            name="SAR and interpretable error analysis",
            command=[py, "src/analysis/egfr_sar_error_analysis.py"],
            expected_outputs=[
                "reports/egfr_sar_interpretability_report.md",
                "reports/metrics/egfr_sar_interpretability_metrics.json",
                "reports/egfr_activity_cliffs.csv",
                "reports/egfr_scaffold_error_table.csv",
                "reports/figures/top_descriptor_importances.png",
                "reports/figures/scaffold_level_error.png",
                "reports/figures/activity_cliff_similarity_vs_delta.png",
            ],
            quality_gates=["descriptor importance reported", "activity cliff count reported", "scaffold-level error table created"],
            repair_strategies=["repair index alignment", "reduce cliff search output if needed"],
            max_attempts=2,
            criticality="noncritical",
        ),
        Stage(
            stage_id="stage_5",
            name="Redocking evidence hardening",
            command=[py, "src/structure/harden_redocking_evidence.py"],
            expected_outputs=[
                "reports/egfr_redocking_audit_report.md",
                "reports/metrics/egfr_redocking_audit_metrics.json",
                "reports/egfr_redocking_report.md",
            ],
            quality_gates=["redocking status preserved", "score and RMSD preserved", "overlay artifact or script created"],
            repair_strategies=["fallback to overlay script if figure backend fails"],
            max_attempts=2,
            criticality="noncritical",
        ),
        Stage(
            stage_id="stage_6",
            name="Final report and portfolio patch",
            command=[py, "src/analysis/patch_final_egfr_reports.py"],
            expected_outputs=[
                "reports/final_egfr_cadd_qsar_report.md",
                "reports/final_egfr_cv_bullets.md",
                "portfolio_assets/egfr_project_card.md",
                "README.md",
            ],
            quality_gates=["new evidence layers mentioned", "claims remain retrospective", "redocking score/RMSD included"],
            repair_strategies=["repair missing metric fallbacks", "compile changed scripts"],
            max_attempts=2,
            criticality="noncritical",
        ),
        Stage(
            stage_id="stage_7",
            name="Smoke tests and reproducibility",
            command=[py, "scripts/run_smoke_tests.py"],
            expected_outputs=["reports/metrics/egfr_smoke_test_metrics.json", "scripts/reproduce_egfr_final_reports.sh"],
            quality_gates=["required artifacts exist", "metrics JSON parses", "public markdown hygiene checked"],
            repair_strategies=["repair stale wording", "fallback to direct test runner"],
            fallback_strategies=[[py, "scripts/run_smoke_tests.py"]],
            max_attempts=2,
            criticality="noncritical",
        ),
        Stage(
            stage_id="stage_8",
            name="Final hardening status",
            command=[py, "src/analysis/build_final_hardening_status.py"],
            expected_outputs=["reports/egfr_final_hardening_status.md", "reports/metrics/egfr_final_hardening_status.json"],
            quality_gates=["final status written", "degraded items reported exactly"],
            repair_strategies=["repair missing metric fallback"],
            max_attempts=2,
            criticality="noncritical",
        ),
    ]


def run_command(command: list[str], log_path: Path) -> tuple[int, float]:
    """Run a command, saving stdout/stderr to log."""
    start = time.time()
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    runtime = time.time() - start
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    return completed.returncode, runtime


def compile_scripts() -> None:
    """Compile hardening scripts as a repair check."""
    py = PYTHON
    scripts = [
        "scripts/preflight_egfr_hardening.py",
        "src/analysis/egfr_hardening_inventory.py",
        "src/validation/assay_aware_validation.py",
        "src/validation/conformal_uncertainty.py",
        "src/analysis/egfr_sar_error_analysis.py",
        "src/structure/harden_redocking_evidence.py",
        "src/analysis/patch_final_egfr_reports.py",
        "src/analysis/build_final_hardening_status.py",
    ]
    subprocess.run([py, "-m", "py_compile", *scripts], cwd=PROJECT_ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)


def output_status(expected_outputs: list[str]) -> tuple[list[str], list[str]]:
    """Return found and missing expected outputs."""
    found = []
    missing = []
    for output in expected_outputs:
        path = PROJECT_ROOT / output
        if path.exists() and path.stat().st_size > 0:
            found.append(output)
        else:
            missing.append(output)
    return found, missing


def save_state(state: dict) -> None:
    """Persist runner state."""
    STATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")


def load_json(path: Path) -> dict:
    """Read JSON with fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def print_final_summary(state: dict) -> None:
    """Print only the requested final hardening summary."""
    status = load_json(METRICS_DIR / "egfr_final_hardening_status.json")
    updated_files = sorted(
        {
            output
            for stage in state.get("stages", {}).values()
            for output in stage.get("expected_outputs_found", [])
        }
    )
    print(f"FINAL_HARDENING_STATUS: {status.get('FINAL_HARDENING_STATUS')}")
    print(f"assay-aware validation status: {status.get('assay_aware_validation_status')}")
    print(f"document-aware validation status: {status.get('document_aware_validation_status')}")
    print(f"conformal uncertainty status: {status.get('conformal_uncertainty_status')}")
    print(f"SAR analysis status: {status.get('SAR_analysis_status')}")
    print(f"redocking audit status: {status.get('redocking_audit_status')}")
    print(f"overlay artifact status: {status.get('overlay_artifact_status')}")
    print(f"smoke-test status: {status.get('smoke_test_status')}")
    print("files updated:")
    for path in updated_files:
        print(f"- {path}")
    degraded = status.get("degraded_items", [])
    print("degraded items:")
    if degraded:
        for item in degraded:
            print(f"- {item.get('stage')}: {item.get('reason')}")
    else:
        print("- none")


def main() -> None:
    """Run hardening stages."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--harden", action="store_true", help="Run EGFR hardening stages.")
    args = parser.parse_args()
    if not args.harden:
        raise SystemExit("Use --harden")

    ensure_dirs()
    state = {"started_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "stages": {}, "next_recommended_action": None}
    for stage in stages():
        stage_state = {
            **asdict(stage),
            "attempts": 0,
            "runtime_seconds": 0.0,
            "status": "PENDING",
            "expected_outputs_found": [],
            "expected_outputs_missing": [],
            "quality_gate_results": {},
            "fallback_used": None,
            "warnings": [],
            "next_recommended_action": None,
        }
        for attempt in range(1, stage.max_attempts + 1):
            stage_state["attempts"] = attempt
            log_path = RUN_LOGS_DIR / f"hardening_{stage.stage_id}_attempt_{attempt}.log"
            returncode, runtime = run_command(stage.command, log_path)
            stage_state["runtime_seconds"] += runtime
            found, missing = output_status(stage.expected_outputs)
            stage_state["expected_outputs_found"] = found
            stage_state["expected_outputs_missing"] = missing
            stage_state["quality_gate_results"] = {gate: returncode == 0 and not missing for gate in stage.quality_gates}
            if returncode == 0 and not missing:
                stage_state["status"] = "PASS"
                break
            compile_scripts()
            stage_state["warnings"].append(f"attempt_{attempt}_failed_returncode_{returncode}")
            if attempt == stage.max_attempts:
                for fallback in stage.fallback_strategies:
                    fallback_log = RUN_LOGS_DIR / f"hardening_{stage.stage_id}_fallback.log"
                    fallback_return, fallback_runtime = run_command(fallback, fallback_log)
                    stage_state["runtime_seconds"] += fallback_runtime
                    found, missing = output_status(stage.expected_outputs)
                    stage_state["expected_outputs_found"] = found
                    stage_state["expected_outputs_missing"] = missing
                    stage_state["fallback_used"] = fallback
                    if fallback_return == 0 and not missing:
                        stage_state["status"] = "FALLBACK_PASS"
                        break
                if stage_state["status"] == "PENDING":
                    stage_state["status"] = "BLOCKED" if stage.criticality == "critical" else "DEGRADED"
        state["stages"][stage.stage_id] = stage_state
        save_state(state)
        if stage_state["status"] == "BLOCKED" and stage.criticality == "critical":
            state["next_recommended_action"] = f"Resolve critical stage: {stage.name}"
            save_state(state)
            break

    # Refresh final status now that runner state exists.
    if (PROJECT_ROOT / "src/analysis/build_final_hardening_status.py").exists():
        py = PYTHON
        run_command([py, "src/analysis/build_final_hardening_status.py"], RUN_LOGS_DIR / "hardening_final_status_refresh.log")
    state["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    state["next_recommended_action"] = "Inspect reports/egfr_final_hardening_status.md"
    save_state(state)
    print_final_summary(state)


if __name__ == "__main__":
    main()
