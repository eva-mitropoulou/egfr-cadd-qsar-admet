"""Run the EGFR medicinal-chemistry alert evidence layer."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPORTS_DIR = PROJECT_ROOT / "reports"
METRICS_DIR = REPORTS_DIR / "metrics"
RUN_LOGS_DIR = REPORTS_DIR / "run_logs"
STATE_PATH = METRICS_DIR / "agentic_medchem_alert_state.json"


@dataclass
class Stage:
    """One runner stage."""

    stage_id: str
    name: str
    command: list[str]
    expected_outputs: list[Path]
    quality_gates: list[str]
    repair_strategies: list[str]
    fallback_strategies: list[str]
    max_attempts: int = 1
    criticality: str = "required"
    status: str = "PENDING"
    attempts: int = 0
    runtime_seconds: float = 0.0
    warnings: list[str] = field(default_factory=list)


def ensure_dirs() -> None:
    """Create runner output directories."""
    for path in [REPORTS_DIR, METRICS_DIR, RUN_LOGS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def run_stage(stage: Stage) -> dict:
    """Run a stage and return state details."""
    start = time.time()
    last_returncode = None
    stdout_path = RUN_LOGS_DIR / f"{stage.stage_id}.stdout.log"
    stderr_path = RUN_LOGS_DIR / f"{stage.stage_id}.stderr.log"

    for attempt in range(1, stage.max_attempts + 1):
        stage.attempts = attempt
        process = subprocess.run(
            stage.command,
            cwd=PROJECT_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        last_returncode = process.returncode
        stdout_path.write_text(process.stdout, encoding="utf-8")
        stderr_path.write_text(process.stderr, encoding="utf-8")
        if process.returncode == 0:
            break

    missing = [str(path.relative_to(PROJECT_ROOT)) for path in stage.expected_outputs if not path.exists()]
    found = [str(path.relative_to(PROJECT_ROOT)) for path in stage.expected_outputs if path.exists()]
    if last_returncode == 0 and not missing:
        stage.status = "PASS"
    elif stage.criticality == "required":
        stage.status = "FAIL"
    else:
        stage.status = "DEGRADED"

    stage.runtime_seconds = time.time() - start
    return {
        "stage_id": stage.stage_id,
        "name": stage.name,
        "command": stage.command,
        "attempts": stage.attempts,
        "runtime_seconds": round(stage.runtime_seconds, 3),
        "returncode": last_returncode,
        "status": stage.status,
        "expected_outputs_found": found,
        "expected_outputs_missing": missing,
        "quality_gates": stage.quality_gates,
        "repair_strategies": stage.repair_strategies,
        "fallback_strategies": stage.fallback_strategies,
        "warnings": stage.warnings,
        "stdout_log": str(stdout_path.relative_to(PROJECT_ROOT)),
        "stderr_log": str(stderr_path.relative_to(PROJECT_ROOT)),
    }


def read_json(path: Path) -> dict:
    """Read JSON with empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    """Run medchem-alert annotation, triage, sensitivity, report patching, and tests."""
    ensure_dirs()
    python = sys.executable
    stages = [
        Stage(
            stage_id="stage_01_medchem_alert_annotation",
            name="Annotate PAINS/Brenk/external SMARTS alerts",
            command=[python, "src/chem/medchem_alerts.py"],
            expected_outputs=[
                PROJECT_ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv",
                REPORTS_DIR / "egfr_medchem_alerts_report.md",
                METRICS_DIR / "egfr_medchem_alerts_metrics.json",
            ],
            quality_gates=[
                "row count unchanged from model-ready input",
                "PAINS, Brenk, unwanted, and combined alert columns exist",
                "external CSV absence is degraded-but-valid",
            ],
            repair_strategies=["retry after path/import check"],
            fallback_strategies=["continue with RDKit PAINS/Brenk/NIH only if external CSV is absent"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_02_triage_refresh",
            name="Refresh alert-aware candidate triage",
            command=[python, "src/triage/rank_egfr_existing_molecules.py"],
            expected_outputs=[
                REPORTS_DIR / "egfr_candidate_triage_report.md",
                REPORTS_DIR / "egfr_ranked_existing_molecules.csv",
                METRICS_DIR / "egfr_candidate_triage_metrics.json",
            ],
            quality_gates=[
                "ranked table keeps all existing molecules",
                "ranked table contains alert-risk columns",
                "top-20 alert composition is reported",
            ],
            repair_strategies=["retry after column-name check"],
            fallback_strategies=["preserve previous ranked table and mark triage degraded"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_03_sensitivity_analysis",
            name="Run alert-exclusion sensitivity analysis",
            command=[python, "src/analysis/egfr_medchem_alert_sensitivity.py"],
            expected_outputs=[
                REPORTS_DIR / "egfr_medchem_alert_sensitivity_report.md",
                METRICS_DIR / "egfr_medchem_alert_sensitivity_metrics.json",
                REPORTS_DIR / "figures" / "medchem_alert_activity_distribution.png",
                REPORTS_DIR / "figures" / "medchem_alert_subset_model_performance.png",
                REPORTS_DIR / "figures" / "top_ranked_medchem_alert_composition.png",
            ],
            quality_gates=[
                "all four subsets attempted or explicitly marked unavailable",
                "random and scaffold metrics reported for full vs PAINS-excluded",
                "main benchmark files are not overwritten",
            ],
            repair_strategies=["retry after path/import check"],
            fallback_strategies=["report composition-only sensitivity if RF retraining fails"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_04_report_patch",
            name="Patch final report, README, and project card",
            command=[python, "src/analysis/patch_medchem_alert_reports.py"],
            expected_outputs=[
                REPORTS_DIR / "final_egfr_cadd_qsar_report.md",
                REPORTS_DIR / "final_egfr_cv_bullets.md",
                PROJECT_ROOT / "portfolio_assets" / "egfr_project_card.md",
                PROJECT_ROOT / "README.md",
            ],
            quality_gates=[
                "required wording about annotations and sensitivity filters present",
                "counts and fractions included",
                "no overclaiming language added",
            ],
            repair_strategies=["retry after metrics-file check"],
            fallback_strategies=["leave standalone medchem reports if patching fails"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_05_medchem_tests",
            name="Run medchem-alert tests",
            command=[python, "-m", "pytest", "-q", "tests/test_egfr_medchem_alerts.py"],
            expected_outputs=[],
            quality_gates=["pytest medchem-alert tests pass"],
            repair_strategies=["retry after test assertion review"],
            fallback_strategies=["write state with tests failed"],
            max_attempts=1,
        ),
    ]

    state_rows: list[dict] = []
    for stage in stages:
        row = run_stage(stage)
        state_rows.append(row)
        if row["status"] == "FAIL" and stage.criticality == "required":
            break

    alerts = read_json(METRICS_DIR / "egfr_medchem_alerts_metrics.json")
    sensitivity = read_json(METRICS_DIR / "egfr_medchem_alert_sensitivity_metrics.json")
    triage = read_json(METRICS_DIR / "egfr_candidate_triage_metrics.json")
    tests_status = next((row["status"] for row in state_rows if row["stage_id"] == "stage_05_medchem_tests"), "NOT_RUN")
    any_fail = any(row["status"] == "FAIL" for row in state_rows)
    final_status = "DONE" if not any_fail else "BLOCKED"
    if sensitivity.get("subset_availability") and any(
        value != "evaluated" for value in sensitivity.get("subset_availability", {}).values()
    ):
        final_status = "DONE_WITH_WARNINGS" if final_status == "DONE" else final_status

    output_files = [
        "src/chem/medchem_alerts.py",
        "src/triage/admet_risk_scoring.py",
        "src/triage/rank_egfr_existing_molecules.py",
        "src/analysis/egfr_medchem_alert_sensitivity.py",
        "src/analysis/patch_medchem_alert_reports.py",
        "reports/egfr_medchem_alerts_report.md",
        "reports/egfr_medchem_alert_sensitivity_report.md",
        "reports/metrics/egfr_medchem_alerts_metrics.json",
        "reports/metrics/egfr_medchem_alert_sensitivity_metrics.json",
        "reports/egfr_ranked_existing_molecules.csv",
        "reports/final_egfr_cadd_qsar_report.md",
        "reports/final_egfr_cv_bullets.md",
        "portfolio_assets/egfr_project_card.md",
        "README.md",
        "tests/test_egfr_medchem_alerts.py",
    ]
    state = {
        "FINAL_MEDCHEM_ALERT_STATUS": final_status,
        "stages": state_rows,
        "next_recommended_action": "Review reports/egfr_medchem_alerts_report.md and reports/egfr_medchem_alert_sensitivity_report.md.",
        "files_updated": output_files,
    }
    STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")

    input_rows = alerts.get("input_row_count", "unavailable")
    output_rows = alerts.get("output_row_count", "unavailable")
    top20_clean = triage.get("top20_clean_medchem_alert_count", "unavailable")
    print(f"FINAL_MEDCHEM_ALERT_STATUS: {final_status}")
    print(f"input model-ready row count: {input_rows}")
    print(f"annotated row count: {output_rows}")
    print(f"PAINS-flagged count/fraction: {alerts.get('pains_flagged_count', 'unavailable')}/{alerts.get('pains_flagged_fraction', 'unavailable')}")
    print(f"Brenk-flagged count/fraction: {alerts.get('brenk_flagged_count', 'unavailable')}/{alerts.get('brenk_flagged_fraction', 'unavailable')}")
    print(
        "unwanted-substructure-flagged count/fraction: "
        f"{alerts.get('unwanted_substructure_flagged_count', 'unavailable')}/"
        f"{alerts.get('unwanted_substructure_flagged_fraction', 'unavailable')}"
    )
    print(f"combined medchem-alert count/fraction: {alerts.get('combined_medchem_alert_count', 'unavailable')}/{alerts.get('combined_medchem_alert_fraction', 'unavailable')}")
    print(f"external CSV status: {alerts.get('external_csv_found', 'unavailable')}")
    print(f"sensitivity analysis status: {'available' if sensitivity else 'unavailable'}")
    print(f"main benchmark preserved yes/no: {sensitivity.get('primary_benchmark_preserved', 'unavailable')}")
    print(f"top-20 clean count: {top20_clean}")
    print(f"tests status: {tests_status}")
    print("files updated:")
    for path in output_files:
        print(f"- {path}")


if __name__ == "__main__":
    main()
