"""Agentic runner for EGFR top-5 structure sanity docking."""

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
STATE_PATH = METRICS_DIR / "agentic_top5_structure_sanity_state.json"
INVENTORY_PATH = METRICS_DIR / "egfr_top5_structure_inventory.json"


@dataclass
class Stage:
    """One structure sanity stage."""

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
    """Create output directories."""
    for path in [REPORTS_DIR, METRICS_DIR, RUN_LOGS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def table_info(path: Path) -> dict:
    """Return shape/columns only for a CSV table."""
    if not path.exists():
        return {"found": False, "path": str(path.relative_to(PROJECT_ROOT))}
    import pandas as pd

    df = pd.read_csv(path)
    return {
        "found": True,
        "path": str(path.relative_to(PROJECT_ROOT)),
        "shape": [int(df.shape[0]), int(df.shape[1])],
        "columns": list(df.columns),
    }


def artifact_inventory() -> dict:
    """Create top-5 structure sanity artifact inventory."""
    ranked = PROJECT_ROOT / "reports" / "egfr_ranked_existing_molecules.csv"
    molecule = PROJECT_ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv"
    medchem = PROJECT_ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv"
    uncertainty = PROJECT_ROOT / "reports" / "egfr_uncertainty_predictions.csv"
    ad = PROJECT_ROOT / "reports" / "egfr_applicability_domain_predictions.csv"
    receptor = PROJECT_ROOT / "data" / "structure_prepared" / "5UG9_receptor.pdbqt"
    ref_ligand = PROJECT_ROOT / "data" / "structure_prepared" / "5UG9_8AM_ligand.pdbqt"
    redocking = PROJECT_ROOT / "reports" / "metrics" / "egfr_redocking_metrics.json"
    case = PROJECT_ROOT / "data" / "processed" / "egfr_redocking_case.json"
    reference_contacts = PROJECT_ROOT / "data" / "processed" / "egfr_binding_site_residues.csv"

    tables = {
        "ranked_existing_molecules": table_info(ranked),
        "model_ready_with_medchem_alerts": table_info(molecule),
        "medchem_alert_table": table_info(medchem),
        "uncertainty_predictions": table_info(uncertainty),
        "applicability_domain_predictions": table_info(ad),
        "reference_contacts": table_info(reference_contacts),
    }
    required_ranked_columns = {
        "molecule_chembl_id",
        "predicted_pIC50",
        "final_score",
        "scaffold_hash",
        "medchem_alert_flag",
    }
    ranked_columns = set(tables["ranked_existing_molecules"].get("columns", []))
    molecule_columns = set(tables["model_ready_with_medchem_alerts"].get("columns", []))
    inventory = {
        "tables": tables,
        "files": {
            "receptor_pdbqt": {"found": receptor.exists(), "path": str(receptor.relative_to(PROJECT_ROOT)), "bytes": receptor.stat().st_size if receptor.exists() else 0},
            "reference_ligand_pdbqt": {"found": ref_ligand.exists(), "path": str(ref_ligand.relative_to(PROJECT_ROOT)), "bytes": ref_ligand.stat().st_size if ref_ligand.exists() else 0},
            "redocking_metrics": {"found": redocking.exists(), "path": str(redocking.relative_to(PROJECT_ROOT))},
            "redocking_case": {"found": case.exists(), "path": str(case.relative_to(PROJECT_ROOT))},
        },
        "quality_gates": {
            "ranked_table_found": ranked.exists(),
            "molecule_id_column_found": "molecule_chembl_id" in ranked_columns,
            "smiles_source_found": bool({"standardized_smiles", "canonical_smiles"} & molecule_columns),
            "receptor_pdbqt_found": receptor.exists() and receptor.stat().st_size > 0,
            "redocking_metrics_found": redocking.exists(),
            "required_selection_fields_found": required_ranked_columns.issubset(ranked_columns),
        },
        "available_selection_fields": sorted(ranked_columns),
    }
    INVENTORY_PATH.write_text(json.dumps(inventory, indent=2, sort_keys=True), encoding="utf-8")
    return inventory


def run_internal_inventory(stage: Stage) -> dict:
    """Run inventory stage and log output."""
    start = time.time()
    stdout_path = RUN_LOGS_DIR / f"{stage.stage_id}.stdout.log"
    stderr_path = RUN_LOGS_DIR / f"{stage.stage_id}.stderr.log"
    try:
        inventory = artifact_inventory()
        missing_gates = [key for key, ok in inventory["quality_gates"].items() if not ok]
        stdout_path.write_text(json.dumps({"missing_gates": missing_gates}, indent=2), encoding="utf-8")
        stderr_path.write_text("", encoding="utf-8")
        status = "PASS" if not missing_gates else "FAIL"
        returncode = 0 if not missing_gates else 1
    except Exception as exc:
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text(str(exc), encoding="utf-8")
        status = "FAIL"
        returncode = 1
    stage.runtime_seconds = time.time() - start
    return {
        "stage_id": stage.stage_id,
        "name": stage.name,
        "command": stage.command,
        "attempts": 1,
        "runtime_seconds": round(stage.runtime_seconds, 3),
        "returncode": returncode,
        "status": status,
        "expected_outputs_found": [str(INVENTORY_PATH.relative_to(PROJECT_ROOT))] if INVENTORY_PATH.exists() else [],
        "expected_outputs_missing": [] if INVENTORY_PATH.exists() else [str(INVENTORY_PATH.relative_to(PROJECT_ROOT))],
        "quality_gates": stage.quality_gates,
        "repair_strategies": stage.repair_strategies,
        "fallback_strategies": stage.fallback_strategies,
        "warnings": stage.warnings,
        "stdout_log": str(stdout_path.relative_to(PROJECT_ROOT)),
        "stderr_log": str(stderr_path.relative_to(PROJECT_ROOT)),
    }


def run_stage(stage: Stage) -> dict:
    """Run one command stage and return a state row."""
    if stage.command and stage.command[0] == "internal:inventory":
        return run_internal_inventory(stage)

    start = time.time()
    stdout_path = RUN_LOGS_DIR / f"{stage.stage_id}.stdout.log"
    stderr_path = RUN_LOGS_DIR / f"{stage.stage_id}.stderr.log"
    last_returncode = None
    for attempt in range(1, stage.max_attempts + 1):
        stage.attempts = attempt
        completed = subprocess.run(
            stage.command,
            cwd=PROJECT_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        last_returncode = completed.returncode
        stdout_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
        stderr_path.write_text(completed.stderr, encoding="utf-8", errors="replace")
        if completed.returncode == 0:
            break

    found = [str(path.relative_to(PROJECT_ROOT)) for path in stage.expected_outputs if path.exists()]
    missing = [str(path.relative_to(PROJECT_ROOT)) for path in stage.expected_outputs if not path.exists()]
    if last_returncode == 0 and not missing:
        stage.status = "PASS"
    elif stage.criticality == "optional":
        stage.status = "DEGRADED"
    else:
        stage.status = "FAIL"
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
    """Run top-5 structure sanity workflow."""
    ensure_dirs()
    python = sys.executable
    stages = [
        Stage(
            stage_id="stage_01_artifact_inventory",
            name="Artifact inventory",
            command=["internal:inventory"],
            expected_outputs=[INVENTORY_PATH],
            quality_gates=["ranked table, molecule table, receptor PDBQT, and redocking metrics found"],
            repair_strategies=["check expected paths and columns"],
            fallback_strategies=["stop if receptor or ranked table is absent"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_02_top5_selection",
            name="Select top-5 clean diverse molecules",
            command=[python, "src/structure/select_top5_for_structure_sanity.py"],
            expected_outputs=[
                REPORTS_DIR / "egfr_top5_structure_selection.csv",
                METRICS_DIR / "egfr_top5_structure_selection_metrics.json",
            ],
            quality_gates=["exactly five selected if possible", "no raw SMILES in public selection table"],
            repair_strategies=["retry after column-name/path check"],
            fallback_strategies=["relax selection while keeping medchem-alert clean"],
            max_attempts=2,
        ),
        Stage(
            stage_id="stage_03_ligand_prep_and_docking",
            name="Prepare ligands and dock into 5UG9 pocket",
            command=[python, "src/structure/dock_top5_ranked_molecules.py"],
            expected_outputs=[
                REPORTS_DIR / "egfr_top5_docking_status.csv",
                METRICS_DIR / "egfr_top5_docking_metrics.json",
            ],
            quality_gates=["each selected molecule has status", "at least three successful dockings if feasible"],
            repair_strategies=["retry after prep-tool/path check"],
            fallback_strategies=["continue with per-molecule degraded status"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_04_report_figures_docs",
            name="Build Vina-score report, figure, and patch docs",
            command=[python, "src/structure/build_top5_structure_sanity_report.py"],
            expected_outputs=[
                REPORTS_DIR / "egfr_top5_structure_sanity_report.md",
                REPORTS_DIR / "egfr_top5_docking_scores.csv",
                REPORTS_DIR / "figures" / "egfr_top5_vina_scores.png",
                METRICS_DIR / "egfr_top5_structure_sanity_metrics.json",
            ],
            quality_gates=["score table exists", "score-only figure exists", "limitations explicit"],
            repair_strategies=["retry after metrics/table path check"],
            fallback_strategies=["write report from docking status only"],
            max_attempts=1,
        ),
        Stage(
            stage_id="stage_05_tests",
            name="Run top-5 structure sanity tests",
            command=[python, "-m", "pytest", "-q", "tests/test_egfr_top5_structure_sanity.py"],
            expected_outputs=[],
            quality_gates=["pytest top-5 tests pass"],
            repair_strategies=["retry after assertion/path wording review"],
            fallback_strategies=["state records failed tests"],
            max_attempts=1,
        ),
    ]

    state_rows: list[dict] = []
    for stage in stages:
        row = run_stage(stage)
        state_rows.append(row)
        if row["status"] == "FAIL" and stage.stage_id in {"stage_01_artifact_inventory", "stage_02_top5_selection"}:
            break

    selection = read_json(METRICS_DIR / "egfr_top5_structure_selection_metrics.json")
    docking = read_json(METRICS_DIR / "egfr_top5_docking_metrics.json")
    sanity = read_json(METRICS_DIR / "egfr_top5_structure_sanity_metrics.json")
    tests_status = next((row["status"] for row in state_rows if row["stage_id"] == "stage_05_tests"), "NOT_RUN")
    successful_docking = int(sanity.get("successful_docking_count", docking.get("successful_docking_count", 0)) or 0)
    stage_failed = any(row["status"] == "FAIL" for row in state_rows)
    if stage_failed:
        final_status = "BLOCKED"
    elif successful_docking >= 3:
        final_status = "DONE"
    else:
        final_status = "DONE_WITH_WARNINGS"
    files_updated = [
        "src/structure/select_top5_for_structure_sanity.py",
        "src/structure/dock_top5_ranked_molecules.py",
        "src/structure/build_top5_structure_sanity_report.py",
        "scripts/agentic_top5_structure_sanity.py",
        "reports/egfr_top5_structure_sanity_report.md",
        "reports/egfr_top5_docking_scores.csv",
        "reports/metrics/egfr_top5_structure_sanity_metrics.json",
        "reports/figures/egfr_top5_vina_scores.png",
        "reports/final_egfr_cadd_qsar_report.md",
        "reports/final_egfr_cv_bullets.md",
        "portfolio_assets/egfr_project_card.md",
        "README.md",
        "tests/test_egfr_top5_structure_sanity.py",
    ]
    state = {
        "FINAL_TOP5_DOCKING_SIMPLIFICATION_STATUS": final_status,
        "stages": state_rows,
        "selected_molecule_count": selection.get("selected_count"),
        "successful_ligand_preparation_count": docking.get("successful_ligand_preparation_count"),
        "successful_docking_count": successful_docking,
        "best_docking_score_kcal_mol": sanity.get("best_docking_score_kcal_mol"),
        "worst_docking_score_kcal_mol": sanity.get("worst_docking_score_kcal_mol"),
        "removed_8am_overlap_columns": sanity.get("removed_8am_overlap_columns"),
        "final_report_patched": sanity.get("final_report_patched"),
        "README_patched": sanity.get("README_patched"),
        "tests_status": tests_status,
        "files_updated": files_updated,
    }
    STATE_PATH.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")

    print(f"FINAL_TOP5_DOCKING_SIMPLIFICATION_STATUS: {final_status}")
    print(f"selected molecule count: {selection.get('selected_count')}")
    print(f"successful ligand preparation count: {docking.get('successful_ligand_preparation_count')}")
    print(f"successful docking count: {successful_docking}")
    print(f"Vina score range: {sanity.get('best_docking_score_kcal_mol')} to {sanity.get('worst_docking_score_kcal_mol')} kcal/mol")
    print(f"removed 8AM-overlap columns: {sanity.get('removed_8am_overlap_columns')}")
    print(f"final report patched: {sanity.get('final_report_patched')}")
    print(f"README patched: {sanity.get('README_patched')}")
    print(f"tests status: {tests_status}")
    print("files updated:")
    for path in files_updated:
        print(f"- {path}")


if __name__ == "__main__":
    main()
