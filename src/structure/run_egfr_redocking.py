"""Run EGFR redocking with AutoDock Vina if receptor/ligand prep succeeds."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import DATA_DIR, FIGURES_DIR, METRICS_DIR, REPORTS_DIR, read_json, save_json, write_text  # noqa: E402


CASE_PATH = DATA_DIR / "processed" / "egfr_redocking_case.json"
DOCKING_DIR = DATA_DIR / "structure_prepared"


def run_vina_binary(case: dict) -> dict:
    """Run local vina binary on prepared PDBQT files."""
    vina = shutil.which("vina") or shutil.which("autodock_vina")
    if vina is None:
        return {"redocking_status": "failed_vina_unavailable", "reason": "No local vina/autodock_vina executable found."}

    protein = PROJECT_ROOT / case.get("protein_pdbqt", "")
    ligand = PROJECT_ROOT / case.get("ligand_pdbqt", "")
    if not protein.exists() or not ligand.exists():
        return {
            "redocking_status": "failed_missing_pdbqt_preparation",
            "reason": "Receptor and/or ligand PDBQT files were not available.",
            "vina_available": True,
        }

    output = DOCKING_DIR / f"{case['pdb_id']}_{case['ligand_id']}_redocked_out.pdbqt"
    log_path = REPORTS_DIR / "run_logs" / "egfr_vina_redocking.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    center = case["box_center"]
    size = case["box_size"]
    command = [
        vina,
        "--receptor",
        str(protein),
        "--ligand",
        str(ligand),
        "--center_x",
        str(center["x"]),
        "--center_y",
        str(center["y"]),
        "--center_z",
        str(center["z"]),
        "--size_x",
        str(size["x"]),
        "--size_y",
        str(size["y"]),
        "--size_z",
        str(size["z"]),
        "--out",
        str(output),
        "--exhaustiveness",
        "8",
    ]
    completed = subprocess.run(command, cwd=PROJECT_ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    if completed.returncode != 0:
        return {
            "redocking_status": "failed_vina_runtime_error",
            "reason": f"Vina exited with return code {completed.returncode}.",
            "vina_available": True,
            "vina_log": str(log_path.relative_to(PROJECT_ROOT)),
        }

    score = None
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        if stripped.startswith("1 "):
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    score = float(parts[1])
                except ValueError:
                    pass
                break
    return {
        "redocking_status": "completed_redocking",
        "reason": "Vina completed successfully.",
        "vina_available": True,
        "docking_score_kcal_mol": score,
        "redocked_pose": str(output.relative_to(PROJECT_ROOT)) if output.exists() else "",
        "pose_recovery_rmsd_angstrom": None,
        "pose_recovery_note": "RMSD was not computed because robust PDBQT-to-reference atom mapping was unavailable in this environment.",
    }


def main() -> None:
    """Run redocking if all docking prerequisites are available."""
    if not CASE_PATH.exists():
        raise FileNotFoundError(f"Missing redocking case: {CASE_PATH}")
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))

    if case.get("case_status") != "prepared":
        result = {
            "redocking_status": "failed_missing_pdbqt_preparation",
            "reason": "Redocking case exists, but receptor/ligand PDBQT preparation did not complete.",
            "vina_available": bool(shutil.which("vina") or shutil.which("autodock_vina")),
            "pdb_id": case.get("pdb_id"),
            "ligand_id": case.get("ligand_id"),
        }
    else:
        result = run_vina_binary(case)
        result["pdb_id"] = case.get("pdb_id")
        result["ligand_id"] = case.get("ligand_id")

    metrics = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
    metrics.update(result)
    if result["redocking_status"] == "completed_redocking":
        metrics["structure_module_status"] = "completed_redocking"
    elif metrics.get("parsed_cocrystal_count", 0) or metrics.get("interaction_residue_count", 0):
        metrics["structure_module_status"] = "structure_analysis_completed_redocking_failed"
    save_json(METRICS_DIR / "egfr_structure_module_metrics.json", metrics)

    report = [
        "# EGFR Redocking Report",
        "",
        f"- Redocking status: {result['redocking_status']}",
        f"- PDB ID: {result.get('pdb_id')}",
        f"- Ligand ID: {result.get('ligand_id')}",
        f"- Vina available: {result.get('vina_available')}",
        f"- Docking score kcal/mol: {result.get('docking_score_kcal_mol')}",
        f"- Pose recovery RMSD angstrom: {result.get('pose_recovery_rmsd_angstrom')}",
        f"- Reason/status note: {result.get('reason')}",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_redocking_report.md", "\n".join(report))

    print(f"Redocking status: {result['redocking_status']}")
    print(f"Vina available: {result.get('vina_available')}")


if __name__ == "__main__":
    main()
