"""Run AutoDock Vina redocking for the prepared EGFR co-crystal case."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CASE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_redocking_case.json"
METRICS_PATH = PROJECT_ROOT / "reports" / "metrics" / "egfr_structure_module_metrics.json"
REDOCKING_METRICS_PATH = PROJECT_ROOT / "reports" / "metrics" / "egfr_redocking_metrics.json"
REPORT_PATH = PROJECT_ROOT / "reports" / "egfr_redocking_report.md"
POSE_PATH_TEMPLATE = PROJECT_ROOT / "data" / "structure_prepared" / "{pdb_id}_{ligand_id}_redocked_out.pdbqt"
RUN_LOG_PATH = PROJECT_ROOT / "reports" / "run_logs" / "egfr_vina_redocking.log"


def read_json(path: Path) -> dict:
    """Read JSON or return an empty dictionary if the file is absent."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    """Write indented JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def tool_path(name: str) -> str | None:
    """Find a tool on PATH or in the active Python environment."""
    found = shutil.which(name)
    if found:
        return found
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.exists():
        return str(candidate)
    return None


def parse_vina_score_from_text(text: str) -> float | None:
    """Parse the first Vina score from CLI or pose output text."""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("REMARK VINA RESULT:"):
            parts = stripped.split()
            try:
                return float(parts[3])
            except (IndexError, ValueError):
                continue
        if stripped.startswith("1 "):
            parts = stripped.split()
            try:
                return float(parts[1])
            except (IndexError, ValueError):
                continue
    return None


def run_vina_cli(case: dict, receptor: Path, ligand: Path, pose_path: Path) -> dict:
    """Run a Vina executable if available."""
    vina = tool_path("vina") or tool_path("autodock_vina")
    if not vina:
        return {"ok": False, "status": "vina_cli_unavailable"}
    center = case["box_center"]
    size = case["box_size"]
    command = [
        vina,
        "--receptor",
        str(receptor),
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
        str(pose_path),
        "--exhaustiveness",
        "8",
        "--seed",
        "42",
    ]
    RUN_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    RUN_LOG_PATH.write_text(completed.stdout, encoding="utf-8", errors="replace")
    if completed.returncode != 0:
        return {"ok": False, "status": f"vina_cli_returncode_{completed.returncode}"}
    pose_text = pose_path.read_text(encoding="utf-8", errors="replace") if pose_path.exists() else ""
    return {
        "ok": pose_path.exists() and pose_path.stat().st_size > 0,
        "status": "vina_cli_completed",
        "docking_score_kcal_mol": parse_vina_score_from_text(completed.stdout) or parse_vina_score_from_text(pose_text),
    }


def run_vina_python(case: dict, receptor: Path, ligand: Path, pose_path: Path) -> dict:
    """Run the Python Vina API if available."""
    try:
        from vina import Vina
    except Exception as exc:
        return {"ok": False, "status": f"vina_python_import_failed_{exc.__class__.__name__}"}

    try:
        center = case["box_center"]
        size = case["box_size"]
        vina = Vina(sf_name="vina", seed=42)
        vina.set_receptor(str(receptor))
        vina.set_ligand_from_file(str(ligand))
        vina.compute_vina_maps(
            center=[float(center["x"]), float(center["y"]), float(center["z"])],
            box_size=[float(size["x"]), float(size["y"]), float(size["z"])],
        )
        vina.dock(exhaustiveness=8, n_poses=9)
        vina.write_poses(str(pose_path), n_poses=9, overwrite=True)
    except Exception as exc:
        return {"ok": False, "status": f"vina_python_failed_{exc.__class__.__name__}"}

    pose_text = pose_path.read_text(encoding="utf-8", errors="replace") if pose_path.exists() else ""
    return {
        "ok": pose_path.exists() and pose_path.stat().st_size > 0,
        "status": "vina_python_completed",
        "docking_score_kcal_mol": parse_vina_score_from_text(pose_text),
    }


def write_report(result: dict) -> None:
    """Write or replace the redocking report."""
    lines = [
        "# EGFR Redocking Report",
        "",
        "## Redocking Status",
        "",
        f"- Redocking status: {result['redocking_status']}",
        f"- PDB ID: {result.get('pdb_id')}",
        f"- Ligand ID: {result.get('ligand_id')}",
        f"- Receptor PDBQT: `{result.get('receptor_pdbqt') or ''}`",
        f"- Ligand PDBQT: `{result.get('ligand_pdbqt') or ''}`",
        f"- Vina backend: {result.get('vina_backend')}",
        f"- Docking score kcal/mol: {result.get('docking_score_kcal_mol')}",
        f"- Docked pose: `{result.get('redocked_pose') or ''}`",
        "",
        "Pose-recovery RMSD is evaluated in a separate stage when atom mapping is feasible.",
        "",
    ]
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run Vina redocking if PDBQT inputs are available."""
    if not CASE_PATH.exists():
        raise FileNotFoundError(f"Missing redocking case: {CASE_PATH}")
    case = read_json(CASE_PATH)
    pdb_id = case.get("pdb_id")
    ligand_id = case.get("ligand_id")
    receptor = PROJECT_ROOT / case.get("protein_pdbqt", "")
    ligand = PROJECT_ROOT / case.get("ligand_pdbqt", "")
    pose_path = Path(str(POSE_PATH_TEMPLATE).format(pdb_id=pdb_id, ligand_id=ligand_id))

    receptor_ok = receptor.exists() and receptor.stat().st_size > 0
    ligand_ok = ligand.exists() and ligand.stat().st_size > 0
    if not receptor_ok or not ligand_ok:
        result = {
            "redocking_status": "failed_missing_pdbqt_preparation",
            "pdb_id": pdb_id,
            "ligand_id": ligand_id,
            "receptor_pdbqt": case.get("protein_pdbqt", ""),
            "ligand_pdbqt": case.get("ligand_pdbqt", ""),
            "vina_backend": "not_run",
            "docking_score_kcal_mol": None,
            "redocked_pose": "",
            "reason": "Receptor and ligand PDBQT files must both exist before Vina redocking.",
        }
    else:
        cli_result = run_vina_cli(case, receptor, ligand, pose_path)
        if cli_result["ok"]:
            backend = cli_result["status"]
            docking_score = cli_result.get("docking_score_kcal_mol")
        else:
            python_result = run_vina_python(case, receptor, ligand, pose_path)
            backend = python_result["status"]
            docking_score = python_result.get("docking_score_kcal_mol")

        pose_ok = pose_path.exists() and pose_path.stat().st_size > 0
        score_ok = docking_score is not None
        result = {
            "redocking_status": "completed_redocking_no_rmsd" if pose_ok and score_ok else "failed_vina_runtime_error",
            "pdb_id": pdb_id,
            "ligand_id": ligand_id,
            "receptor_pdbqt": str(receptor.relative_to(PROJECT_ROOT)),
            "ligand_pdbqt": str(ligand.relative_to(PROJECT_ROOT)),
            "vina_backend": backend,
            "docking_score_kcal_mol": docking_score,
            "redocked_pose": str(pose_path.relative_to(PROJECT_ROOT)) if pose_ok else "",
            "reason": "Vina produced a pose and score; RMSD has not yet been evaluated."
            if pose_ok and score_ok
            else "Vina did not produce both a pose and a docking score.",
        }

    redocking_metrics = read_json(REDOCKING_METRICS_PATH)
    redocking_metrics.update(result)
    write_json(REDOCKING_METRICS_PATH, redocking_metrics)

    structure_metrics = read_json(METRICS_PATH)
    structure_metrics.update(result)
    if result["redocking_status"] == "completed_redocking_no_rmsd":
        structure_metrics["structure_module_status"] = "completed_redocking_no_rmsd"
        structure_metrics["docking_ready"] = True
    else:
        structure_metrics["structure_module_status"] = "structure_analysis_completed_redocking_failed"
    write_json(METRICS_PATH, structure_metrics)
    write_report(result)

    print(f"Redocking status: {result['redocking_status']}")
    print(f"Docking score available: {result.get('docking_score_kcal_mol') is not None}")
    print(f"Pose output available: {bool(result.get('redocked_pose'))}")


if __name__ == "__main__":
    main()
