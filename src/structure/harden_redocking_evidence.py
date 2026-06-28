"""Harden and audit the completed EGFR Vina redocking evidence."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import REPORTS_DIR, METRICS_DIR, save_json, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


REDOCKING_METRICS_PATH = METRICS_DIR / "egfr_redocking_metrics.json"
STRUCTURE_METRICS_PATH = METRICS_DIR / "egfr_structure_module_metrics.json"
CASE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_redocking_case.json"
AUDIT_METRICS_PATH = METRICS_DIR / "egfr_redocking_audit_metrics.json"
AUDIT_REPORT_PATH = REPORTS_DIR / "egfr_redocking_audit_report.md"
REDOCKING_REPORT_PATH = REPORTS_DIR / "egfr_redocking_report.md"
OVERLAY_PATH = REPORTS_DIR / "figures" / "5UG9_8AM_redocking_pose_overlay.png"
OVERLAY_SCRIPT_PATH = REPORTS_DIR / "structure_visualization" / "5UG9_8AM_overlay.pml"


def read_json(path: Path) -> dict:
    """Read JSON with empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def command_version(command: list[str]) -> str:
    """Return a concise tool version string or unavailability reason."""
    executable = shutil.which(command[0]) or str(Path(sys.executable).resolve().parent / command[0])
    if not Path(executable).exists() and shutil.which(command[0]) is None:
        return "unavailable"
    try:
        completed = subprocess.run(
            [executable, *command[1:]],
            cwd=PROJECT_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=20,
        )
    except Exception as exc:
        return f"unavailable_{exc.__class__.__name__}"
    first_line = completed.stdout.strip().splitlines()[0] if completed.stdout.strip() else f"returncode_{completed.returncode}"
    return first_line[:160]


def package_version(package_name: str) -> str:
    """Return installed package version if available."""
    try:
        import importlib.metadata as metadata

        return metadata.version(package_name)
    except Exception:
        return "unavailable"


def parse_pdbqt_coords(path: Path, first_pose_only: bool) -> np.ndarray:
    """Parse heavy-atom PDBQT coordinates."""
    coords: list[list[float]] = []
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            if saw_model:
                break
            saw_model = True
            continue
        if first_pose_only and line.startswith("ENDMDL") and saw_model:
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_type = line.split()[-1].upper() if line.split() else ""
        if atom_type in {"H", "HD", "HS"}:
            continue
        try:
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            continue
    return np.asarray(coords, dtype=float)


def create_overlay(reference_path: Path, docked_path: Path) -> str:
    """Create a simple coordinate overlay figure with matplotlib."""
    reference = parse_pdbqt_coords(reference_path, first_pose_only=False)
    docked = parse_pdbqt_coords(docked_path, first_pose_only=True)
    if reference.size == 0 or docked.size == 0:
        write_overlay_script(reference_path, docked_path)
        return "overlay_script_created_no_coordinates"

    OVERLAY_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(6, 5))
    axis = fig.add_subplot(111, projection="3d")
    axis.scatter(reference[:, 0], reference[:, 1], reference[:, 2], s=20, label="co-crystal ligand", alpha=0.85)
    axis.scatter(docked[:, 0], docked[:, 1], docked[:, 2], s=20, label="redocked pose", alpha=0.85)
    axis.set_title("5UG9 with ligand 8AM Redocking Pose Overlay")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    axis.legend(loc="best")
    fig.tight_layout()
    fig.savefig(OVERLAY_PATH, dpi=220)
    plt.close(fig)
    write_overlay_script(reference_path, docked_path)
    return "overlay_figure_created"


def write_overlay_script(reference_path: Path, docked_path: Path) -> None:
    """Write a PyMOL overlay helper script."""
    OVERLAY_SCRIPT_PATH.parent.mkdir(parents=True, exist_ok=True)
    script = [
        "reinitialize",
        f"load {reference_path.resolve()}, reference_ligand",
        f"load {docked_path.resolve()}, redocked_pose",
        "show sticks, reference_ligand",
        "show sticks, redocked_pose",
        "color cyan, reference_ligand",
        "color magenta, redocked_pose",
        "set ray_opaque_background, off",
        "zoom",
    ]
    OVERLAY_SCRIPT_PATH.write_text("\n".join(script) + "\n", encoding="utf-8")


def main() -> None:
    """Write redocking audit metrics/report and overlay artifact."""
    redocking = read_json(REDOCKING_METRICS_PATH)
    structure = read_json(STRUCTURE_METRICS_PATH)
    case = read_json(CASE_PATH)
    receptor = PROJECT_ROOT / redocking.get("receptor_pdbqt", "")
    ligand = PROJECT_ROOT / redocking.get("ligand_pdbqt", "")
    pose = PROJECT_ROOT / redocking.get("redocked_pose", "")
    overlay_status = create_overlay(ligand, pose) if ligand.exists() and pose.exists() else "overlay_not_created_missing_files"

    from rdkit import rdBase

    audit = {
        "status": "completed" if redocking.get("redocking_status") == "completed_redocking" else "degraded_redocking_not_completed",
        "pdb_id": redocking.get("pdb_id"),
        "ligand_id": redocking.get("ligand_id"),
        "redocking_status": redocking.get("redocking_status"),
        "docking_score_kcal_mol": redocking.get("docking_score_kcal_mol"),
        "pose_recovery_rmsd_angstrom": redocking.get("pose_recovery_rmsd_angstrom"),
        "pose_recovery_status": redocking.get("pose_recovery_status"),
        "vina_backend": redocking.get("vina_backend"),
        "vina_cli_version": command_version(["vina", "--version"]),
        "vina_python_version": package_version("vina"),
        "meeko_version": package_version("meeko"),
        "openbabel_cli_version": command_version(["obabel", "-V"]),
        "openbabel_wheel_version": package_version("openbabel-wheel"),
        "rdkit_version": rdBase.rdkitVersion,
        "receptor_pdbqt": redocking.get("receptor_pdbqt"),
        "ligand_pdbqt": redocking.get("ligand_pdbqt"),
        "redocked_pose": redocking.get("redocked_pose"),
        "co_crystal_ligand_reference_path": case.get("ligand_pdb"),
        "receptor_pdbqt_bytes": receptor.stat().st_size if receptor.exists() else 0,
        "ligand_pdbqt_bytes": ligand.stat().st_size if ligand.exists() else 0,
        "redocked_pose_bytes": pose.stat().st_size if pose.exists() else 0,
        "docking_box_center": case.get("box_center"),
        "docking_box_size": case.get("box_size"),
        "exhaustiveness": 8,
        "number_of_poses_requested": 9,
        "rmsd_method": "direct fixed-frame coordinate RMSD from prepared ligand PDBQT original coordinates to first Vina pose",
        "rmsd_atom_scope": "heavy atom",
        "hydrogens_ignored_for_rmsd": True,
        "atom_mapping": "direct prepared-ligand PDBQT atom order",
        "overlay_artifact_status": overlay_status,
        "overlay_figure": str(OVERLAY_PATH.relative_to(PROJECT_ROOT)) if OVERLAY_PATH.exists() else "",
        "overlay_script": str(OVERLAY_SCRIPT_PATH.relative_to(PROJECT_ROOT)) if OVERLAY_SCRIPT_PATH.exists() else "",
    }
    save_json(AUDIT_METRICS_PATH, audit)
    structure.update(
        {
            "redocking_audit_status": audit["status"],
            "redocking_overlay_artifact_status": overlay_status,
            "redocking_audit_metrics": str(AUDIT_METRICS_PATH.relative_to(PROJECT_ROOT)),
        }
    )
    save_json(STRUCTURE_METRICS_PATH, structure)

    lines = [
        "# EGFR Redocking Audit Report",
        "",
        "This audit hardens the completed retrospective Vina redocking evidence for the 5UG9 with ligand 8AM co-crystal case.",
        "",
        "## Result",
        "",
        f"- Redocking status: {audit['redocking_status']}",
        f"- Docking score: {audit['docking_score_kcal_mol']} kcal/mol",
        f"- Pose recovery RMSD: {audit['pose_recovery_rmsd_angstrom']} angstrom",
        f"- Pose recovery status: {audit['pose_recovery_status']}",
        "",
        "## Tooling",
        "",
        f"- Vina Python package: {audit['vina_python_version']}",
        f"- Vina CLI: {audit['vina_cli_version']}",
        f"- Meeko: {audit['meeko_version']}",
        f"- OpenBabel CLI: {audit['openbabel_cli_version']}",
        f"- RDKit: {audit['rdkit_version']}",
        "",
        "## Input/Output Audit",
        "",
        f"- Receptor PDBQT: `{audit['receptor_pdbqt']}` ({audit['receptor_pdbqt_bytes']} bytes)",
        f"- Ligand PDBQT: `{audit['ligand_pdbqt']}` ({audit['ligand_pdbqt_bytes']} bytes)",
        f"- Docked pose: `{audit['redocked_pose']}` ({audit['redocked_pose_bytes']} bytes)",
        f"- Overlay artifact status: {audit['overlay_artifact_status']}",
        f"- Overlay figure: `{audit['overlay_figure']}`",
        f"- Overlay script: `{audit['overlay_script']}`",
        "",
        "## RMSD Policy",
        "",
        f"- Method: {audit['rmsd_method']}",
        f"- Atom scope: {audit['rmsd_atom_scope']}",
        f"- Hydrogens ignored: {audit['hydrogens_ignored_for_rmsd']}",
        f"- Atom mapping: {audit['atom_mapping']}",
        "",
        "This is retrospective redocking validation of a known co-crystallized ligand and pose-recovery audit.",
        "",
    ]
    write_text(AUDIT_REPORT_PATH, "\n".join(lines))

    redocking_report = [
        "# EGFR Redocking Report",
        "",
        "## Redocking Result",
        "",
        f"- Redocking status: {audit['redocking_status']}",
        f"- PDB ID: {audit['pdb_id']}",
        f"- Ligand ID: {audit['ligand_id']}",
        f"- Receptor PDBQT: `{audit['receptor_pdbqt']}`",
        f"- Ligand PDBQT: `{audit['ligand_pdbqt']}`",
        f"- Vina backend: {audit['vina_backend']}",
        f"- Docking score kcal/mol: {audit['docking_score_kcal_mol']}",
        f"- Docked pose: `{audit['redocked_pose']}`",
        f"- Pose recovery RMSD angstrom: {audit['pose_recovery_rmsd_angstrom']}",
        f"- Pose recovery status: {audit['pose_recovery_status']}",
        f"- Audit report: `{AUDIT_REPORT_PATH.relative_to(PROJECT_ROOT)}`",
        "",
        "This is retrospective redocking of a known co-crystallized ligand for pose-recovery review.",
        "",
    ]
    write_text(REDOCKING_REPORT_PATH, "\n".join(redocking_report))

    print(f"Redocking audit status: {audit['status']}")
    print(f"Overlay artifact status: {overlay_status}")


if __name__ == "__main__":
    main()
