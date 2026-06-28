"""Prepare an EGFR redocking case from extracted co-crystal protein/ligand files."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import DATA_DIR, METRICS_DIR, read_json, save_json, write_text  # noqa: E402


STRUCTURE_TABLE = DATA_DIR / "processed" / "egfr_structure_candidates.csv"
CASE_PATH = DATA_DIR / "processed" / "egfr_redocking_case.json"


def pdb_coords(path: Path) -> np.ndarray:
    """Read coordinates from PDB atom records."""
    coords = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except Exception:
                continue
    return np.asarray(coords, dtype=float)


def try_obabel(input_path: Path, output_path: Path) -> tuple[bool, str]:
    """Try converting a PDB file to PDBQT with Open Babel if available."""
    obabel = shutil.which("obabel")
    if obabel is None:
        return False, "obabel_unavailable"
    completed = subprocess.run(
        [obabel, str(input_path), "-O", str(output_path)],
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if completed.returncode == 0 and output_path.exists() and output_path.stat().st_size > 0:
        return True, "converted_with_obabel"
    return False, f"obabel_failed_returncode_{completed.returncode}"


def main() -> None:
    """Prepare redocking box and optional PDBQT files."""
    if not STRUCTURE_TABLE.exists():
        raise FileNotFoundError(f"Missing structure table: {STRUCTURE_TABLE}")

    structures = pd.read_csv(STRUCTURE_TABLE)
    usable = structures[structures["parse_status"] == "parsed_with_ligand"].copy()
    if usable.empty:
        payload = {
            "case_status": "failed_no_parsed_cocrystal",
            "reason": "No parsed EGFR co-crystal with selected ligand was available.",
        }
        save_json(CASE_PATH, payload)
        print("Redocking case status: failed_no_parsed_cocrystal")
        return

    selected = usable.sort_values(["resolution_angstrom", "ligand_heavy_atom_count"], ascending=[True, False]).iloc[0]
    protein_pdb = PROJECT_ROOT / selected["protein_pdb"]
    ligand_pdb = PROJECT_ROOT / selected["ligand_pdb"]
    ligand_xyz = pdb_coords(ligand_pdb)
    if ligand_xyz.size == 0:
        payload = {
            "case_status": "failed_no_ligand_coordinates",
            "pdb_id": selected["pdb_id"],
            "ligand_id": selected["ligand_id"],
        }
        save_json(CASE_PATH, payload)
        print("Redocking case status: failed_no_ligand_coordinates")
        return

    center = ligand_xyz.mean(axis=0)
    span = ligand_xyz.max(axis=0) - ligand_xyz.min(axis=0)
    size = np.maximum(span + 12.0, np.array([18.0, 18.0, 18.0]))

    protein_pdbqt = protein_pdb.with_suffix(".pdbqt")
    ligand_pdbqt = ligand_pdb.with_suffix(".pdbqt")
    protein_converted, protein_conversion_status = try_obabel(protein_pdb, protein_pdbqt)
    ligand_converted, ligand_conversion_status = try_obabel(ligand_pdb, ligand_pdbqt)

    payload = {
        "case_status": "prepared" if protein_converted and ligand_converted else "prepared_without_pdbqt",
        "pdb_id": selected["pdb_id"],
        "ligand_id": selected["ligand_id"],
        "resolution_angstrom": None if pd.isna(selected["resolution_angstrom"]) else float(selected["resolution_angstrom"]),
        "protein_pdb": selected["protein_pdb"],
        "ligand_pdb": selected["ligand_pdb"],
        "protein_pdbqt": str(protein_pdbqt.relative_to(PROJECT_ROOT)) if protein_converted else "",
        "ligand_pdbqt": str(ligand_pdbqt.relative_to(PROJECT_ROOT)) if ligand_converted else "",
        "protein_conversion_status": protein_conversion_status,
        "ligand_conversion_status": ligand_conversion_status,
        "box_center": {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])},
        "box_size": {"x": float(size[0]), "y": float(size[1]), "z": float(size[2])},
    }
    save_json(CASE_PATH, payload)

    metrics = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
    metrics.update(
        {
            "redocking_case_status": payload["case_status"],
            "redocking_case_pdb_id": payload["pdb_id"],
            "redocking_case_ligand_id": payload["ligand_id"],
            "protein_conversion_status": protein_conversion_status,
            "ligand_conversion_status": ligand_conversion_status,
        }
    )
    save_json(METRICS_DIR / "egfr_structure_module_metrics.json", metrics)

    report = [
        "# EGFR Redocking Case Preparation",
        "",
        f"- Case status: {payload['case_status']}",
        f"- PDB ID: {payload['pdb_id']}",
        f"- Ligand ID: {payload['ligand_id']}",
        f"- Protein conversion status: {protein_conversion_status}",
        f"- Ligand conversion status: {ligand_conversion_status}",
        f"- Case file: `{CASE_PATH.relative_to(PROJECT_ROOT)}`",
        "",
    ]
    write_text(PROJECT_ROOT / "reports" / "egfr_redocking_case_report.md", "\n".join(report))

    print(f"Redocking case status: {payload['case_status']}")
    print(f"Redocking case PDB: {payload['pdb_id']}")


if __name__ == "__main__":
    main()
