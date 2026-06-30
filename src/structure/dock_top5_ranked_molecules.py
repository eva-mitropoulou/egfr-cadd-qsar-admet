"""Prepare and dock top-5 ranked EGFR molecules into the validated 5UG9 site."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, save_json  # noqa: E402


SELECTION_PATH = REPORTS_DIR / "egfr_top5_structure_selection.csv"
MOLECULE_TABLE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv"
CASE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_redocking_case.json"
REDOCKING_METRICS_PATH = METRICS_DIR / "egfr_redocking_metrics.json"
PREP_DIR = PROJECT_ROOT / "data" / "structure_prepared" / "top5"
DOCKED_DIR = PROJECT_ROOT / "data" / "structure_prepared" / "top5_docked"
RUN_LOG_DIR = REPORTS_DIR / "run_logs"
DOCKING_STATUS_PATH = REPORTS_DIR / "egfr_top5_docking_status.csv"
METRICS_PATH = METRICS_DIR / "egfr_top5_docking_metrics.json"


def read_json(path: Path) -> dict:
    """Read JSON with an empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def safe_name(value: object) -> str:
    """Return a filesystem-safe molecule ID."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))[:80]


def tool_path(name: str) -> str | None:
    """Find an executable on PATH or in the active Python environment."""
    found = shutil.which(name)
    if found:
        return found
    candidate = Path(sys.executable).resolve().parent / name
    if candidate.exists():
        return str(candidate)
    return None


def parse_vina_score(path: Path) -> float | None:
    """Parse the first Vina score from PDBQT output."""
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if stripped.startswith("REMARK VINA RESULT:"):
            parts = stripped.split()
            try:
                return float(parts[3])
            except (IndexError, ValueError):
                continue
    return None


def generate_3d_sdf(smiles: str, output_path: Path) -> tuple[bool, str]:
    """Generate an RDKit 3D conformer and write SDF."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return False, "rdkit_parse_failed"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        return False, f"rdkit_embed_failed_status_{status}"
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=250)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=250)
    except Exception:
        pass
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_path))
    writer.write(mol)
    writer.close()
    if output_path.exists() and output_path.stat().st_size > 0:
        return True, "rdkit_3d_sdf_prepared"
    return False, "rdkit_sdf_empty"


def prepare_ligand_with_meeko(sdf_path: Path, pdbqt_path: Path) -> tuple[bool, str]:
    """Prepare PDBQT with Python Meeko."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except Exception as exc:
        return False, f"meeko_import_failed_{exc.__class__.__name__}"

    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=True)
    if mol is None:
        return False, "meeko_rdkit_sdf_parse_failed"
    try:
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return False, "meeko_no_setups"
        written = PDBQTWriterLegacy.write_string(setups[0])
        if isinstance(written, tuple):
            text, ok, error_message = written
            if not ok:
                return False, f"meeko_writer_failed_{error_message or 'unknown'}"
        else:
            text = written
        pdbqt_path.write_text(text, encoding="utf-8")
    except Exception as exc:
        return False, f"meeko_failed_{exc.__class__.__name__}"
    if pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
        return True, "prepared_with_python_meeko"
    return False, "meeko_empty_output"


def run_command(command: list[str], log_path: Path) -> tuple[bool, str]:
    """Run command and save combined output."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        command,
        cwd=PROJECT_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    log_path.write_text(completed.stdout, encoding="utf-8", errors="replace")
    return completed.returncode == 0, f"returncode_{completed.returncode}"


def prepare_ligand_with_cli(sdf_path: Path, pdbqt_path: Path, molecule_id: str) -> tuple[bool, str]:
    """Try installed ligand preparation CLIs."""
    for tool in ["mk_prepare_ligand.py", "prepare_ligand4.py"]:
        executable = tool_path(tool)
        if not executable:
            continue
        ok, status = run_command(
            [executable, "-i", str(sdf_path), "-o", str(pdbqt_path)],
            RUN_LOG_DIR / f"top5_{molecule_id}_{tool}.log",
        )
        if ok and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            return True, f"prepared_with_{tool}"
        if pdbqt_path.exists() and pdbqt_path.stat().st_size == 0:
            pdbqt_path.unlink(missing_ok=True)
    obabel = tool_path("obabel")
    if obabel:
        ok, status = run_command(
            [obabel, str(sdf_path), "-O", str(pdbqt_path), "-xh"],
            RUN_LOG_DIR / f"top5_{molecule_id}_obabel.log",
        )
        if ok and pdbqt_path.exists() and pdbqt_path.stat().st_size > 0:
            return True, "prepared_with_obabel"
    return False, "no_cli_pdbqt_preparation_succeeded"


def prepare_pdbqt(sdf_path: Path, pdbqt_path: Path, molecule_id: str) -> tuple[bool, str]:
    """Prepare ligand PDBQT through available tools."""
    ok, status = prepare_ligand_with_cli(sdf_path, pdbqt_path, molecule_id)
    if ok:
        return ok, status
    ok, status = prepare_ligand_with_meeko(sdf_path, pdbqt_path)
    if ok:
        return ok, status
    return False, status


def dock_with_vina_python(receptor: Path, ligand: Path, pose_path: Path, case: dict, molecule_id: str) -> tuple[bool, str, float | None]:
    """Run Vina Python API for one ligand."""
    try:
        from vina import Vina
    except Exception as exc:
        return False, f"vina_import_failed_{exc.__class__.__name__}", None
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
        vina.dock(exhaustiveness=int(case.get("exhaustiveness", 8)), n_poses=9)
        pose_path.parent.mkdir(parents=True, exist_ok=True)
        vina.write_poses(str(pose_path), n_poses=9, overwrite=True)
    except Exception as exc:
        log_path = RUN_LOG_DIR / f"top5_{molecule_id}_vina_python_error.log"
        log_path.write_text(str(exc), encoding="utf-8", errors="replace")
        return False, f"vina_python_failed_{exc.__class__.__name__}", None
    score = parse_vina_score(pose_path)
    return pose_path.exists() and pose_path.stat().st_size > 0 and score is not None, "vina_python_completed", score


def dock_with_vina_cli(receptor: Path, ligand: Path, pose_path: Path, case: dict, molecule_id: str) -> tuple[bool, str, float | None]:
    """Run Vina CLI if installed."""
    vina = tool_path("vina") or tool_path("autodock_vina")
    if not vina:
        return False, "vina_cli_unavailable", None
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
        str(case.get("exhaustiveness", 8)),
        "--seed",
        "42",
    ]
    ok, status = run_command(command, RUN_LOG_DIR / f"top5_{molecule_id}_vina_cli.log")
    score = parse_vina_score(pose_path)
    return ok and pose_path.exists() and pose_path.stat().st_size > 0 and score is not None, f"vina_cli_{status}", score


def main() -> None:
    """Prepare and dock selected top-5 molecules."""
    if not SELECTION_PATH.exists():
        raise FileNotFoundError(f"Missing top-5 selection: {SELECTION_PATH}")
    selection = pd.read_csv(SELECTION_PATH)
    molecule_table = pd.read_csv(MOLECULE_TABLE_PATH)
    case = read_json(CASE_PATH)
    redocking = read_json(REDOCKING_METRICS_PATH)
    receptor = PROJECT_ROOT / redocking.get("receptor_pdbqt", case.get("protein_pdbqt", ""))
    if not receptor.exists():
        raise FileNotFoundError(f"Missing 5UG9 receptor PDBQT: {receptor}")

    smiles_column = "standardized_smiles" if "standardized_smiles" in molecule_table.columns else "canonical_smiles"
    molecule_lookup = molecule_table.set_index("molecule_chembl_id")
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    DOCKED_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for _, selected in selection.iterrows():
        molecule_id = str(selected["molecule_id"])
        safe_id = safe_name(molecule_id)
        sdf_path = PREP_DIR / f"{safe_id}.sdf"
        pdbqt_path = PREP_DIR / f"{safe_id}.pdbqt"
        pose_path = DOCKED_DIR / f"{safe_id}_vina_out.pdbqt"
        log_path = RUN_LOG_DIR / f"top5_{safe_id}_docking_status.json"

        row = {
            "molecule_id": molecule_id,
            "rank_before_docking": int(selected["rank_before_docking"]),
            "scaffold_id": str(selected["scaffold_id"]),
            "sdf_path": str(sdf_path.relative_to(PROJECT_ROOT)),
            "ligand_pdbqt_path": str(pdbqt_path.relative_to(PROJECT_ROOT)),
            "pose_path": str(pose_path.relative_to(PROJECT_ROOT)),
            "ligand_preparation_status": "not_started",
            "docking_status": "not_started",
            "vina_score_kcal_mol": None,
            "failure_reason": "",
        }
        if molecule_id not in molecule_lookup.index:
            row["ligand_preparation_status"] = "failed_missing_molecule_record"
            row["docking_status"] = "not_attempted_missing_molecule_record"
            rows.append(row)
            continue
        smiles = molecule_lookup.loc[molecule_id, smiles_column]
        sdf_ok, sdf_status = generate_3d_sdf(str(smiles), sdf_path)
        row["sdf_preparation_status"] = sdf_status
        if not sdf_ok:
            row["ligand_preparation_status"] = "failed_sdf_preparation"
            row["docking_status"] = "not_attempted_ligand_preparation_failed"
            row["failure_reason"] = sdf_status
            rows.append(row)
            continue
        pdbqt_ok, pdbqt_status = prepare_pdbqt(sdf_path, pdbqt_path, safe_id)
        row["ligand_preparation_status"] = pdbqt_status if pdbqt_ok else f"failed_{pdbqt_status}"
        if not pdbqt_ok:
            row["docking_status"] = "not_attempted_pdbqt_preparation_failed"
            row["failure_reason"] = pdbqt_status
            rows.append(row)
            continue

        docking_ok, docking_status, score = dock_with_vina_python(receptor, pdbqt_path, pose_path, case, safe_id)
        if not docking_ok:
            docking_ok, docking_status, score = dock_with_vina_cli(receptor, pdbqt_path, pose_path, case, safe_id)
        row["docking_status"] = docking_status if docking_ok else f"failed_{docking_status}"
        row["vina_score_kcal_mol"] = score
        if not docking_ok:
            row["failure_reason"] = docking_status
        log_path.write_text(json.dumps(row, indent=2, sort_keys=True), encoding="utf-8")
        rows.append(row)

    status = pd.DataFrame(rows)
    status.to_csv(DOCKING_STATUS_PATH, index=False)
    successful_preparation = int(status["ligand_preparation_status"].astype(str).str.contains("prepared_with").sum())
    successful_docking = int(status["docking_status"].astype(str).str.contains("completed").sum())
    scores = pd.to_numeric(status["vina_score_kcal_mol"], errors="coerce").dropna()
    metrics = {
        "selected_molecule_count": int(len(selection)),
        "successful_ligand_preparation_count": successful_preparation,
        "successful_docking_count": successful_docking,
        "pdb_id": case.get("pdb_id", "5UG9"),
        "reference_ligand_id": case.get("ligand_id", "8AM"),
        "receptor_pdbqt": str(receptor.relative_to(PROJECT_ROOT)),
        "docking_box_center": case.get("box_center"),
        "docking_box_size": case.get("box_size"),
        "docking_status_table": str(DOCKING_STATUS_PATH.relative_to(PROJECT_ROOT)),
        "best_docking_score_kcal_mol": float(scores.min()) if not scores.empty else None,
        "worst_docking_score_kcal_mol": float(scores.max()) if not scores.empty else None,
        "docking_backend_preference": "vina_python_then_cli",
    }
    save_json(METRICS_PATH, metrics)

    print(f"Docking status table: {DOCKING_STATUS_PATH}")
    print(f"Successful ligand preparation count: {successful_preparation}")
    print(f"Successful docking count: {successful_docking}")


if __name__ == "__main__":
    main()
