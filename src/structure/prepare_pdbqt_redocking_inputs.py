"""Prepare PDBQT inputs for the selected EGFR redocking case."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CASE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_redocking_case.json"
METRICS_PATH = PROJECT_ROOT / "reports" / "metrics" / "egfr_structure_module_metrics.json"
REDOCKING_METRICS_PATH = PROJECT_ROOT / "reports" / "metrics" / "egfr_redocking_metrics.json"
REPORT_PATH = PROJECT_ROOT / "reports" / "egfr_redocking_report.md"
PREPARED_DIR = PROJECT_ROOT / "data" / "structure_prepared"
RUN_LOG_DIR = PROJECT_ROOT / "reports" / "run_logs"


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


def run_command(command: list[str], log_path: Path) -> tuple[bool, str]:
    """Run a command and write combined output to a log file."""
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


def pdb_coordinates(path: Path) -> list[tuple[float, float, float]]:
    """Read atom coordinates from a PDB-like file."""
    coords: list[tuple[float, float, float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            continue
    return coords


def infer_element(line: str) -> str:
    """Infer element from a PDB atom line."""
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    name = line[12:16].strip()
    letters = "".join(char for char in name if char.isalpha())
    if not letters:
        return "C"
    if len(letters) >= 2 and letters[:2].upper() in {"CL", "BR", "NA", "MG", "ZN", "FE", "MN", "CA"}:
        return letters[:2].upper()
    return letters[0].upper()


def autodock_type(element: str) -> str:
    """Map a PDB element to a conservative AutoDock atom type."""
    element = element.upper()
    if element in {"C", "N", "O", "S", "P", "F", "CL", "BR", "I", "H"}:
        return element
    if element in {"NA", "MG", "ZN", "FE", "MN", "CA"}:
        return element
    return element[:1] or "C"


def write_minimal_receptor_pdbqt(protein_pdb: Path, receptor_pdbqt: Path) -> tuple[bool, str]:
    """Create a minimal receptor PDBQT when dedicated receptor prep tools are unavailable."""
    lines: list[str] = []
    serial = 1
    for line in protein_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        element = infer_element(line)
        if element == "H":
            continue
        atom_type = autodock_type(element)
        base = f"{line[:6]}{serial:5d}{line[11:66]}".ljust(66)
        lines.append(f"{base}    0.000 {atom_type:>2}")
        serial += 1
    if not lines:
        return False, "minimal_receptor_failed_no_atoms"
    receptor_pdbqt.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return receptor_pdbqt.stat().st_size > 0, "minimal_receptor_pdbqt_fallback"


def fetch_ligand_sdf(ligand_id: str, output_path: Path) -> tuple[bool, str]:
    """Download an RCSB ligand SDF for bond-aware ligand preparation."""
    if output_path.exists() and output_path.stat().st_size > 0:
        return True, "existing_ligand_sdf"
    url = f"https://files.rcsb.org/ligands/download/{ligand_id}_ideal.sdf"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            payload = response.read()
    except Exception as exc:
        return False, f"ligand_sdf_download_failed_{exc.__class__.__name__}"
    if len(payload) < 100:
        return False, "ligand_sdf_download_failed_empty_payload"
    output_path.write_bytes(payload)
    return True, "downloaded_rcsb_ligand_sdf"


def prepare_ligand_with_meeko(input_path: Path, ligand_pdbqt: Path) -> tuple[bool, str]:
    """Prepare ligand PDBQT with RDKit and Meeko."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except Exception as exc:
        return False, f"meeko_import_failed_{exc.__class__.__name__}"

    if input_path.suffix.lower() == ".sdf":
        mol = Chem.MolFromMolFile(str(input_path), removeHs=False, sanitize=True)
    else:
        mol = Chem.MolFromPDBFile(str(input_path), removeHs=False, sanitize=False)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                pass
    if mol is None:
        return False, "meeko_failed_rdkit_molecule_parse"
    mol = Chem.AddHs(mol, addCoords=True)
    conformer_status = "existing_conformer"
    if mol.GetNumConformers() == 0:
        status = AllChem.EmbedMolecule(mol, randomSeed=42)
        conformer_status = f"embedded_conformer_status_{status}"
        if status != 0:
            return False, conformer_status
    try:
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return False, "meeko_failed_no_setups"
        written = PDBQTWriterLegacy.write_string(setups[0])
        if isinstance(written, tuple):
            pdbqt_text, ok, error_message = written
            if not ok:
                return False, f"meeko_writer_failed_{error_message or 'unknown'}"
        else:
            pdbqt_text = written
        ligand_pdbqt.write_text(pdbqt_text, encoding="utf-8")
    except Exception as exc:
        return False, f"meeko_failed_{exc.__class__.__name__}"
    if ligand_pdbqt.exists() and ligand_pdbqt.stat().st_size > 0:
        return True, f"prepared_with_meeko_{conformer_status}"
    return False, "meeko_failed_empty_output"


def prepare_with_obabel(input_path: Path, output_path: Path, is_receptor: bool) -> tuple[bool, str]:
    """Prepare PDBQT with Open Babel CLI if available."""
    obabel = tool_path("obabel")
    if not obabel:
        return False, "obabel_unavailable"
    flags = ["-xr"] if is_receptor else ["-xh"]
    commands = [
        [obabel, str(input_path), "-O", str(output_path), *flags],
        [obabel, f"-i{input_path.suffix.lstrip('.')}", str(input_path), "-opdbqt", "-O", str(output_path), *flags],
    ]
    for index, command in enumerate(commands, start=1):
        ok, status = run_command(command, RUN_LOG_DIR / f"obabel_pdbqt_{output_path.stem}_{index}.log")
        if ok and output_path.exists() and output_path.stat().st_size > 0:
            return True, "prepared_with_obabel"
    return False, status


def prepare_with_meeko_script(script_name: str, input_path: Path, output_path: Path) -> tuple[bool, str]:
    """Run a Meeko command-line preparation script if installed."""
    script = tool_path(script_name)
    if not script:
        return False, f"{script_name}_unavailable"
    command = [script, "-i", str(input_path), "-o", str(output_path)]
    ok, status = run_command(command, RUN_LOG_DIR / f"{script_name}_{output_path.stem}.log")
    if ok and output_path.exists() and output_path.stat().st_size > 0:
        return True, f"prepared_with_{script_name}"
    return False, f"{script_name}_{status}"


def prepare_receptor(protein_pdb: Path, receptor_pdbqt: Path) -> tuple[bool, str]:
    """Prepare receptor PDBQT through available tools."""
    for script in ("mk_prepare_receptor.py", "prepare_receptor4.py"):
        ok, status = prepare_with_meeko_script(script, protein_pdb, receptor_pdbqt)
        if ok:
            return ok, status
    ok, status = prepare_with_obabel(protein_pdb, receptor_pdbqt, is_receptor=True)
    if ok:
        return ok, status
    return write_minimal_receptor_pdbqt(protein_pdb, receptor_pdbqt)


def prepare_ligand(ligand_id: str, ligand_pdb: Path, ligand_pdbqt: Path) -> tuple[bool, str, str]:
    """Prepare ligand PDBQT through available tools."""
    ligand_sdf = ligand_pdbqt.with_suffix(".sdf")
    sdf_ok, sdf_status = fetch_ligand_sdf(ligand_id, ligand_sdf)
    ligand_input = ligand_sdf if sdf_ok else ligand_pdb

    for script in ("mk_prepare_ligand.py", "prepare_ligand4.py"):
        ok, status = prepare_with_meeko_script(script, ligand_input, ligand_pdbqt)
        if ok:
            return ok, status, str(ligand_input.relative_to(PROJECT_ROOT))

    ok, status = prepare_with_obabel(ligand_input, ligand_pdbqt, is_receptor=False)
    if ok:
        return ok, status, str(ligand_input.relative_to(PROJECT_ROOT))

    ok, status = prepare_ligand_with_meeko(ligand_input, ligand_pdbqt)
    if ok:
        return ok, status, str(ligand_input.relative_to(PROJECT_ROOT))
    return False, f"{status}; sdf_status={sdf_status}", str(ligand_input.relative_to(PROJECT_ROOT))


def write_report(result: dict) -> None:
    """Write the redocking preparation report section."""
    lines = [
        "# EGFR Redocking Report",
        "",
        "## PDBQT Preparation",
        "",
        f"- Preparation status: {result['pdbqt_preparation_status']}",
        f"- PDB ID: {result.get('pdb_id')}",
        f"- Ligand ID: {result.get('ligand_id')}",
        f"- Receptor PDBQT: `{result.get('receptor_pdbqt') or ''}`",
        f"- Ligand PDBQT: `{result.get('ligand_pdbqt') or ''}`",
        f"- Receptor preparation method: {result.get('receptor_pdbqt_method')}",
        f"- Ligand preparation method: {result.get('ligand_pdbqt_method')}",
        "",
        "Redocking has not been claimed until Vina produces a docking score.",
        "",
    ]
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Prepare receptor and ligand PDBQT files for the current redocking case."""
    if not CASE_PATH.exists():
        raise FileNotFoundError(f"Missing redocking case: {CASE_PATH}")
    case = read_json(CASE_PATH)
    pdb_id = case.get("pdb_id", "5UG9")
    ligand_id = case.get("ligand_id", "8AM")
    protein_pdb = PROJECT_ROOT / case.get("protein_pdb", f"data/structure_prepared/{pdb_id}_protein.pdb")
    ligand_pdb = PROJECT_ROOT / case.get("ligand_pdb", f"data/structure_prepared/{pdb_id}_{ligand_id}_ligand.pdb")
    receptor_pdbqt = PREPARED_DIR / f"{pdb_id}_receptor.pdbqt"
    ligand_pdbqt = PREPARED_DIR / f"{pdb_id}_{ligand_id}_ligand.pdbqt"

    if not protein_pdb.exists() or not ligand_pdb.exists():
        result = {
            "pdbqt_preparation_status": "failed_missing_extracted_structure_files",
            "pdb_id": pdb_id,
            "ligand_id": ligand_id,
            "receptor_pdbqt": "",
            "ligand_pdbqt": "",
        }
    else:
        receptor_ok, receptor_status = prepare_receptor(protein_pdb, receptor_pdbqt)
        ligand_ok, ligand_status, ligand_input = prepare_ligand(ligand_id, ligand_pdb, ligand_pdbqt)
        result = {
            "pdbqt_preparation_status": "prepared" if receptor_ok and ligand_ok else "failed_pdbqt_preparation",
            "pdb_id": pdb_id,
            "ligand_id": ligand_id,
            "receptor_pdbqt": str(receptor_pdbqt.relative_to(PROJECT_ROOT)) if receptor_ok else "",
            "ligand_pdbqt": str(ligand_pdbqt.relative_to(PROJECT_ROOT)) if ligand_ok else "",
            "ligand_preparation_input": ligand_input,
            "receptor_pdbqt_method": receptor_status,
            "ligand_pdbqt_method": ligand_status,
            "receptor_pdbqt_bytes": receptor_pdbqt.stat().st_size if receptor_pdbqt.exists() else 0,
            "ligand_pdbqt_bytes": ligand_pdbqt.stat().st_size if ligand_pdbqt.exists() else 0,
        }
        coords = pdb_coordinates(ligand_pdb)
        if coords and not case.get("box_center"):
            xs, ys, zs = zip(*coords)
            center = {"x": sum(xs) / len(xs), "y": sum(ys) / len(ys), "z": sum(zs) / len(zs)}
            size = {
                "x": max(max(xs) - min(xs) + 12.0, 18.0),
                "y": max(max(ys) - min(ys) + 12.0, 18.0),
                "z": max(max(zs) - min(zs) + 12.0, 18.0),
            }
            case["box_center"] = center
            case["box_size"] = size

        case.update(
            {
                "case_status": "prepared" if receptor_ok and ligand_ok else "prepared_without_pdbqt",
                "protein_pdbqt": result["receptor_pdbqt"],
                "ligand_pdbqt": result["ligand_pdbqt"],
                "protein_conversion_status": result["receptor_pdbqt_method"],
                "ligand_conversion_status": result["ligand_pdbqt_method"],
            }
        )
        write_json(CASE_PATH, case)

    metrics = read_json(METRICS_PATH)
    metrics.update(
        {
            "redocking_case_status": "prepared" if result["pdbqt_preparation_status"] == "prepared" else "prepared_without_pdbqt",
            "redocking_case_pdb_id": result.get("pdb_id"),
            "redocking_case_ligand_id": result.get("ligand_id"),
            "protein_conversion_status": result.get("receptor_pdbqt_method"),
            "ligand_conversion_status": result.get("ligand_pdbqt_method"),
            "docking_ready": result["pdbqt_preparation_status"] == "prepared",
        }
    )
    if result["pdbqt_preparation_status"] != "prepared":
        metrics["structure_module_status"] = "structure_analysis_completed_redocking_failed"
        metrics["redocking_status"] = "failed_missing_pdbqt_preparation"
    write_json(METRICS_PATH, metrics)

    redocking_metrics = read_json(REDOCKING_METRICS_PATH)
    redocking_metrics.update(result)
    write_json(REDOCKING_METRICS_PATH, redocking_metrics)
    write_report(result)

    print(f"PDBQT preparation status: {result['pdbqt_preparation_status']}")
    print(f"Receptor PDBQT path: {result.get('receptor_pdbqt') or 'unavailable'}")
    print(f"Ligand PDBQT path: {result.get('ligand_pdbqt') or 'unavailable'}")


if __name__ == "__main__":
    main()
