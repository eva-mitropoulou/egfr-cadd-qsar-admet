"""Retrieve and parse real EGFR co-crystal structures for structure-based analysis."""

from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import DATA_DIR, METRICS_DIR, REPORTS_DIR, save_json, write_text  # noqa: E402


STRUCTURE_DIR = DATA_DIR / "structures"
PREP_DIR = DATA_DIR / "structure_prepared"
STRUCTURE_TABLE = DATA_DIR / "processed" / "egfr_structure_candidates.csv"
BINDING_SITE_TABLE = DATA_DIR / "processed" / "egfr_binding_site_residues.csv"


PDB_CANDIDATES = [
    {"pdb_id": "1M17", "description": "EGFR kinase domain with co-crystallized inhibitor"},
    {"pdb_id": "2ITY", "description": "EGFR kinase domain with co-crystallized inhibitor"},
    {"pdb_id": "4HJO", "description": "EGFR kinase domain with co-crystallized inhibitor"},
    {"pdb_id": "5UG9", "description": "EGFR kinase domain with co-crystallized inhibitor"},
]

EXCLUDED_HET = {
    "HOH",
    "WAT",
    "DOD",
    "NA",
    "K",
    "CL",
    "MG",
    "MN",
    "ZN",
    "CA",
    "SO4",
    "PO4",
    "GOL",
    "EDO",
    "PEG",
    "ACT",
    "DMS",
    "TRS",
    "TPO",
    "PTR",
    "SEP",
    "MSE",
}

HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}
CHARGED_RESIDUES = {"ASP", "GLU", "LYS", "ARG", "HIS"}


def fetch_pdb(pdb_id: str, output_path: Path) -> str:
    """Fetch a PDB file from RCSB if not already present."""
    if output_path.exists() and output_path.stat().st_size > 0:
        return "local_file_present"
    try:
        import requests

        response = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=20)
        if response.ok and ("ATOM" in response.text or "HETATM" in response.text):
            output_path.write_text(response.text, encoding="utf-8")
            return "fetched"
        return f"fetch_failed_http_{response.status_code}"
    except Exception as exc:  # network can be unavailable in sandboxed runs
        return f"fetch_failed_{type(exc).__name__}"


def parse_resolution(lines: list[str]) -> float | None:
    """Parse PDB resolution from REMARK 2 where available."""
    for line in lines:
        if line.startswith("REMARK   2 RESOLUTION."):
            match = re.search(r"([0-9]+\.[0-9]+)\s+ANGSTROMS", line)
            if match:
                return float(match.group(1))
    return None


def atom_record(line: str) -> dict[str, object] | None:
    """Parse coordinates from ATOM/HETATM PDB records."""
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None
    try:
        return {
            "record": line[:6].strip(),
            "atom_name": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21].strip(),
            "resseq": line[22:26].strip(),
            "icode": line[26].strip(),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
            "line": line.rstrip("\n"),
        }
    except Exception:
        return None


def parse_atoms(path: Path) -> tuple[list[str], pd.DataFrame]:
    """Parse PDB lines and atom coordinate table."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = [record for line in lines if (record := atom_record(line)) is not None]
    return lines, pd.DataFrame(atoms)


def choose_ligand(atoms: pd.DataFrame) -> dict[str, object] | None:
    """Choose the largest non-solvent HET ligand."""
    het = atoms[(atoms["record"] == "HETATM") & (~atoms["resname"].isin(EXCLUDED_HET))].copy()
    if het.empty:
        return None
    grouped = (
        het.groupby(["resname", "chain", "resseq", "icode"], dropna=False)
        .agg(atom_count=("atom_name", "size"), heavy_atom_count=("element", lambda values: int((values != "H").sum())))
        .reset_index()
        .sort_values(["heavy_atom_count", "atom_count"], ascending=False)
    )
    best = grouped.iloc[0]
    return {
        "ligand_id": str(best["resname"]),
        "ligand_chain": str(best["chain"]),
        "ligand_resseq": str(best["resseq"]),
        "ligand_icode": str(best["icode"]),
        "ligand_atom_count": int(best["atom_count"]),
        "ligand_heavy_atom_count": int(best["heavy_atom_count"]),
    }


def ligand_mask(atoms: pd.DataFrame, ligand: dict[str, object]) -> pd.Series:
    """Return mask for selected ligand atoms."""
    return (
        (atoms["record"] == "HETATM")
        & (atoms["resname"] == ligand["ligand_id"])
        & (atoms["chain"] == ligand["ligand_chain"])
        & (atoms["resseq"] == ligand["ligand_resseq"])
        & (atoms["icode"] == ligand["ligand_icode"])
    )


def binding_site_contacts(atoms: pd.DataFrame, ligand: dict[str, object], cutoff: float = 4.0) -> pd.DataFrame:
    """Calculate residue contacts within a cutoff of ligand atoms."""
    lig_atoms = atoms[ligand_mask(atoms, ligand)].reset_index(drop=True)
    protein = atoms[atoms["record"] == "ATOM"].reset_index(drop=True)
    if lig_atoms.empty or protein.empty:
        return pd.DataFrame()

    lig_xyz = lig_atoms[["x", "y", "z"]].to_numpy(float)
    prot_xyz = protein[["x", "y", "z"]].to_numpy(float)
    distances = np.sqrt(((prot_xyz[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(axis=2))
    close_pairs = np.argwhere(distances <= cutoff)
    residue_data: dict[tuple[str, str, str, str], dict[str, object]] = {}

    for protein_idx, ligand_idx in close_pairs:
        prot = protein.iloc[int(protein_idx)]
        lig = lig_atoms.iloc[int(ligand_idx)]
        key = (str(prot["chain"]), str(prot["resname"]), str(prot["resseq"]), str(prot["icode"]))
        entry = residue_data.setdefault(
            key,
            {
                "chain": key[0],
                "resname": key[1],
                "resseq": key[2],
                "icode": key[3],
                "min_distance_angstrom": float("inf"),
                "contact_atom_count": 0,
                "hydrophobic_contact": False,
                "hydrogen_bond_candidate": False,
                "aromatic_contact": False,
                "charged_contact": False,
            },
        )
        distance = float(distances[int(protein_idx), int(ligand_idx)])
        entry["min_distance_angstrom"] = min(float(entry["min_distance_angstrom"]), distance)
        entry["contact_atom_count"] = int(entry["contact_atom_count"]) + 1
        entry["hydrophobic_contact"] = bool(entry["hydrophobic_contact"] or prot["resname"] in HYDROPHOBIC_RESIDUES)
        entry["aromatic_contact"] = bool(entry["aromatic_contact"] or (prot["resname"] in AROMATIC_RESIDUES and distance <= 5.0))
        entry["charged_contact"] = bool(entry["charged_contact"] or prot["resname"] in CHARGED_RESIDUES)
        donor_acceptor_elements = {"N", "O", "S"}
        entry["hydrogen_bond_candidate"] = bool(
            entry["hydrogen_bond_candidate"]
            or (
                distance <= 3.5
                and str(prot["element"]).upper() in donor_acceptor_elements
                and str(lig["element"]).upper() in donor_acceptor_elements
            )
        )

    table = pd.DataFrame(residue_data.values()).sort_values(["min_distance_angstrom", "chain", "resseq"])
    return table.reset_index(drop=True)


def write_extracted_files(pdb_id: str, lines: list[str], atoms: pd.DataFrame, ligand: dict[str, object]) -> dict[str, str]:
    """Write protein and selected ligand PDB files."""
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    protein_path = PREP_DIR / f"{pdb_id}_protein.pdb"
    ligand_path = PREP_DIR / f"{pdb_id}_{ligand['ligand_id']}_ligand.pdb"

    protein_lines = [line for line in lines if line.startswith("ATOM")]
    protein_path.write_text("\n".join(protein_lines + ["END", ""]), encoding="utf-8")

    selected = atoms[ligand_mask(atoms, ligand)]
    ligand_lines = selected["line"].tolist()
    ligand_path.write_text("\n".join(ligand_lines + ["END", ""]), encoding="utf-8")

    return {
        "protein_pdb": str(protein_path.relative_to(PROJECT_ROOT)),
        "ligand_pdb": str(ligand_path.relative_to(PROJECT_ROOT)),
    }


def main() -> None:
    """Fetch structures, select ligands, extract files, and save metadata."""
    STRUCTURE_DIR.mkdir(parents=True, exist_ok=True)
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    binding_rows: list[pd.DataFrame] = []

    for candidate in PDB_CANDIDATES:
        pdb_id = candidate["pdb_id"]
        pdb_path = STRUCTURE_DIR / f"{pdb_id}.pdb"
        fetch_status = fetch_pdb(pdb_id, pdb_path)
        row: dict[str, object] = {
            "pdb_id": pdb_id,
            "description": candidate["description"],
            "fetch_status": fetch_status,
            "local_file": str(pdb_path.relative_to(PROJECT_ROOT)) if pdb_path.exists() else "",
            "resolution_angstrom": None,
            "chain_ids": "",
            "ligand_id": "",
            "ligand_atom_count": 0,
            "binding_site_residue_count": 0,
            "protein_pdb": "",
            "ligand_pdb": "",
            "parse_status": "not_parsed",
        }
        if not pdb_path.exists():
            rows.append(row)
            continue

        lines, atoms = parse_atoms(pdb_path)
        if atoms.empty:
            row["parse_status"] = "no_atoms_parsed"
            rows.append(row)
            continue

        row["resolution_angstrom"] = parse_resolution(lines)
        row["chain_ids"] = ",".join(sorted(str(item) for item in atoms.loc[atoms["record"] == "ATOM", "chain"].dropna().unique()))
        ligand = choose_ligand(atoms)
        if ligand is None:
            row["parse_status"] = "no_small_molecule_ligand_found"
            rows.append(row)
            continue

        contacts = binding_site_contacts(atoms, ligand)
        extracted = write_extracted_files(pdb_id, lines, atoms, ligand)
        row.update(ligand)
        row.update(extracted)
        row["binding_site_residue_count"] = int(len(contacts))
        row["parse_status"] = "parsed_with_ligand"

        if not contacts.empty:
            contacts.insert(0, "pdb_id", pdb_id)
            contacts.insert(1, "ligand_id", ligand["ligand_id"])
            binding_rows.append(contacts)
        rows.append(row)

    metadata = pd.DataFrame(rows)
    STRUCTURE_TABLE.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(STRUCTURE_TABLE, index=False)
    binding_site = pd.concat(binding_rows, ignore_index=True) if binding_rows else pd.DataFrame()
    binding_site.to_csv(BINDING_SITE_TABLE, index=False)

    usable = metadata[metadata["parse_status"] == "parsed_with_ligand"]
    status_counts = Counter(metadata["fetch_status"])
    parse_counts = Counter(metadata["parse_status"])
    metrics = {
        "candidate_structure_count": int(len(metadata)),
        "available_structure_count": int((metadata["local_file"].fillna("") != "").sum()),
        "parsed_cocrystal_count": int(len(usable)),
        "pdb_ids_used": usable["pdb_id"].tolist(),
        "ligand_ids_used": usable["ligand_id"].tolist(),
        "fetch_status_counts": dict(status_counts),
        "parse_status_counts": dict(parse_counts),
        "structure_table": str(STRUCTURE_TABLE.relative_to(PROJECT_ROOT)),
        "binding_site_table": str(BINDING_SITE_TABLE.relative_to(PROJECT_ROOT)),
    }
    save_json(METRICS_DIR / "egfr_structure_retrieval_metrics.json", metrics)

    lines_report = [
        "# EGFR Structure Retrieval and Extraction",
        "",
        f"- Candidate structures: {metrics['candidate_structure_count']}",
        f"- Parsed co-crystal structures with ligand: {metrics['parsed_cocrystal_count']}",
        f"- PDB IDs used: {', '.join(metrics['pdb_ids_used']) if metrics['pdb_ids_used'] else 'none'}",
        f"- Binding-site table: `{metrics['binding_site_table']}`",
        "",
        "No raw coordinate blocks are printed in this report.",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_structure_retrieval_report.md", "\n".join(lines_report))

    print(f"Candidate structures: {metrics['candidate_structure_count']}")
    print(f"Parsed co-crystals: {metrics['parsed_cocrystal_count']}")
    print(f"PDB IDs used: {','.join(metrics['pdb_ids_used']) if metrics['pdb_ids_used'] else 'none'}")


if __name__ == "__main__":
    main()
