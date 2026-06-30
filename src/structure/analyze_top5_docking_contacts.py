"""Analyze contacts for top-5 docked EGFR molecules against the 5UG9/8AM site."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, save_json  # noqa: E402


SELECTION_PATH = REPORTS_DIR / "egfr_top5_structure_selection.csv"
DOCKING_STATUS_PATH = REPORTS_DIR / "egfr_top5_docking_status.csv"
CASE_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_redocking_case.json"
REFERENCE_CONTACTS_PATH = PROJECT_ROOT / "data" / "processed" / "egfr_binding_site_residues.csv"
SANITY_TABLE_PATH = REPORTS_DIR / "egfr_top5_structure_sanity_table.csv"
METRICS_PATH = METRICS_DIR / "egfr_top5_structure_sanity_metrics.json"

HYDROPHOBIC_RESIDUES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
AROMATIC_RESIDUES = {"PHE", "TRP", "TYR", "HIS"}
CHARGED_RESIDUES = {"ASP", "GLU", "LYS", "ARG", "HIS"}
HBOND_RESIDUES = {"SER", "THR", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS", "CYS", "MET", "TRP"}


def read_json(path: Path) -> dict:
    """Read JSON with an empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def parse_pdbqt_coords(path: Path, first_pose_only: bool = True) -> np.ndarray:
    """Parse heavy-atom coordinates from a PDBQT file."""
    coords: list[list[float]] = []
    saw_model = False
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            if saw_model and first_pose_only:
                break
            saw_model = True
            continue
        if first_pose_only and saw_model and line.startswith("ENDMDL"):
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


def residue_id(chain: str, resname: str, resseq: str, icode: str = "") -> str:
    """Return a compact residue identifier."""
    suffix = str(icode).strip()
    return f"{str(chain).strip() or '_'}:{str(resname).strip()}{str(resseq).strip()}{suffix}"


def parse_protein_atoms(path: Path) -> pd.DataFrame:
    """Parse protein atom coordinates from PDB."""
    rows: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        chain = line[21:22].strip()
        resname = line[17:20].strip()
        resseq = line[22:26].strip()
        icode = line[26:27].strip()
        rows.append(
            {
                "chain": chain,
                "resname": resname,
                "resseq": resseq,
                "icode": icode,
                "residue_id": residue_id(chain, resname, resseq, icode),
                "element": infer_element(line),
                "x": x,
                "y": y,
                "z": z,
            }
        )
    return pd.DataFrame(rows)


def reference_contact_set() -> set[str]:
    """Load reference 5UG9/8AM binding-site contact residues."""
    if not REFERENCE_CONTACTS_PATH.exists():
        return set()
    contacts = pd.read_csv(REFERENCE_CONTACTS_PATH)
    contacts = contacts[(contacts["pdb_id"].astype(str) == "5UG9")]
    return {
        residue_id(row["chain"], row["resname"], row["resseq"], row.get("icode", ""))
        for _, row in contacts.iterrows()
    }


def contact_residues(protein: pd.DataFrame, ligand_coords: np.ndarray, cutoff: float = 4.0) -> pd.DataFrame:
    """Return residue contacts within cutoff of ligand atoms."""
    if protein.empty or ligand_coords.size == 0:
        return pd.DataFrame()
    protein_coords = protein[["x", "y", "z"]].to_numpy(dtype=float)
    distances = np.linalg.norm(protein_coords[:, None, :] - ligand_coords[None, :, :], axis=2)
    min_dist = distances.min(axis=1)
    contact_atoms = protein.loc[min_dist <= cutoff].copy()
    if contact_atoms.empty:
        return pd.DataFrame()
    contact_atoms["min_distance_angstrom"] = min_dist[min_dist <= cutoff]
    grouped = (
        contact_atoms.groupby(["residue_id", "chain", "resname", "resseq", "icode"], dropna=False)
        .agg(
            min_distance_angstrom=("min_distance_angstrom", "min"),
            contact_atom_count=("element", "size"),
        )
        .reset_index()
    )
    grouped["hydrophobic_contact_candidate"] = grouped["resname"].isin(HYDROPHOBIC_RESIDUES)
    grouped["hydrogen_bond_candidate"] = grouped["resname"].isin(HBOND_RESIDUES)
    grouped["aromatic_contact_candidate"] = grouped["resname"].isin(AROMATIC_RESIDUES)
    grouped["charged_contact_candidate"] = grouped["resname"].isin(CHARGED_RESIDUES)
    return grouped


def plausibility_label(
    docking_status: str,
    contact_count: int,
    shared_fraction: float | None,
    centroid_distance: float | None,
    medchem_alert_flag: bool,
) -> tuple[str, str, str]:
    """Assign pose plausibility and structure sanity labels."""
    if "completed" not in str(docking_status):
        return "not_evaluable", "structure_sanity_fail", "Docking or preparation failed."
    if centroid_distance is None or contact_count == 0:
        return "not_evaluable", "structure_sanity_fail", "No evaluable binding-site contacts."
    inside_site = centroid_distance <= 8.0
    has_overlap = shared_fraction is not None and shared_fraction > 0.20
    reasonable_contacts = contact_count >= 5
    if inside_site and has_overlap and reasonable_contacts and not medchem_alert_flag:
        return "inside_reference_site_with_shared_contacts", "structure_sanity_pass", "Docking pose localizes to the reference site and shares 8AM contacts."
    if inside_site and reasonable_contacts:
        return "inside_reference_site_borderline_overlap", "structure_sanity_warning", "Docking pose is in the reference site but contact overlap is limited or risk flags exist."
    return "weak_or_far_pose", "structure_sanity_warning", "Docking pose is evaluable but has weak pocket/contact support."


def main() -> None:
    """Analyze top-5 docking contacts and build structure sanity table."""
    if not SELECTION_PATH.exists() or not DOCKING_STATUS_PATH.exists():
        raise FileNotFoundError("Selection and docking status outputs are required before contact analysis.")
    case = read_json(CASE_PATH)
    protein_path = PROJECT_ROOT / case.get("protein_pdb", "data/structure_prepared/5UG9_protein.pdb")
    reference_ligand = PROJECT_ROOT / case.get("ligand_pdbqt", "data/structure_prepared/5UG9_8AM_ligand.pdbqt")
    protein = parse_protein_atoms(protein_path)
    ref_coords = parse_pdbqt_coords(reference_ligand, first_pose_only=False)
    ref_centroid = ref_coords.mean(axis=0) if ref_coords.size else None
    ref_contacts = reference_contact_set()

    selection = pd.read_csv(SELECTION_PATH)
    docking = pd.read_csv(DOCKING_STATUS_PATH)
    merged = selection.merge(docking, on=["molecule_id", "rank_before_docking", "scaffold_id"], how="left", validate="one_to_one")

    rows: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        pose_path = PROJECT_ROOT / str(row.get("pose_path", ""))
        docking_status = str(row.get("docking_status", "missing"))
        score = pd.to_numeric(row.get("vina_score_kcal_mol"), errors="coerce")
        ligand_coords = parse_pdbqt_coords(pose_path, first_pose_only=True) if pose_path.exists() else np.asarray([])
        centroid_distance = None
        contacts = pd.DataFrame()
        if ligand_coords.size and ref_centroid is not None:
            centroid_distance = float(np.linalg.norm(ligand_coords.mean(axis=0) - ref_centroid))
            contacts = contact_residues(protein, ligand_coords)
        contact_ids = set(contacts["residue_id"]) if not contacts.empty else set()
        shared = sorted(contact_ids & ref_contacts)
        shared_fraction = float(len(shared) / len(ref_contacts)) if ref_contacts else None
        medchem_alert = bool(row.get("medchem_alert_flag", False))
        pose_plausibility, sanity_label, note = plausibility_label(
            docking_status=docking_status,
            contact_count=len(contact_ids),
            shared_fraction=shared_fraction,
            centroid_distance=centroid_distance,
            medchem_alert_flag=medchem_alert,
        )
        inside_site = bool(centroid_distance is not None and centroid_distance <= 8.0 and len(contact_ids) >= 3)
        rows.append(
            {
                "molecule_id": row["molecule_id"],
                "rank_before_docking": int(row["rank_before_docking"]),
                "scaffold_id": row["scaffold_id"],
                "predicted_pIC50": row["predicted_pIC50"],
                "conformal_interval_width": row["conformal_interval_width"],
                "applicability_domain_bin": row["applicability_domain_bin"],
                "medchem_alert_flag": medchem_alert,
                "docking_status": docking_status,
                "vina_score_kcal_mol": None if pd.isna(score) else float(score),
                "contact_residue_count": int(len(contact_ids)),
                "contact_residue_list": ";".join(sorted(contact_ids)),
                "shared_contact_count_with_8AM": int(len(shared)),
                "shared_contact_fraction_with_8AM": shared_fraction,
                "distance_to_8AM_centroid": centroid_distance,
                "inside_reference_binding_site": inside_site,
                "pose_plausibility": pose_plausibility,
                "structure_sanity_label": sanity_label,
                "final_triage_note": note,
                "hydrophobic_contact_candidate_count": int(contacts["hydrophobic_contact_candidate"].sum()) if not contacts.empty else 0,
                "hydrogen_bond_candidate_count": int(contacts["hydrogen_bond_candidate"].sum()) if not contacts.empty else 0,
                "aromatic_contact_candidate_count": int(contacts["aromatic_contact_candidate"].sum()) if not contacts.empty else 0,
                "charged_contact_candidate_count": int(contacts["charged_contact_candidate"].sum()) if not contacts.empty else 0,
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(SANITY_TABLE_PATH, index=False)
    scores = pd.to_numeric(table["vina_score_kcal_mol"], errors="coerce").dropna()
    shared_fractions = pd.to_numeric(table["shared_contact_fraction_with_8AM"], errors="coerce").dropna()
    label_counts = table["structure_sanity_label"].value_counts().to_dict()
    metrics = {
        "top5_structure_sanity_status": "completed" if len(table) == 5 else "degraded_missing_selected_rows",
        "pdb_id": "5UG9",
        "reference_ligand_id": "8AM",
        "selected_molecule_count": int(len(table)),
        "successful_docking_count": int(table["docking_status"].astype(str).str.contains("completed").sum()),
        "best_docking_score_kcal_mol": float(scores.min()) if not scores.empty else None,
        "worst_docking_score_kcal_mol": float(scores.max()) if not scores.empty else None,
        "mean_shared_contact_fraction_with_8AM": float(shared_fractions.mean()) if not shared_fractions.empty else None,
        "structure_sanity_label_counts": {str(key): int(value) for key, value in label_counts.items()},
        "reference_contact_count_8AM": int(len(ref_contacts)),
        "sanity_table": str(SANITY_TABLE_PATH.relative_to(PROJECT_ROOT)),
        "contact_typing_policy": "heuristic residue-category counts from protein residues within 4 angstrom of docked ligand atoms.",
    }
    save_json(METRICS_PATH, metrics)

    print(f"Sanity table: {SANITY_TABLE_PATH}")
    print(f"Selected molecule rows: {len(table)}")
    print(f"Successful docking count: {metrics['successful_docking_count']}")


if __name__ == "__main__":
    main()
