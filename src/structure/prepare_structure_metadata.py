"""Prepare EGFR structure metadata report from parsed co-crystal structures."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import DATA_DIR, FIGURES_DIR, METRICS_DIR, REPORTS_DIR, markdown_table, read_json, save_figure, save_json, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


STRUCTURE_TABLE = DATA_DIR / "processed" / "egfr_structure_candidates.csv"
BINDING_SITE_TABLE = DATA_DIR / "processed" / "egfr_binding_site_residues.csv"


def main() -> None:
    """Create structure metadata report and workflow figure."""
    if not STRUCTURE_TABLE.exists():
        raise FileNotFoundError(f"Missing structure table: {STRUCTURE_TABLE}")

    structures = pd.read_csv(STRUCTURE_TABLE)
    binding = pd.read_csv(BINDING_SITE_TABLE) if BINDING_SITE_TABLE.exists() else pd.DataFrame()
    usable = structures[structures["parse_status"] == "parsed_with_ligand"].copy()
    retrieval_metrics = read_json(METRICS_DIR / "egfr_structure_retrieval_metrics.json")

    existing_metrics = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
    if existing_metrics.get("redocking_status") == "completed_redocking":
        status = "completed_redocking"
    elif existing_metrics.get("redocking_status") and len(usable) > 0 and not binding.empty:
        status = "structure_analysis_completed_redocking_failed"
    elif len(usable) > 0 and not binding.empty:
        status = "structure_analysis_completed_redocking_pending"
    else:
        status = "metadata_only_degraded"

    metrics = {
        "structure_candidates": int(len(structures)),
        "available_structures": int((structures["local_file"].fillna("") != "").sum()),
        "parsed_cocrystal_count": int(len(usable)),
        "pdb_ids_used": usable["pdb_id"].tolist(),
        "ligand_ids_used": usable["ligand_id"].tolist(),
        "fetch_status_counts": retrieval_metrics.get("fetch_status_counts", {}),
        "parse_status_counts": retrieval_metrics.get("parse_status_counts", {}),
        "structure_module_status": status,
        "binding_site_residue_rows": int(len(binding)),
        "docking_ready": False,
    }
    for key in [
        "redocking_status",
        "redocking_reason",
        "vina_available",
        "redocking_case_status",
        "redocking_case_pdb_id",
        "redocking_case_ligand_id",
        "protein_conversion_status",
        "ligand_conversion_status",
        "interaction_analysis_status",
        "interaction_fingerprint_status",
        "interaction_residue_count",
        "interaction_category_counts",
    ]:
        if key in existing_metrics:
            metrics[key] = existing_metrics[key]
    save_json(METRICS_DIR / "egfr_structure_module_metrics.json", metrics)

    steps = ["fetched", "ligand parsed", "contacts", "redocking"]
    values = [
        int(metrics["available_structures"] > 0),
        int(metrics["parsed_cocrystal_count"] > 0),
        int(metrics["binding_site_residue_rows"] > 0),
        0,
    ]
    plt.figure(figsize=(7, 4.5))
    plt.bar(steps, values, color=["#4C78A8" if value else "#BAB0AC" for value in values])
    plt.ylim(0, 1.2)
    plt.ylabel("Completed")
    plt.title("EGFR Structure Workflow Status")
    save_figure(FIGURES_DIR / "egfr_structure_workflow.png")

    display = usable[
        [
            "pdb_id",
            "ligand_id",
            "chain_ids",
            "resolution_angstrom",
            "ligand_atom_count",
            "binding_site_residue_count",
            "parse_status",
        ]
    ].copy()
    lines = [
        "# EGFR Structure Metadata Report",
        "",
        "This module uses real EGFR co-crystal structures with selected bound small-molecule ligands when available.",
        "",
        f"- Structure candidates considered: {metrics['structure_candidates']}",
        f"- Available structures: {metrics['available_structures']}",
        f"- Parsed co-crystals with ligand: {metrics['parsed_cocrystal_count']}",
        f"- Binding-site residue rows: {metrics['binding_site_residue_rows']}",
        f"- Structure module status: {metrics['structure_module_status']}",
        "",
        "## Parsed Structure Metadata",
        "",
        markdown_table(display) if not display.empty else "No parsed co-crystal structure was available.",
        "",
        "No raw coordinate blocks are printed in this report.",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_structure_metadata_report.md", "\n".join(lines))

    print(f"Structure module status: {metrics['structure_module_status']}")
    print(f"Parsed co-crystals: {metrics['parsed_cocrystal_count']}")


if __name__ == "__main__":
    main()
