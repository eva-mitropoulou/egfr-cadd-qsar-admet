"""Analyze EGFR ligand binding-site contacts from parsed co-crystal structures."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import DATA_DIR, FIGURES_DIR, METRICS_DIR, REPORTS_DIR, markdown_table, read_json, save_figure, save_json, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


BINDING_SITE_TABLE = DATA_DIR / "processed" / "egfr_binding_site_residues.csv"


def main() -> None:
    """Summarize binding-site residue contacts and heuristic interaction classes."""
    if not BINDING_SITE_TABLE.exists():
        raise FileNotFoundError(f"Missing binding-site table: {BINDING_SITE_TABLE}")
    contacts = pd.read_csv(BINDING_SITE_TABLE)

    if contacts.empty:
        metrics = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
        metrics.update(
            {
                "interaction_analysis_status": "failed_no_binding_site_contacts",
                "interaction_residue_count": 0,
                "interaction_fingerprint_status": "failed_no_contacts",
            }
        )
        save_json(METRICS_DIR / "egfr_structure_module_metrics.json", metrics)
        write_text(REPORTS_DIR / "egfr_interaction_fingerprint_report.md", "# EGFR Interaction Fingerprint Report\n\nNo binding-site contacts were available.\n")
        print("Interaction analysis status: failed_no_binding_site_contacts")
        return

    category_columns = [
        "hydrophobic_contact",
        "hydrogen_bond_candidate",
        "aromatic_contact",
        "charged_contact",
    ]
    category_counts = {column: int(contacts[column].sum()) for column in category_columns if column in contacts.columns}
    pdb_counts = contacts.groupby("pdb_id").size().reset_index(name="binding_site_residue_count")
    residue_counts = contacts["resname"].value_counts().reset_index()
    residue_counts.columns = ["resname", "count"]

    metrics = read_json(METRICS_DIR / "egfr_structure_module_metrics.json")
    metrics.update(
        {
            "interaction_analysis_status": "completed",
            "interaction_fingerprint_status": "completed_heuristic_contacts",
            "interaction_residue_count": int(len(contacts)),
            "interaction_category_counts": category_counts,
            "binding_site_table": str(BINDING_SITE_TABLE.relative_to(PROJECT_ROOT)),
        }
    )
    if metrics.get("redocking_status") == "completed_redocking":
        metrics["structure_module_status"] = "completed_redocking"
    elif metrics.get("parsed_cocrystal_count", 0) or len(contacts) > 0:
        metrics["structure_module_status"] = "structure_analysis_completed_redocking_failed"
    save_json(METRICS_DIR / "egfr_structure_module_metrics.json", metrics)

    plt.figure(figsize=(7, 4.5))
    plt.bar([key.replace("_", " ") for key in category_counts], list(category_counts.values()), color="#4C78A8")
    plt.ylabel("Residue-contact count")
    plt.title("EGFR Interaction Category Frequency")
    plt.xticks(rotation=20, ha="right")
    save_figure(FIGURES_DIR / "interaction_frequency.png")

    display_residues = residue_counts.head(12).copy()
    report = [
        "# EGFR Interaction Fingerprint Report",
        "",
        "Binding-site contacts were calculated from protein atoms within 4 angstrom of the selected co-crystal ligand.",
        "Contact classes are heuristic and should not be interpreted as rigorous interaction-energy decomposition.",
        "",
        "## Structure Contact Counts",
        "",
        markdown_table(pdb_counts),
        "",
        "## Heuristic Contact Categories",
        "",
    ]
    for key, value in category_counts.items():
        report.append(f"- {key}: {value}")
    report.extend(
        [
            "",
            "## Most Frequent Binding-Site Residue Names",
            "",
            markdown_table(display_residues),
            "",
            f"Detailed contact table: `{BINDING_SITE_TABLE.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    write_text(REPORTS_DIR / "egfr_interaction_fingerprint_report.md", "\n".join(report))

    print("Interaction analysis status: completed")
    print(f"Interaction residue count: {len(contacts)}")


if __name__ == "__main__":
    main()
