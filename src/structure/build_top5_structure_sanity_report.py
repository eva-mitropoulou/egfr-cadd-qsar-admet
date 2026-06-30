"""Build report/figures and patch project docs for top-5 structure sanity docking."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import FIGURES_DIR, METRICS_DIR, REPORTS_DIR, markdown_table, save_figure, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


SANITY_TABLE_PATH = REPORTS_DIR / "egfr_top5_structure_sanity_table.csv"
SELECTION_PATH = REPORTS_DIR / "egfr_top5_structure_selection.csv"
REPORT_PATH = REPORTS_DIR / "egfr_top5_structure_sanity_report.md"
METRICS_PATH = METRICS_DIR / "egfr_top5_structure_sanity_metrics.json"
DOCKING_METRICS_PATH = METRICS_DIR / "egfr_top5_docking_metrics.json"
REDOCKING_METRICS_PATH = METRICS_DIR / "egfr_redocking_metrics.json"
FIGURE_PATH = FIGURES_DIR / "egfr_top5_docking_scores_and_contact_overlap.png"
OVERLAY_SCRIPT_PATH = REPORTS_DIR / "structure_visualization" / "top5_pose_overlay.pml"
OVERLAY_INSTRUCTIONS_PATH = FIGURES_DIR / "egfr_top5_pose_overlay_instructions.md"
FINAL_REPORT_PATH = REPORTS_DIR / "final_egfr_cadd_qsar_report.md"
CV_BULLETS_PATH = REPORTS_DIR / "final_egfr_cv_bullets.md"
PROJECT_CARD_PATH = PROJECT_ROOT / "portfolio_assets" / "egfr_project_card.md"
README_PATH = PROJECT_ROOT / "README.md"

REQUIRED_LIMITATION = (
    "Docking of top-ranked molecules was used as a structure-aware sanity check, "
    "not as proof of binding affinity, therapeutic efficacy, or prospective discovery."
)
DOC_CLAIM = (
    "Docked the top 5 clean, diverse, high-ranked existing EGFR molecules into the validated 5UG9 binding site as a structure-aware sanity check, "
    "reporting Vina scores, pocket localization, and shared contact residues with the 8AM reference ligand."
)
CV_ADDITION = (
    "Extended the EGFR triage workflow with a structure-aware sanity check by docking the top 5 clean, scaffold-diverse, high-ranked existing molecules "
    "into the validated 5UG9 binding site and comparing Vina scores, pocket localization, and shared contacts against the 8AM reference ligand."
)


def read_json(path: Path) -> dict:
    """Read JSON with an empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def write_overlay_fallback_script(table: pd.DataFrame) -> str:
    """Write a PyMOL helper script for manual pose review."""
    OVERLAY_SCRIPT_PATH.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "reinitialize",
        f"load {(PROJECT_ROOT / 'data/structure_prepared/5UG9_protein.pdb').resolve()}, receptor_5UG9",
        f"load {(PROJECT_ROOT / 'data/structure_prepared/5UG9_8AM_ligand.pdb').resolve()}, ligand_8AM_reference",
        "show cartoon, receptor_5UG9",
        "show sticks, ligand_8AM_reference",
        "color cyan, ligand_8AM_reference",
    ]
    for molecule_id in table["molecule_id"]:
        pose = PROJECT_ROOT / "data" / "structure_prepared" / "top5_docked" / f"{molecule_id}_vina_out.pdbqt"
        if pose.exists():
            object_name = f"pose_{molecule_id}".replace("-", "_")
            lines.extend([f"load {pose.resolve()}, {object_name}", f"show sticks, {object_name}"])
    lines.extend(["zoom ligand_8AM_reference", "set ray_opaque_background, off"])
    OVERLAY_SCRIPT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_text(
        OVERLAY_INSTRUCTIONS_PATH,
        "\n".join(
            [
                "# EGFR Top-5 Pose Overlay Instructions",
                "",
                "A static overlay grid was not generated automatically in this run.",
                f"Use `{OVERLAY_SCRIPT_PATH.relative_to(PROJECT_ROOT)}` in PyMOL to inspect top-5 docked poses against the 5UG9/8AM reference ligand.",
                "",
            ]
        ),
    )
    return str(OVERLAY_SCRIPT_PATH.relative_to(PROJECT_ROOT))


def save_primary_figure(table: pd.DataFrame) -> None:
    """Save score/contact-overlap comparison figure."""
    display = table.copy()
    display["molecule_label"] = display["molecule_id"].astype(str)
    scores = pd.to_numeric(display["vina_score_kcal_mol"], errors="coerce")
    overlap = pd.to_numeric(display["shared_contact_fraction_with_8AM"], errors="coerce")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].bar(display["molecule_label"], scores, color="#4C78A8")
    axes[0].axhline(-9.471, color="black", linestyle="--", linewidth=1, label="8AM redocking score")
    axes[0].set_ylabel("Vina score (kcal/mol)")
    axes[0].set_title("Top-5 Docking Scores")
    axes[0].tick_params(axis="x", rotation=30)
    axes[0].legend(fontsize=8)

    axes[1].bar(display["molecule_label"], overlap, color="#59A14F")
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Shared contact fraction with 8AM")
    axes[1].set_title("5UG9 Contact Overlap")
    axes[1].tick_params(axis="x", rotation=30)
    fig.tight_layout()
    save_figure(FIGURE_PATH)


def replace_or_append_section(text: str, title: str, section: str, before_title: str | None = None) -> str:
    """Replace an existing section or insert before a known section."""
    if title in text:
        start = text.index(title)
        next_start = text.find("\n## ", start + len(title))
        if next_start == -1:
            return text[:start].rstrip() + "\n\n" + section.rstrip() + "\n"
        return text[:start].rstrip() + "\n\n" + section.rstrip() + "\n" + text[next_start:]
    if before_title and before_title in text:
        return text.replace(before_title, section.rstrip() + "\n\n" + before_title)
    return text.rstrip() + "\n\n" + section.rstrip() + "\n"


def patch_final_docs(metrics: dict) -> list[str]:
    """Patch final project report, CV bullets, project card, and README."""
    updated: list[str] = []
    pass_count = metrics.get("structure_sanity_label_counts", {}).get("structure_sanity_pass", 0)
    warning_count = metrics.get("structure_sanity_label_counts", {}).get("structure_sanity_warning", 0)
    fail_count = metrics.get("structure_sanity_label_counts", {}).get("structure_sanity_fail", 0)
    summary_section = "\n".join(
        [
            "## Top-5 Structure Sanity Docking",
            "",
            DOC_CLAIM,
            "",
            f"- Successful dockings: {metrics.get('successful_docking_count')}/{metrics.get('selected_molecule_count')}",
            f"- Best/worst Vina score among successful dockings: {metrics.get('best_docking_score_kcal_mol')} / {metrics.get('worst_docking_score_kcal_mol')} kcal/mol",
            f"- Mean shared contact fraction with 8AM: {metrics.get('mean_shared_contact_fraction_with_8AM')}",
            f"- Structure sanity labels: pass {pass_count}, warning {warning_count}, fail {fail_count}",
            "",
            REQUIRED_LIMITATION,
            "",
        ]
    )

    if FINAL_REPORT_PATH.exists():
        text = FINAL_REPORT_PATH.read_text(encoding="utf-8")
        text = replace_or_append_section(text, "## Top-5 Structure Sanity Docking", summary_section, "## Limitations")
        FINAL_REPORT_PATH.write_text(text, encoding="utf-8")
        updated.append(str(FINAL_REPORT_PATH.relative_to(PROJECT_ROOT)))

    if CV_BULLETS_PATH.exists():
        text = CV_BULLETS_PATH.read_text(encoding="utf-8")
        if CV_ADDITION not in text:
            text = text.rstrip() + f"\n- {CV_ADDITION}\n"
            CV_BULLETS_PATH.write_text(text, encoding="utf-8")
        updated.append(str(CV_BULLETS_PATH.relative_to(PROJECT_ROOT)))

    if PROJECT_CARD_PATH.exists():
        text = PROJECT_CARD_PATH.read_text(encoding="utf-8")
        if "Top-5 structure sanity check" not in text:
            text = text.replace(
                "## Positioning",
                f"- Top-5 structure sanity check: {metrics.get('successful_docking_count')}/5 molecules docked in the 5UG9 pocket; mean shared contact fraction with 8AM {metrics.get('mean_shared_contact_fraction_with_8AM')}\n\n## Positioning",
            )
            PROJECT_CARD_PATH.write_text(text, encoding="utf-8")
        updated.append(str(PROJECT_CARD_PATH.relative_to(PROJECT_ROOT)))

    if README_PATH.exists():
        text = README_PATH.read_text(encoding="utf-8")
        if "Top-5 structure sanity docking" not in text:
            text = text.replace(
                "| Structure module | Runs co-crystal contact analysis and one retrospective redocking pose-recovery audit. |",
                "| Structure module | Runs co-crystal contact analysis, one retrospective redocking pose-recovery audit, and a top-5 structure-aware sanity check. |",
            )
            text = text.replace(
                "- Retrospective Vina redocking pose-recovery audit on 5UG9 with ligand 8AM with a -9.471 kcal/mol score and 0.968 A\n  pose-recovery RMSD.",
                "- Retrospective Vina redocking pose-recovery audit on 5UG9 with ligand 8AM with a -9.471 kcal/mol score and 0.968 A\n  pose-recovery RMSD.\n- Top-5 structure sanity docking of clean, diverse, high-ranked existing molecules into the validated 5UG9 pocket.",
            )
            text = text.replace(
                "| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 A |",
                "| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 A |\n"
                f"| Top-5 structure sanity docking | {metrics.get('successful_docking_count')}/5 docked; mean shared contact fraction {metrics.get('mean_shared_contact_fraction_with_8AM')} |",
            )
            text = text.replace(
                "| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 angstrom |",
                "| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 angstrom |\n"
                f"| Top-5 structure sanity check | {metrics.get('successful_docking_count')}/5 docked; scores {metrics.get('best_docking_score_kcal_mol')} to {metrics.get('worst_docking_score_kcal_mol')} kcal/mol |",
            )
            text = text.replace(
                "The redocking result is a retrospective pose-recovery check on a known co-crystal case. It is useful structure-based context, but it is not a prospective docking campaign.",
                "The redocking result is a retrospective pose-recovery check on a known co-crystal case. The top-5 docking step uses that validated 5UG9 pocket to sanity-check existing ranked molecules by Vina score, pocket localization, and shared contacts with the 8AM reference ligand.",
            )
            text = text.replace(
                "ADMET-style triage uses simple drug-likeness and model-risk proxy rules.",
                "ADMET-style triage uses simple drug-likeness and model-risk proxy rules.\n- Top-5 docking is a structure-aware sanity check over existing ranked molecules, not binding confirmation.",
            )
            text = text.replace(
                "- `reports/egfr_redocking_audit_report.md`",
                "- `reports/egfr_redocking_audit_report.md`\n- `reports/egfr_top5_structure_sanity_report.md`\n- `reports/egfr_top5_structure_sanity_table.csv`",
            )
            README_PATH.write_text(text, encoding="utf-8")
        updated.append(str(README_PATH.relative_to(PROJECT_ROOT)))
    return updated


def main() -> None:
    """Build top-5 structure sanity report and patch project docs."""
    if not SANITY_TABLE_PATH.exists():
        raise FileNotFoundError(f"Missing top-5 sanity table: {SANITY_TABLE_PATH}")
    table = pd.read_csv(SANITY_TABLE_PATH)
    metrics = read_json(METRICS_PATH)
    docking_metrics = read_json(DOCKING_METRICS_PATH)
    redocking = read_json(REDOCKING_METRICS_PATH)

    save_primary_figure(table)
    overlay_script = write_overlay_fallback_script(table)
    patched_files = patch_final_docs(metrics)

    display_columns = [
        "molecule_id",
        "rank_before_docking",
        "predicted_pIC50",
        "docking_status",
        "vina_score_kcal_mol",
        "shared_contact_fraction_with_8AM",
        "distance_to_8AM_centroid",
        "structure_sanity_label",
    ]
    report_table = table[display_columns].copy()
    lines = [
        "# EGFR Top-5 Structure Sanity Report",
        "",
        "## Purpose",
        "",
        "This module links the final ranked EGFR molecule table to the validated 5UG9 structure workflow.",
        REQUIRED_LIMITATION,
        "",
        "## Selection Criteria",
        "",
        "The top-5 molecules were selected from existing ranked EGFR records using high triage score, applicability-domain support, low/acceptable uncertainty, no PAINS/Brenk/unwanted-substructure alert, acceptable drug-likeness, and scaffold diversity where possible.",
        f"Selection table: `{SELECTION_PATH.relative_to(PROJECT_ROOT)}`",
        "",
        "## Reference Redocking Context",
        "",
        "- PDB ID: 5UG9",
        "- Reference ligand: 8AM",
        f"- Reference redocking score: {redocking.get('docking_score_kcal_mol', -9.471)} kcal/mol",
        f"- Reference redocking RMSD: {redocking.get('pose_recovery_rmsd_angstrom', 0.968)} A",
        "",
        "## Docking Setup",
        "",
        f"- Receptor PDBQT: `{docking_metrics.get('receptor_pdbqt')}`",
        f"- Docking box center: {docking_metrics.get('docking_box_center')}",
        f"- Docking box size: {docking_metrics.get('docking_box_size')}",
        "- Backend preference: Python Vina API, then Vina CLI fallback",
        "",
        "## Top-5 Results",
        "",
        markdown_table(report_table),
        "",
        "## Contact Overlap With 8AM",
        "",
        f"- Mean shared contact fraction with 8AM: {metrics.get('mean_shared_contact_fraction_with_8AM')}",
        f"- Reference 8AM contact count: {metrics.get('reference_contact_count_8AM')}",
        f"- Primary figure: `{FIGURE_PATH.relative_to(PROJECT_ROOT)}`",
        f"- Overlay helper script: `{overlay_script}`",
        "",
        "## Structure-Aware Sanity Labels",
        "",
    ]
    for label, count in metrics.get("structure_sanity_label_counts", {}).items():
        lines.append(f"- {label}: {count}")
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- Docking scores are Vina scoring-function outputs and should not be interpreted as physical binding energies.",
            "- Contact classes are heuristic residue-contact annotations.",
            "- This is a retrospective sanity check over existing ranked molecules.",
            "- Molecules that pass this check remain computationally prioritized existing records, not experimentally validated binding or inhibition claims.",
            "",
            "## Reproducibility",
            "",
            "```bash",
            "python scripts/agentic_top5_structure_sanity.py",
            "```",
            "",
        ]
    )
    write_text(REPORT_PATH, "\n".join(lines))

    metrics.update(
        {
            "top5_structure_sanity_report": str(REPORT_PATH.relative_to(PROJECT_ROOT)),
            "top5_structure_sanity_figure": str(FIGURE_PATH.relative_to(PROJECT_ROOT)),
            "top5_pose_overlay_script": overlay_script,
            "patched_files": patched_files,
        }
    )
    METRICS_PATH.write_text(json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8")

    print(f"Top-5 structure report: {REPORT_PATH}")
    print(f"Primary figure: {FIGURE_PATH}")
    print(f"Patched files: {len(patched_files)}")


if __name__ == "__main__":
    main()
