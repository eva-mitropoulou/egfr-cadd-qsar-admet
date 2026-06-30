"""Patch final public-facing EGFR reports with medchem-alert evidence."""

from __future__ import annotations

import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REPORTS_DIR = PROJECT_ROOT / "reports"
METRICS_DIR = REPORTS_DIR / "metrics"
DOCS_DIR = PROJECT_ROOT / "docs"


REQUIRED_WORDING = (
    "PAINS, Brenk, and external unwanted-substructure SMARTS alerts were used as "
    "medicinal-chemistry risk annotations and sensitivity-analysis filters, not "
    "automatic exclusions from the primary EGFR QSAR benchmark."
)


def read_json(path: Path) -> dict:
    """Read JSON with an empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def pct(value: object) -> str:
    """Format a fraction as a percentage."""
    try:
        return f"{float(value):.1%}"
    except Exception:
        return "unavailable"


def fmt_int(value: object) -> str:
    """Format integer-like values."""
    try:
        return f"{int(value):,}"
    except Exception:
        return "unavailable"


def scaffold_rmse_delta(sensitivity: dict) -> tuple[float | None, str]:
    """Return strict-clean scaffold RMSE delta versus full sensitivity baseline."""
    rows = sensitivity.get("performance_rows", [])
    full = next(
        (
            row
            for row in rows
            if row.get("subset") == "main_all_model_ready" and row.get("split") == "scaffold_split"
        ),
        None,
    )
    strict = next(
        (
            row
            for row in rows
            if row.get("subset") == "strict_medchem_clean" and row.get("split") == "scaffold_split"
        ),
        None,
    )
    if not full or not strict:
        return None, "Strict medchem-clean scaffold-split sensitivity was unavailable."
    delta = float(strict["RMSE"]) - float(full["RMSE"])
    if abs(delta) < 0.05:
        interpretation = "Strict medchem-clean sensitivity did not materially change the scaffold-split RMSE conclusion."
    elif delta > 0:
        interpretation = "Strict medchem-clean sensitivity produced a higher scaffold-split RMSE than the full sensitivity baseline."
    else:
        interpretation = "Strict medchem-clean sensitivity produced a lower scaffold-split RMSE than the full sensitivity baseline."
    return delta, interpretation


def replace_or_insert_section(text: str, title: str, section: str, before_title: str) -> str:
    """Replace an existing top-level section or insert it before another section."""
    if title in text:
        start = text.index(title)
        next_start = text.find("\n## ", start + len(title))
        if next_start == -1:
            return text[:start].rstrip() + "\n\n" + section.rstrip() + "\n"
        return text[:start].rstrip() + "\n\n" + section.rstrip() + "\n" + text[next_start:]
    if before_title in text:
        return text.replace(before_title, section.rstrip() + "\n\n" + before_title)
    return text.rstrip() + "\n\n" + section.rstrip() + "\n"


def patch_final_report(alerts: dict, sensitivity: dict, triage: dict) -> None:
    """Patch the final report with alert evidence."""
    path = REPORTS_DIR / "final_egfr_cadd_qsar_report.md"
    text = path.read_text(encoding="utf-8")
    delta, interpretation = scaffold_rmse_delta(sensitivity)
    composition = sensitivity.get("top20_candidate_composition", {})
    section = "\n".join(
        [
            "## Medicinal-Chemistry Alerts And Sensitivity",
            "",
            REQUIRED_WORDING,
            "",
            f"- PAINS-flagged molecules: {fmt_int(alerts.get('pains_flagged_count'))} ({pct(alerts.get('pains_flagged_fraction'))})",
            f"- Brenk-flagged molecules: {fmt_int(alerts.get('brenk_flagged_count'))} ({pct(alerts.get('brenk_flagged_fraction'))})",
            f"- External unwanted-substructure flagged molecules: {fmt_int(alerts.get('unwanted_substructure_flagged_count'))} ({pct(alerts.get('unwanted_substructure_flagged_fraction'))})",
            f"- Any medicinal-chemistry alert: {fmt_int(alerts.get('combined_medchem_alert_count'))} ({pct(alerts.get('combined_medchem_alert_fraction'))})",
            f"- External SMARTS catalog found: {alerts.get('external_csv_found')}",
            f"- Top-20 molecules without medchem alerts: {composition.get('top20_clean_count', triage.get('top20_clean_medchem_alert_count', 'unavailable'))}/20",
            f"- Diverse top-20 molecules without medchem alerts: {composition.get('diverse_top20_clean_count', triage.get('diverse_top20_clean_medchem_alert_count', 'unavailable'))}/20",
            (
                f"- Strict medchem-clean scaffold RMSE delta versus full sensitivity baseline: {delta:.3f}"
                if delta is not None
                else "- Strict medchem-clean scaffold RMSE delta: unavailable"
            ),
            f"- Sensitivity interpretation: {interpretation}",
            "",
            "Alert-containing molecules remain visible in triage outputs and are not treated as proven false positives or assay artifacts.",
            "",
        ]
    )
    text = replace_or_insert_section(text, "## Medicinal-Chemistry Alerts And Sensitivity", section, "## Structure-Based Module")
    path.write_text(text, encoding="utf-8")


def patch_cv_bullets() -> None:
    """Write conservative CV bullets with the new alert evidence."""
    path = REPORTS_DIR / "final_egfr_cv_bullets.md"
    bullet = (
        "Built a retrospective EGFR CADD/QSAR decision workflow from ChEMBL, "
        "curating 26,600 IC50 records into 10,593 model-ready molecules; "
        "benchmarked RDKit descriptor, Morgan fingerprint, and GPU PyTorch GNN "
        "models under random, scaffold, assay-aware, and document-aware validation, "
        "with Morgan RF achieving scaffold-split RMSE 0.871/R2 0.550; added "
        "applicability-domain analysis, split-conformal uncertainty, PAINS/Brenk/"
        "unwanted-substructure risk annotation and sensitivity analysis, "
        "ADMET-style/model-risk triage, SAR/error analysis, active-learning "
        "simulation, CLI prediction, ligand-contact analysis across four EGFR PDB "
        "structures, and 5UG9 redocking validation recovering the co-crystal ligand "
        "pose at 0.968 A RMSD."
    )
    path.write_text("# EGFR Project CV Bullets\n\n- " + bullet + "\n", encoding="utf-8")


def patch_project_card(alerts: dict, sensitivity: dict, triage: dict) -> None:
    """Patch the project card in a factual, organic tone."""
    composition = sensitivity.get("top20_candidate_composition", {})
    path = DOCS_DIR / "project_card.md"
    content = "\n".join(
        [
            "# EGFR CADD and QSAR Decision Workflow",
            "",
            "Retrospective EGFR inhibitor-like molecule prioritization using ChEMBL, RDKit, Morgan fingerprints, scaffold validation, uncertainty, applicability-domain analysis, medicinal-chemistry alerts, and ADMET-style/model-risk triage.",
            "",
            "## Snapshot",
            "",
            "- 26,600 raw EGFR IC50 activity rows",
            "- 10,593 model-ready molecules",
            "- Best scaffold-split QSAR model: Morgan Random Forest, R2 0.550",
            "- Applicability-domain MAE changed from 0.957 at low similarity to 0.513 at high similarity",
            f"- PAINS-flagged molecules: {fmt_int(alerts.get('pains_flagged_count'))} ({pct(alerts.get('pains_flagged_fraction'))})",
            f"- Brenk-flagged molecules: {fmt_int(alerts.get('brenk_flagged_count'))} ({pct(alerts.get('brenk_flagged_fraction'))})",
            f"- Top-20 without medchem alerts: {composition.get('top20_clean_count', triage.get('top20_clean_medchem_alert_count', 'unavailable'))}/20",
            "- Structure module: four EGFR co-crystals parsed; retrospective 5UG9/8AM Vina redocking pose-recovery RMSD 0.968 A",
            "- Exploratory custom PyTorch dense GCN benchmark did not beat the Morgan RF baseline",
            "",
            "## Positioning",
            "",
            "A complete, model-risk-aware CADD and QSAR workflow for existing public EGFR records. No molecule generation or efficacy claim.",
            "",
            REQUIRED_WORDING,
            "",
        ]
    )
    path.write_text(content, encoding="utf-8")


def patch_readme(alerts: dict, sensitivity: dict, triage: dict) -> None:
    """Patch README with the medchem-alert layer."""
    path = PROJECT_ROOT / "README.md"
    text = path.read_text(encoding="utf-8")
    if "PAINS, Brenk, and optional external unwanted-substructure SMARTS alert annotations." not in text:
        text = text.replace(
            "- ADMET-style and model-risk-aware ranking over existing molecules.",
            "- ADMET-style and model-risk-aware ranking over existing molecules.\n"
            "- PAINS, Brenk, and optional external unwanted-substructure SMARTS alert annotations.",
        )
    rows = [
        f"| PAINS-flagged molecules | {fmt_int(alerts.get('pains_flagged_count'))} ({pct(alerts.get('pains_flagged_fraction'))}) |",
        f"| Brenk-flagged molecules | {fmt_int(alerts.get('brenk_flagged_count'))} ({pct(alerts.get('brenk_flagged_fraction'))}) |",
        f"| External unwanted-substructure flagged molecules | {fmt_int(alerts.get('unwanted_substructure_flagged_count'))} ({pct(alerts.get('unwanted_substructure_flagged_fraction'))}) |",
        f"| Top-20 without medchem alerts | {sensitivity.get('top20_candidate_composition', {}).get('top20_clean_count', triage.get('top20_clean_medchem_alert_count', 'unavailable'))}/20 |",
    ]
    if "| PAINS-flagged molecules |" not in text:
        text = text.replace(
            "| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 A |",
            "| Redocking case | 5UG9 with ligand 8AM, RMSD 0.968 A |\n" + "\n".join(rows),
        )
    if REQUIRED_WORDING not in text:
        text = text.replace(
            "## Interpretation Context",
            "## Medicinal-Chemistry Alert Context\n\n" + REQUIRED_WORDING + "\n\n## Interpretation Context",
        )
    path.write_text(text, encoding="utf-8")


def main() -> None:
    """Patch all public-facing project summaries."""
    alerts = read_json(METRICS_DIR / "egfr_medchem_alerts_metrics.json")
    sensitivity = read_json(METRICS_DIR / "egfr_medchem_alert_sensitivity_metrics.json")
    triage = read_json(METRICS_DIR / "egfr_candidate_triage_metrics.json")
    if not alerts:
        raise FileNotFoundError("Missing egfr_medchem_alerts_metrics.json")

    patch_final_report(alerts, sensitivity, triage)
    patch_cv_bullets()
    patch_project_card(alerts, sensitivity, triage)
    patch_readme(alerts, sensitivity, triage)

    print(f"Updated: {REPORTS_DIR / 'final_egfr_cadd_qsar_report.md'}")
    print(f"Updated: {REPORTS_DIR / 'final_egfr_cv_bullets.md'}")
    print(f"Updated: {DOCS_DIR / 'project_card.md'}")
    print(f"Updated: {PROJECT_ROOT / 'README.md'}")


if __name__ == "__main__":
    main()
