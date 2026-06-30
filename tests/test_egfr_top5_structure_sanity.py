from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCORES_TABLE = ROOT / "reports" / "egfr_top5_docking_scores.csv"
TOP5_REPORT = ROOT / "reports" / "egfr_top5_structure_sanity_report.md"
FINAL_REPORT = ROOT / "reports" / "final_egfr_cadd_qsar_report.md"
FIGURE = ROOT / "reports" / "figures" / "egfr_top5_vina_scores.png"

REQUIRED_COLUMNS = {
    "molecule_id",
    "rank_before_docking",
    "scaffold_id",
    "predicted_pIC50",
    "conformal_interval_width",
    "applicability_domain_bin",
    "medchem_alert_flag",
    "docking_status",
    "vina_score_kcal_mol",
    "docking_note",
}

FORBIDDEN_COLUMN_PATTERNS = [
    "shared_contact",
    "8am",
    "contact_fraction",
    "pose_plausibility",
    "structure_sanity_warning",
]

BANNED_REPORT_PHRASES = [
    "shared contact fraction",
    "contact overlap with 8am",
    "confirmed binding",
    "validated top candidates",
    "therapeutic efficacy",
    "prospective discovery",
]


def public_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace").lower().replace("-", " ")


def test_top5_docking_score_table_exists_and_has_five_rows():
    assert SCORES_TABLE.exists()
    table = pd.read_csv(SCORES_TABLE)
    assert len(table) == 5
    assert REQUIRED_COLUMNS.issubset(set(table.columns))


def test_score_table_has_no_raw_smiles_or_overlap_columns():
    table = pd.read_csv(SCORES_TABLE, nrows=0)
    lowered = [column.lower() for column in table.columns]
    assert "smiles" not in lowered
    assert "canonical_smiles" not in lowered
    assert "standardized_smiles" not in lowered
    offenders = [
        column
        for column in lowered
        if any(pattern in column for pattern in FORBIDDEN_COLUMN_PATTERNS)
    ]
    assert not offenders


def test_successful_dockings_have_numeric_vina_scores():
    table = pd.read_csv(SCORES_TABLE)
    successful = table[table["docking_status"].astype(str).str.contains("completed")]
    assert not successful.empty
    assert pd.to_numeric(successful["vina_score_kcal_mol"], errors="coerce").notna().all()


def test_reports_and_vina_score_figure_exist():
    assert TOP5_REPORT.exists()
    assert FINAL_REPORT.exists()
    assert FIGURE.exists()


def test_report_mentions_score_only_triage_language():
    text = public_text(TOP5_REPORT)
    assert "vina scores" in text
    assert "structure aware triage annotations" in text


def test_public_reports_do_not_contain_removed_overlap_or_overclaims():
    offenders = []
    for path in [TOP5_REPORT, FINAL_REPORT, ROOT / "README.md", ROOT / "portfolio_assets" / "egfr_project_card.md"]:
        if not path.exists():
            continue
        text = public_text(path)
        for phrase in BANNED_REPORT_PHRASES:
            if phrase in text:
                offenders.append(f"{path.relative_to(ROOT)}:{phrase}")
    assert not offenders
