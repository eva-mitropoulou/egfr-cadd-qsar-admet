from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "reports" / "egfr_top5_structure_selection.csv"
SANITY_TABLE = ROOT / "reports" / "egfr_top5_structure_sanity_table.csv"
TOP5_REPORT = ROOT / "reports" / "egfr_top5_structure_sanity_report.md"
FINAL_REPORT = ROOT / "reports" / "final_egfr_cadd_qsar_report.md"
FIGURE = ROOT / "reports" / "figures" / "egfr_top5_docking_scores_and_contact_overlap.png"

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
    "shared_contact_count_with_8AM",
    "shared_contact_fraction_with_8AM",
    "structure_sanity_label",
}

ALLOWED_NEGATED_LIMITATION = (
    "Docking of top-ranked molecules was used as a structure-aware sanity check, "
    "not as proof of binding affinity, therapeutic efficacy, or prospective discovery."
).lower()

BANNED_OVERCLAIMS = [
    "confirmed binding",
    "confirmed inhibitors",
    "therapeutic efficacy",
    "clinical candidate",
    "prospective discovery",
    "binding free energy prediction",
]


def cleaned_public_text(path: Path) -> str:
    text = path.read_text(encoding="utf-8", errors="replace").lower().replace("-", " ")
    return text.replace(ALLOWED_NEGATED_LIMITATION.replace("-", " "), "")


def test_selection_csv_exists_and_has_five_rows():
    assert SELECTION.exists()
    selection = pd.read_csv(SELECTION)
    assert len(selection) == 5
    assert "molecule_id" in selection.columns


def test_sanity_table_exists_and_has_required_columns():
    assert SANITY_TABLE.exists()
    table = pd.read_csv(SANITY_TABLE)
    assert len(table) == 5
    assert REQUIRED_COLUMNS.issubset(set(table.columns))


def test_public_sanity_table_does_not_include_raw_smiles_column():
    table = pd.read_csv(SANITY_TABLE, nrows=0)
    lowered = {column.lower() for column in table.columns}
    assert "smiles" not in lowered
    assert "canonical_smiles" not in lowered
    assert "standardized_smiles" not in lowered


def test_reports_and_primary_figure_exist():
    assert TOP5_REPORT.exists()
    assert FINAL_REPORT.exists()
    assert FIGURE.exists()


def test_final_report_mentions_structure_aware_sanity_check():
    text = FINAL_REPORT.read_text(encoding="utf-8", errors="replace").lower()
    assert "structure-aware sanity check" in text


def test_top5_public_reports_do_not_make_banned_overclaims():
    offenders = []
    for path in [TOP5_REPORT, FINAL_REPORT, ROOT / "README.md", ROOT / "portfolio_assets" / "egfr_project_card.md"]:
        if not path.exists():
            continue
        text = cleaned_public_text(path)
        for phrase in BANNED_OVERCLAIMS:
            if phrase in text:
                offenders.append(f"{path.relative_to(ROOT)}:{phrase}")
    assert not offenders
