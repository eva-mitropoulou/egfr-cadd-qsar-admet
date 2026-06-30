import re
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
MODEL_READY = ROOT / "data" / "processed" / "egfr_model_ready.csv"
ANNOTATED = ROOT / "data" / "processed" / "egfr_model_ready_with_medchem_alerts.csv"
RANKED = ROOT / "reports" / "egfr_ranked_existing_molecules.csv"
SENSITIVITY_REPORT = ROOT / "reports" / "egfr_medchem_alert_sensitivity_report.md"
FINAL_REPORT = ROOT / "reports" / "final_egfr_cadd_qsar_report.md"

REQUIRED_WORDING = (
    "PAINS, Brenk, and external unwanted-substructure SMARTS alerts were used as "
    "medicinal-chemistry risk annotations and sensitivity-analysis filters, not "
    "automatic exclusions from the primary EGFR QSAR benchmark."
)

ALERT_COLUMNS = {
    "pains_flag",
    "pains_alert_count",
    "brenk_flag",
    "brenk_alert_count",
    "unwanted_substructure_flag",
    "unwanted_substructure_count",
    "medchem_alert_flag",
    "medchem_alert_count",
    "medchem_alert_summary",
}

RANKED_ALERT_COLUMNS = ALERT_COLUMNS | {"triage_risk_bin", "triage_reason"}
SMILES_LIKE = re.compile(r"(?<![A-Za-z0-9_])[BCNOFPSIclbr@+\-\[\]\(\)=#$\\/0-9]{30,}(?![A-Za-z0-9_])")
BANNED_PHRASES = [
    "pains are false positives",
    "removed all bad molecules",
    "therapeutic efficacy",
    "clinical candidate",
    "production grade",
    "prospective discovery",
]


def require_local_processed_data() -> None:
    if not MODEL_READY.exists():
        pytest.skip("Local processed EGFR tables are not present in this checkout.")


def test_annotated_table_exists_and_preserves_row_count():
    require_local_processed_data()
    assert ANNOTATED.exists()
    model_ready_rows = len(pd.read_csv(MODEL_READY, usecols=["molecule_chembl_id"]))
    annotated = pd.read_csv(ANNOTATED)
    assert len(annotated) == model_ready_rows


def test_alert_columns_exist():
    require_local_processed_data()
    annotated = pd.read_csv(ANNOTATED, nrows=0)
    assert ALERT_COLUMNS.issubset(set(annotated.columns))


def test_ranked_table_contains_alert_risk_columns():
    require_local_processed_data()
    assert RANKED.exists()
    ranked = pd.read_csv(RANKED, nrows=0)
    assert RANKED_ALERT_COLUMNS.issubset(set(ranked.columns))


def test_sensitivity_report_exists():
    require_local_processed_data()
    assert SENSITIVITY_REPORT.exists()


def test_final_report_defines_alerts_as_annotations_not_exclusions():
    require_local_processed_data()
    text = FINAL_REPORT.read_text(encoding="utf-8")
    assert REQUIRED_WORDING in text


def test_public_markdown_reports_do_not_expose_obvious_raw_smiles():
    public_paths = list((ROOT / "reports").glob("*.md")) + list((ROOT / "portfolio_assets").glob("*.md")) + [ROOT / "README.md"]
    offenders = []
    for path in public_paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if SMILES_LIKE.search(text):
            offenders.append(str(path.relative_to(ROOT)))
    assert not offenders


def test_public_markdown_reports_do_not_overclaim_medchem_alerts():
    public_paths = list((ROOT / "reports").glob("*.md")) + list((ROOT / "portfolio_assets").glob("*.md")) + [ROOT / "README.md"]
    offenders = []
    for path in public_paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace").lower().replace("-", " ")
        for phrase in BANNED_PHRASES:
            if phrase in text:
                offenders.append(f"{path.relative_to(ROOT)}:{phrase}")
    assert not offenders
