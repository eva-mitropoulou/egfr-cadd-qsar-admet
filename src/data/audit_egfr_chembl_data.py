"""Audit EGFR ChEMBL activity provenance and cleaning policy."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    CLEAN_ACTIVITY_PATH,
    METRICS_DIR,
    MODEL_READY_PATH,
    RAW_ACTIVITY_PATH,
    REPORTS_DIR,
    save_json,
    value_counts_dict,
    write_text,
)


def main() -> None:
    """Create ChEMBL activity provenance audit artifacts."""
    if not RAW_ACTIVITY_PATH.exists():
        fetch_script = PROJECT_ROOT / "src" / "fetch_egfr_ic50_raw.py"
        raise FileNotFoundError(f"Raw activity file missing. Rebuild with {fetch_script}")

    raw = pd.read_csv(RAW_ACTIVITY_PATH)
    clean = pd.read_csv(CLEAN_ACTIVITY_PATH) if CLEAN_ACTIVITY_PATH.exists() else pd.DataFrame()
    model_ready = pd.read_csv(MODEL_READY_PATH) if MODEL_READY_PATH.exists() else pd.DataFrame()

    raw_rows = int(len(raw))
    raw_molecules = int(raw["molecule_chembl_id"].nunique()) if "molecule_chembl_id" in raw.columns else None
    assay_count = int(raw["assay_chembl_id"].nunique()) if "assay_chembl_id" in raw.columns else None
    document_count = int(raw["document_chembl_id"].nunique()) if "document_chembl_id" in raw.columns else None
    target_ids = sorted(str(item) for item in raw["target_chembl_id"].dropna().unique()) if "target_chembl_id" in raw.columns else []

    ic50_rows = int((raw.get("standard_type", pd.Series(dtype=str)) == "IC50").sum())
    exact_rows = int((raw.get("standard_relation", pd.Series(dtype=str)) == "=").sum())
    nm_rows = int((raw.get("standard_units", pd.Series(dtype=str)) == "nM").sum())
    duplicate_molecule_rows = raw_rows - raw_molecules if raw_molecules is not None else None

    metrics = {
        "source": "ChEMBL public activity records",
        "target_ids": target_ids,
        "primary_target_id": "CHEMBL203" if "CHEMBL203" in target_ids else None,
        "raw_activity_row_count": raw_rows,
        "raw_unique_molecule_count": raw_molecules,
        "assay_count": assay_count,
        "document_count": document_count,
        "standard_type_distribution": value_counts_dict(raw["standard_type"]) if "standard_type" in raw.columns else {},
        "standard_units_distribution": value_counts_dict(raw["standard_units"]) if "standard_units" in raw.columns else {},
        "standard_relation_distribution": value_counts_dict(raw["standard_relation"]) if "standard_relation" in raw.columns else {},
        "ic50_rows": ic50_rows,
        "exact_relation_rows": exact_rows,
        "nm_rows": nm_rows,
        "duplicate_molecule_row_count": duplicate_molecule_rows,
        "clean_pIC50_molecule_count": int(len(clean)) if not clean.empty else None,
        "model_ready_molecule_count": int(len(model_ready)) if not model_ready.empty else None,
        "pIC50_conversion_policy": "pIC50 = 9 - log10(IC50_nM), using exact nM IC50 values only for the clean regression target.",
        "duplicate_aggregation_policy": "Duplicate activity measurements are aggregated to one molecule-level target using median pIC50 and median IC50_nM.",
    }
    save_json(METRICS_DIR / "egfr_data_provenance_audit.json", metrics)

    lines = [
        "# EGFR ChEMBL Data Provenance Audit",
        "",
        "- Source: public ChEMBL activity records",
        f"- Target IDs found: {', '.join(target_ids) if target_ids else 'not available'}",
        f"- Raw activity rows: {raw_rows:,}",
        f"- Raw unique molecules: {raw_molecules:,}" if raw_molecules is not None else "- Raw unique molecules: not available",
        f"- Assays: {assay_count:,}" if assay_count is not None else "- Assays: not available",
        f"- Documents: {document_count:,}" if document_count is not None else "- Documents: not available",
        f"- IC50 rows: {ic50_rows:,}",
        f"- Exact `=` rows: {exact_rows:,}",
        f"- nM rows: {nm_rows:,}",
        f"- Clean molecule-level pIC50 rows: {metrics['clean_pIC50_molecule_count']:,}",
        f"- Model-ready molecule rows: {metrics['model_ready_molecule_count']:,}",
        "",
        "## Cleaning Policy",
        "",
        f"- {metrics['pIC50_conversion_policy']}",
        f"- {metrics['duplicate_aggregation_policy']}",
        "",
        "## Distribution Summaries",
        "",
        f"- Standard type counts: {metrics['standard_type_distribution']}",
        f"- Standard units counts: {metrics['standard_units_distribution']}",
        f"- Standard relation counts: {metrics['standard_relation_distribution']}",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_data_provenance_audit.md", "\n".join(lines))

    print(f"Raw rows: {raw_rows}")
    print(f"Clean molecule count: {metrics['clean_pIC50_molecule_count']}")
    print(f"Metrics: {METRICS_DIR / 'egfr_data_provenance_audit.json'}")


if __name__ == "__main__":
    main()
