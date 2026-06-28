"""Inventory artifacts and metadata availability for EGFR hardening."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, save_json, write_text  # noqa: E402


REPORT_PATH = REPORTS_DIR / "egfr_hardening_inventory.md"
METRICS_PATH = METRICS_DIR / "egfr_hardening_inventory.json"


FILES_TO_AUDIT = [
    "data/processed/egfr_model_ready.csv",
    "data/processed/egfr_ic50_clean.csv",
    "data/processed/egfr_standardized_molecules.csv",
    "data/processed/features_morgan_index.csv",
    "data/processed/features_morgan_fingerprints.npz",
    "data/processed/features_rdkit_descriptors.csv",
    "data/raw/egfr_chembl_ic50_raw.csv",
    "reports/metrics/qsar_matched_benchmark_metrics.json",
    "reports/metrics/applicability_domain_metrics.json",
    "reports/metrics/egfr_gnn_benchmark_metrics.json",
    "reports/metrics/egfr_redocking_metrics.json",
    "reports/egfr_ranked_existing_molecules.csv",
    "data/structure_prepared/5UG9_receptor.pdbqt",
    "data/structure_prepared/5UG9_8AM_ligand.pdbqt",
    "data/structure_prepared/5UG9_8AM_redocked_out.pdbqt",
]


def csv_shape_and_columns(path: Path) -> dict:
    """Return row count, columns, and missing counts for CSV."""
    columns = pd.read_csv(path, nrows=0).columns.tolist()
    df = pd.read_csv(path)
    missing = {column: int(df[column].isna().sum()) for column in columns}
    return {"shape": [int(df.shape[0]), int(df.shape[1])], "columns": columns, "missing_counts": missing}


def json_keys(path: Path) -> dict:
    """Return top-level JSON keys."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    return {"json_keys": sorted(payload.keys())}


def audit_file(relative_path: str) -> dict:
    """Audit one file without exposing records."""
    path = PROJECT_ROOT / relative_path
    item = {
        "path": relative_path,
        "exists": path.exists(),
        "bytes": path.stat().st_size if path.exists() else 0,
    }
    if not path.exists():
        return item
    if path.suffix == ".csv":
        item.update(csv_shape_and_columns(path))
    elif path.suffix == ".json":
        item.update(json_keys(path))
    elif path.suffix == ".npz":
        from scipy import sparse

        matrix = sparse.load_npz(path)
        item["shape"] = [int(matrix.shape[0]), int(matrix.shape[1])]
        item["nnz"] = int(matrix.nnz)
    return item


def main() -> None:
    """Create hardening artifact inventory."""
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    METRICS_DIR.mkdir(parents=True, exist_ok=True)
    audited = [audit_file(path) for path in FILES_TO_AUDIT]
    model_ready = next(item for item in audited if item["path"] == "data/processed/egfr_model_ready.csv")
    raw = next(item for item in audited if item["path"] == "data/raw/egfr_chembl_ic50_raw.csv")
    model_columns = set(model_ready.get("columns", []))
    raw_columns = set(raw.get("columns", []))
    metadata = {
        "model_ready_found": bool(model_ready["exists"]),
        "pIC50_label_found": "median_pIC50" in model_columns,
        "molecule_identifier_found": "molecule_chembl_id" in model_columns,
        "molecule_representation_found": bool({"canonical_smiles", "standardized_smiles"}.intersection(model_columns)),
        "assay_metadata_available": "assay_chembl_id" in raw_columns,
        "document_metadata_available": bool({"document_chembl_id", "doc_chembl_id", "pubmed_id"}.intersection(raw_columns)),
        "assay_metadata_source": "data/raw/egfr_chembl_ic50_raw.csv" if "assay_chembl_id" in raw_columns else None,
        "document_metadata_source": "data/raw/egfr_chembl_ic50_raw.csv" if raw_columns else None,
    }
    payload = {"status": "completed", "files": audited, "metadata_availability": metadata}
    save_json(METRICS_PATH, payload)

    lines = [
        "# EGFR Hardening Inventory",
        "",
        "This inventory lists artifact paths, shapes, columns, missing-value counts, and metadata availability without exposing raw molecule structures.",
        "",
        "## Metadata Availability",
        "",
    ]
    for key, value in metadata.items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Artifact Summary", ""])
    for item in audited:
        lines.append(f"- `{item['path']}`: exists={item['exists']}, bytes={item['bytes']}, shape={item.get('shape')}")
        if "columns" in item:
            lines.append(f"  columns: {', '.join(item['columns'])}")
    write_text(REPORT_PATH, "\n".join(lines) + "\n")

    print("Hardening inventory status: completed")
    print(f"Audited files: {len(audited)}")


if __name__ == "__main__":
    main()
