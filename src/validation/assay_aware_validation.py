"""Assay-aware and document-aware validation for the EGFR Morgan RF model."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupShuffleSplit


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    METRICS_DIR,
    PROCESSED_DIR,
    RANDOM_STATE,
    RAW_ACTIVITY_PATH,
    REPORTS_DIR,
    regression_metrics,
    save_json,
    setup_matplotlib,
    write_text,
)


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


REPORT_PATH = REPORTS_DIR / "egfr_assay_aware_validation_report.md"
METRICS_PATH = METRICS_DIR / "egfr_assay_aware_validation_metrics.json"
FIGURE_PATH = REPORTS_DIR / "figures" / "random_vs_scaffold_vs_assay_document_split.png"


def dominant_value(values: pd.Series) -> str | None:
    """Return the most frequent non-missing group label."""
    clean = values.dropna().astype(str)
    if clean.empty:
        return None
    return str(clean.value_counts().idxmax())


def load_molecule_metadata(index: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Trace assay/document metadata from raw activity rows to molecule-level rows."""
    availability = {
        "raw_activity_available": RAW_ACTIVITY_PATH.exists(),
        "assay_metadata_available": False,
        "document_metadata_available": False,
        "metadata_mapping_status": "not_attempted",
    }
    if not RAW_ACTIVITY_PATH.exists():
        availability["metadata_mapping_status"] = "raw_activity_missing"
        return index.copy(), availability

    raw_columns = pd.read_csv(RAW_ACTIVITY_PATH, nrows=0).columns.tolist()
    assay_col = "assay_chembl_id" if "assay_chembl_id" in raw_columns else None
    document_col = None
    for candidate in ("document_chembl_id", "doc_chembl_id", "pubmed_id"):
        if candidate in raw_columns:
            document_col = candidate
            break
    read_columns = ["molecule_chembl_id"]
    if assay_col:
        read_columns.append(assay_col)
    if document_col:
        read_columns.append(document_col)
    raw = pd.read_csv(RAW_ACTIVITY_PATH, usecols=read_columns)
    grouped = raw.groupby("molecule_chembl_id", dropna=False).agg(
        dominant_assay_chembl_id=(assay_col, dominant_value) if assay_col else ("molecule_chembl_id", lambda _: None),
        dominant_document_chembl_id=(document_col, dominant_value) if document_col else ("molecule_chembl_id", lambda _: None),
        raw_activity_record_count=("molecule_chembl_id", "size"),
    )
    metadata = index.merge(grouped.reset_index(), on="molecule_chembl_id", how="left")
    availability.update(
        {
            "assay_metadata_available": bool(assay_col),
            "document_metadata_available": bool(document_col),
            "assay_column": assay_col,
            "document_column": document_col,
            "metadata_mapping_status": "mapped_from_raw_activity",
            "molecules_with_assay_metadata": int(metadata["dominant_assay_chembl_id"].notna().sum()),
            "molecules_with_document_metadata": int(metadata["dominant_document_chembl_id"].notna().sum()),
            "unique_assay_groups": int(metadata["dominant_assay_chembl_id"].nunique(dropna=True)),
            "unique_document_groups": int(metadata["dominant_document_chembl_id"].nunique(dropna=True)),
        }
    )
    return metadata, availability


def evaluate_group_split(x_matrix, metadata: pd.DataFrame, group_column: str, label: str) -> dict:
    """Train/evaluate a Morgan RF with a group-disjoint split."""
    usable = metadata[group_column].notna().to_numpy()
    row_indices = np.flatnonzero(usable)
    groups = metadata.loc[usable, group_column].astype(str).to_numpy()
    unique_groups = pd.Series(groups).nunique()
    if unique_groups < 5 or len(row_indices) < 100:
        return {
            "split": label,
            "status": "degraded_group_metadata_insufficient",
            "usable_rows": int(len(row_indices)),
            "unique_groups": int(unique_groups),
        }

    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_STATE)
    local_train, local_test = next(splitter.split(row_indices, groups=groups))
    train_idx = row_indices[local_train]
    test_idx = row_indices[local_test]
    train_groups = set(metadata.iloc[train_idx][group_column].astype(str))
    test_groups = set(metadata.iloc[test_idx][group_column].astype(str))
    train_ids = set(metadata.iloc[train_idx]["molecule_chembl_id"].astype(str))
    test_ids = set(metadata.iloc[test_idx]["molecule_chembl_id"].astype(str))

    model = RandomForestRegressor(n_estimators=100, random_state=RANDOM_STATE, n_jobs=-1)
    y = metadata["median_pIC50"].astype(float).to_numpy()
    model.fit(x_matrix[train_idx], y[train_idx])
    predictions = model.predict(x_matrix[test_idx])
    metrics = regression_metrics(y[test_idx], predictions)
    return {
        "split": label,
        "status": "completed",
        "model": "Random Forest Morgan",
        "train_size": int(len(train_idx)),
        "test_size": int(len(test_idx)),
        "train_unique_groups": int(len(train_groups)),
        "test_unique_groups": int(len(test_groups)),
        "group_overlap_count": int(len(train_groups.intersection(test_groups))),
        "molecule_overlap_count": int(len(train_ids.intersection(test_ids))),
        **metrics,
    }


def main() -> None:
    """Run assay/document-aware validation and write reports."""
    x_matrix = sparse.load_npz(PROCESSED_DIR / "features_morgan_fingerprints.npz").astype(np.float32)
    index = pd.read_csv(PROCESSED_DIR / "features_morgan_index.csv")
    metadata, availability = load_molecule_metadata(index)

    rows: list[dict] = []
    rows.append(evaluate_group_split(x_matrix, metadata, "dominant_assay_chembl_id", "assay_group_split"))
    rows.append(evaluate_group_split(x_matrix, metadata, "dominant_document_chembl_id", "document_group_split"))

    benchmark_path = METRICS_DIR / "qsar_matched_benchmark_metrics.json"
    reference = {}
    if benchmark_path.exists():
        import json

        reference = json.loads(benchmark_path.read_text(encoding="utf-8"))
    reference_rows = []
    for row in reference.get("matched_benchmark_rows", []):
        if row.get("feature_set") == "morgan_fingerprints" and row.get("model") == "Random Forest":
            reference_rows.append(row)

    metrics_payload = {
        "status": "completed" if any(row.get("status") == "completed" for row in rows) else "degraded_no_group_validation",
        "metadata_availability": availability,
        "validation_rows": rows,
        "reference_rows": reference_rows,
    }
    save_json(METRICS_PATH, metrics_payload)

    completed_rows = [row for row in rows if row.get("status") == "completed"]
    plot_rows = [
        {"split": row["split"], "RMSE": row.get("RMSE"), "R2": row.get("R2")}
        for row in reference_rows + completed_rows
        if row.get("RMSE") is not None
    ]
    if plot_rows:
        plot_df = pd.DataFrame(plot_rows)
        FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        axes[0].bar(plot_df["split"], plot_df["RMSE"], color="#4C78A8")
        axes[0].set_ylabel("RMSE")
        axes[0].tick_params(axis="x", rotation=25)
        axes[1].bar(plot_df["split"], plot_df["R2"], color="#59A14F")
        axes[1].set_ylabel("R2")
        axes[1].tick_params(axis="x", rotation=25)
        fig.suptitle("Morgan RF Validation Contexts")
        fig.tight_layout()
        fig.savefig(FIGURE_PATH, dpi=200)
        plt.close(fig)

    lines = [
        "# EGFR Assay-Aware and Document-Aware Validation",
        "",
        "This stage tests whether the Morgan fingerprint Random Forest remains useful when held-out groups are defined by assay or document context.",
        "",
        "## Metadata Availability",
        "",
        f"- Raw activity table available: {availability['raw_activity_available']}",
        f"- Assay metadata available: {availability['assay_metadata_available']}",
        f"- Document metadata available: {availability['document_metadata_available']}",
        f"- Metadata mapping status: {availability['metadata_mapping_status']}",
        f"- Molecules with assay metadata: {availability.get('molecules_with_assay_metadata')}",
        f"- Molecules with document metadata: {availability.get('molecules_with_document_metadata')}",
        "",
        "## Validation Summary",
        "",
    ]
    for row in rows:
        lines.extend(
            [
                f"### {row['split']}",
                "",
                f"- Status: {row.get('status')}",
                f"- Train molecules: {row.get('train_size')}",
                f"- Test molecules: {row.get('test_size')}",
                f"- Train unique groups: {row.get('train_unique_groups')}",
                f"- Test unique groups: {row.get('test_unique_groups')}",
                f"- Group overlap count: {row.get('group_overlap_count')}",
                f"- Molecule overlap count: {row.get('molecule_overlap_count')}",
                f"- MAE: {row.get('MAE')}",
                f"- RMSE: {row.get('RMSE')}",
                f"- R2: {row.get('R2')}",
                f"- Pearson: {row.get('Pearson')}",
                f"- Spearman: {row.get('Spearman')}",
                "",
            ]
        )
    lines.append("No raw molecule structures or raw assay records are shown in this report.")
    write_text(REPORT_PATH, "\n".join(lines))

    print(f"Assay-aware validation status: {metrics_payload['status']}")
    print(f"Assay metadata available: {availability['assay_metadata_available']}")
    print(f"Document metadata available: {availability['document_metadata_available']}")


if __name__ == "__main__":
    main()
