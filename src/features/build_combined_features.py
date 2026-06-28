"""Build combined descriptor plus Morgan fingerprint feature matrix."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from scipy import sparse
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    DESCRIPTOR_COLUMNS,
    METRICS_DIR,
    PROCESSED_DIR,
    REPORTS_DIR,
    TARGET_COLUMN,
    read_json,
    save_json,
)


DESCRIPTOR_PATH = PROCESSED_DIR / "features_rdkit_descriptors.csv"
FP_PATH = PROCESSED_DIR / "features_morgan_fingerprints.npz"
INDEX_PATH = PROCESSED_DIR / "features_morgan_index.csv"
COMBINED_PATH = PROCESSED_DIR / "features_combined_descriptors_morgan.npz"
COMBINED_INDEX_PATH = PROCESSED_DIR / "features_combined_index.csv"


def main() -> None:
    """Create combined sparse matrix and metadata report."""
    if not DESCRIPTOR_PATH.exists() or not FP_PATH.exists() or not INDEX_PATH.exists():
        raise FileNotFoundError("Descriptor and Morgan feature files must be built before combined features.")

    descriptors = pd.read_csv(DESCRIPTOR_PATH)
    fp_matrix = sparse.load_npz(FP_PATH)
    index = pd.read_csv(INDEX_PATH)

    if len(descriptors) != fp_matrix.shape[0] or len(index) != fp_matrix.shape[0]:
        raise ValueError("Feature matrices are not aligned by row count.")

    x_desc = descriptors[DESCRIPTOR_COLUMNS].astype(float).to_numpy()
    x_desc_scaled = StandardScaler().fit_transform(x_desc)
    combined = sparse.hstack([sparse.csr_matrix(x_desc_scaled), fp_matrix], format="csr")
    sparse.save_npz(COMBINED_PATH, combined)
    index.to_csv(COMBINED_INDEX_PATH, index=False)

    metrics = {
        "combined_rows": int(combined.shape[0]),
        "combined_columns": int(combined.shape[1]),
        "descriptor_columns": int(len(DESCRIPTOR_COLUMNS)),
        "fingerprint_columns": int(fp_matrix.shape[1]),
        "feature_label_alignment": int(combined.shape[0]) == int(index[TARGET_COLUMN].notna().sum()),
        "combined_path": str(COMBINED_PATH.relative_to(PROJECT_ROOT)),
        "combined_index_path": str(COMBINED_INDEX_PATH.relative_to(PROJECT_ROOT)),
    }
    metrics_path = METRICS_DIR / "feature_generation_metrics.json"
    existing = read_json(metrics_path)
    existing["combined_features"] = metrics
    save_json(metrics_path, existing)

    with (REPORTS_DIR / "feature_generation_report.md").open("a", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "",
                    "## Combined Features",
                    "",
                    f"- Combined matrix shape: {combined.shape[0]:,} x {combined.shape[1]:,}",
                    f"- Descriptor columns: {len(DESCRIPTOR_COLUMNS)}",
                    f"- Morgan fingerprint columns: {fp_matrix.shape[1]:,}",
                    f"- Feature-label alignment: {metrics['feature_label_alignment']}",
                    "",
                ]
            )
        )

    print(f"Combined matrix rows: {combined.shape[0]}")
    print(f"Combined columns: {combined.shape[1]}")


if __name__ == "__main__":
    main()
