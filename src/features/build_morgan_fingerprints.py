"""Build Morgan fingerprint feature matrix for standardized EGFR molecules."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    FINGERPRINT_BITS,
    FINGERPRINT_RADIUS,
    METRICS_DIR,
    PROCESSED_DIR,
    REPORTS_DIR,
    TARGET_COLUMN,
    calculate_morgan_matrix,
    load_standardized_or_model_ready,
    read_json,
    save_json,
)


FP_PATH = PROCESSED_DIR / "features_morgan_fingerprints.npz"
INDEX_PATH = PROCESSED_DIR / "features_morgan_index.csv"


def main() -> None:
    """Create sparse Morgan fingerprint matrix and index metadata."""
    df = load_standardized_or_model_ready()
    smiles_column = "standardized_smiles" if "standardized_smiles" in df.columns else "canonical_smiles"

    matrix, _bit_vectors, invalid_count = calculate_morgan_matrix(df[smiles_column])
    sparse.save_npz(FP_PATH, matrix)

    index_columns = ["molecule_chembl_id", TARGET_COLUMN]
    for optional in ["molecule_hash", "scaffold_hash", "scaffold"]:
        if optional in df.columns:
            index_columns.append(optional)
    index = df[index_columns].copy()
    index.to_csv(INDEX_PATH, index=False)

    bit_sums = matrix.sum(axis=0).A1
    metrics = {
        "input_rows": int(len(df)),
        "fingerprint_rows": int(matrix.shape[0]),
        "fingerprint_bits": int(matrix.shape[1]),
        "radius": FINGERPRINT_RADIUS,
        "invalid_molecule_count": int(invalid_count),
        "nonzero_bit_count": int((bit_sums > 0).sum()),
        "constant_zero_bit_count": int((bit_sums == 0).sum()),
        "feature_label_alignment": int(matrix.shape[0]) == int(index[TARGET_COLUMN].notna().sum()),
        "fingerprint_path": str(FP_PATH.relative_to(PROJECT_ROOT)),
        "index_path": str(INDEX_PATH.relative_to(PROJECT_ROOT)),
    }

    metrics_path = METRICS_DIR / "feature_generation_metrics.json"
    existing = read_json(metrics_path)
    existing["morgan_fingerprints"] = metrics
    save_json(metrics_path, existing)

    with (REPORTS_DIR / "feature_generation_report.md").open("a", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "",
                    "## Morgan Fingerprints",
                    "",
                    f"- Fingerprint matrix shape: {matrix.shape[0]:,} x {matrix.shape[1]:,}",
                    f"- Radius: {FINGERPRINT_RADIUS}",
                    f"- Nonzero bit count: {metrics['nonzero_bit_count']:,}",
                    f"- Constant zero bit count: {metrics['constant_zero_bit_count']:,}",
                    f"- Feature-label alignment: {metrics['feature_label_alignment']}",
                    "",
                ]
            )
        )

    print(f"Fingerprint matrix rows: {matrix.shape[0]}")
    print(f"Fingerprint bits: {matrix.shape[1]}")


if __name__ == "__main__":
    main()
