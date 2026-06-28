"""Applicability-domain analysis using existing scaffold-CV prediction evidence."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import FIGURES_DIR, METRICS_DIR, REPORTS_DIR, markdown_table, save_figure, save_json, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


SOURCE_PREDICTIONS = PROJECT_ROOT / "results" / "applicability_domain_predictions.csv"
SOURCE_SUMMARY = PROJECT_ROOT / "results" / "applicability_domain_summary.csv"
PREDICTIONS_PATH = REPORTS_DIR / "egfr_applicability_domain_predictions.csv"
METRICS_PATH = METRICS_DIR / "applicability_domain_metrics.json"
REPORT_PATH = REPORTS_DIR / "applicability_domain_report.md"


def similarity_bin(value: float) -> str:
    """Assign low/medium/high applicability-domain bins."""
    if value < 0.3:
        return "low"
    if value <= 0.7:
        return "medium"
    return "high"


def main() -> None:
    """Create applicability-domain report from existing OOF scaffold-CV predictions."""
    if not SOURCE_PREDICTIONS.exists():
        raise FileNotFoundError(f"Missing existing applicability-domain predictions: {SOURCE_PREDICTIONS}")

    source = pd.read_csv(SOURCE_PREDICTIONS)
    required = {"molecule_chembl_id", "true_pIC50", "predicted_pIC50", "absolute_error", "max_tanimoto_to_train"}
    if not required.issubset(source.columns):
        raise ValueError("Existing applicability-domain prediction file lacks required columns.")

    safe = source[[column for column in source.columns if column in required.union({"fold", "scaffold"})]].copy()
    safe["molecule_hash"] = safe["molecule_chembl_id"].astype(str)
    if "scaffold" in safe.columns:
        import hashlib

        safe["scaffold_hash"] = safe["scaffold"].astype(str).map(lambda value: hashlib.sha256(value.encode("utf-8")).hexdigest()[:12])
        safe = safe.drop(columns=["scaffold"])
    else:
        safe["scaffold_hash"] = "not_available"
    safe["similarity_bin"] = safe["max_tanimoto_to_train"].map(similarity_bin)
    safe["out_of_domain_flag"] = safe["similarity_bin"] == "low"
    safe.to_csv(PREDICTIONS_PATH, index=False)

    summary = (
        safe.groupby("similarity_bin")
        .agg(
            count=("absolute_error", "size"),
            MAE=("absolute_error", "mean"),
            RMSE=("absolute_error", lambda values: float(np.sqrt(np.mean(np.square(values))))),
            mean_similarity=("max_tanimoto_to_train", "mean"),
        )
        .reindex(["low", "medium", "high"])
        .reset_index()
    )
    low_mae = float(summary.loc[summary["similarity_bin"] == "low", "MAE"].iloc[0])
    high_mae = float(summary.loc[summary["similarity_bin"] == "high", "MAE"].iloc[0])

    metrics = {
        "prediction_rows": int(len(safe)),
        "summary_by_similarity_bin": summary.to_dict(orient="records"),
        "low_similarity_mae": low_mae,
        "high_similarity_mae": high_mae,
        "low_vs_high_mae_difference": float(low_mae - high_mae),
        "out_of_domain_count": int(safe["out_of_domain_flag"].sum()),
        "evidence_source": "existing scaffold-CV applicability-domain predictions",
        "confidence_rules": {
            "low": "max Tanimoto < 0.3; out-of-domain warning",
            "medium": "0.3 <= max Tanimoto <= 0.7; moderate reliability",
            "high": "max Tanimoto > 0.7; interpolation-like prediction",
        },
    }
    save_json(METRICS_PATH, metrics)

    plt.figure(figsize=(7, 4.5))
    plt.scatter(safe["max_tanimoto_to_train"], safe["absolute_error"], s=10, alpha=0.25)
    plt.xlabel("Max Tanimoto similarity to training set")
    plt.ylabel("Absolute pIC50 error")
    plt.title("Applicability Domain: Error vs Similarity")
    save_figure(FIGURES_DIR / "error_vs_similarity.png")

    plt.figure(figsize=(7, 4.5))
    plt.hist(safe["max_tanimoto_to_train"], bins=30, color="#4C78A8", alpha=0.85)
    plt.xlabel("Max Tanimoto similarity to training set")
    plt.ylabel("Molecule count")
    plt.title("Similarity Distribution")
    save_figure(FIGURES_DIR / "similarity_distribution.png")

    report = [
        "# Applicability Domain Report",
        "",
        "Applicability domain was estimated from max Tanimoto similarity to the training set during scaffold GroupKFold validation.",
        "",
        markdown_table(summary),
        "",
        f"- Low-similarity MAE: {low_mae:.3f}",
        f"- High-similarity MAE: {high_mae:.3f}",
        f"- Out-of-domain molecules: {metrics['out_of_domain_count']:,}",
        "",
        "Predictions with max Tanimoto < 0.3 are flagged as out-of-domain.",
        "",
    ]
    write_text(REPORT_PATH, "\n".join(report))

    print(f"Applicability rows: {len(safe)}")
    print(f"Out-of-domain count: {metrics['out_of_domain_count']}")


if __name__ == "__main__":
    main()
