"""Conformal-style uncertainty checks for EGFR Morgan RF QSAR."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    METRICS_DIR,
    PROCESSED_DIR,
    RANDOM_STATE,
    REPORTS_DIR,
    regression_metrics,
    save_json,
    setup_matplotlib,
    write_text,
)


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


REPORT_PATH = REPORTS_DIR / "egfr_conformal_uncertainty_report.md"
METRICS_PATH = METRICS_DIR / "egfr_conformal_uncertainty_metrics.json"
FIGURES_DIR = REPORTS_DIR / "figures"
TARGET_COVERAGE = 0.90


def conformal_quantile(abs_residuals: np.ndarray, target_coverage: float = TARGET_COVERAGE) -> float:
    """Return residual quantile with finite-sample correction."""
    n_cal = len(abs_residuals)
    if n_cal == 0:
        raise ValueError("Calibration residuals are empty.")
    rank = math.ceil((n_cal + 1) * target_coverage)
    q = min(rank / n_cal, 1.0)
    return float(np.quantile(abs_residuals, q, method="higher"))


def similarity_bin(value: float) -> str:
    """Assign applicability-domain bin labels."""
    if value < 0.3:
        return "low"
    if value <= 0.7:
        return "medium"
    return "high"


def max_tanimoto_sparse(x_test, x_train, batch_size: int = 256) -> np.ndarray:
    """Calculate max Tanimoto similarity from binary sparse rows to a train matrix."""
    train_counts = np.asarray(x_train.sum(axis=1)).ravel().astype(float)
    results: list[np.ndarray] = []
    for start in range(0, x_test.shape[0], batch_size):
        batch = x_test[start : start + batch_size]
        batch_counts = np.asarray(batch.sum(axis=1)).ravel().astype(float)
        intersections = batch @ x_train.T
        intersections = intersections.astype(float).toarray()
        denominators = batch_counts[:, None] + train_counts[None, :] - intersections
        with np.errstate(divide="ignore", invalid="ignore"):
            similarities = np.divide(intersections, denominators, out=np.zeros_like(intersections), where=denominators > 0)
        results.append(similarities.max(axis=1))
    return np.concatenate(results)


def scaffold_three_way_split(index: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create train/calibration/test indices by whole scaffold groups."""
    rng = np.random.default_rng(RANDOM_STATE)
    groups = index.groupby("scaffold_hash", sort=False).indices
    group_items = [(group, np.asarray(rows, dtype=int)) for group, rows in groups.items()]
    rng.shuffle(group_items)
    n_total = len(index)
    train_target = 0.60 * n_total
    cal_target = 0.20 * n_total
    train: list[int] = []
    cal: list[int] = []
    test: list[int] = []
    for _, rows in group_items:
        if len(train) < train_target:
            train.extend(rows.tolist())
        elif len(cal) < cal_target:
            cal.extend(rows.tolist())
        else:
            test.extend(rows.tolist())
    return np.asarray(train, dtype=int), np.asarray(cal, dtype=int), np.asarray(test, dtype=int)


def evaluate_split(name: str, x_matrix, y: np.ndarray, train_idx: np.ndarray, cal_idx: np.ndarray, test_idx: np.ndarray) -> tuple[dict, pd.DataFrame]:
    """Fit model, calibrate conformal interval, and evaluate held-out coverage."""
    model = RandomForestRegressor(n_estimators=100, random_state=RANDOM_STATE, n_jobs=-1)
    model.fit(x_matrix[train_idx], y[train_idx])
    cal_pred = model.predict(x_matrix[cal_idx])
    cal_abs = np.abs(y[cal_idx] - cal_pred)
    q_hat = conformal_quantile(cal_abs)
    test_pred = model.predict(x_matrix[test_idx])
    lower = test_pred - q_hat
    upper = test_pred + q_hat
    covered = (y[test_idx] >= lower) & (y[test_idx] <= upper)
    max_sim = max_tanimoto_sparse(x_matrix[test_idx], x_matrix[train_idx])
    errors = np.abs(y[test_idx] - test_pred)
    pred_df = pd.DataFrame(
        {
            "split": name,
            "row_index": test_idx,
            "true_pIC50": y[test_idx],
            "predicted_pIC50": test_pred,
            "absolute_error": errors,
            "interval_lower_90": lower,
            "interval_upper_90": upper,
            "interval_width": upper - lower,
            "covered_90": covered,
            "max_tanimoto_to_train": max_sim,
            "similarity_bin": [similarity_bin(value) for value in max_sim],
        }
    )
    bin_rows = (
        pred_df.groupby("similarity_bin", observed=True)
        .agg(
            count=("covered_90", "size"),
            empirical_coverage=("covered_90", "mean"),
            mean_interval_width=("interval_width", "mean"),
            MAE=("absolute_error", "mean"),
        )
        .reset_index()
        .to_dict(orient="records")
    )
    metrics = {
        "split": name,
        "status": "completed",
        "model": "Random Forest Morgan",
        "target_coverage": TARGET_COVERAGE,
        "train_size": int(len(train_idx)),
        "calibration_size": int(len(cal_idx)),
        "test_size": int(len(test_idx)),
        "conformal_q_hat": float(q_hat),
        "empirical_coverage": float(covered.mean()),
        "mean_interval_width": float((upper - lower).mean()),
        "median_interval_width": float(np.median(upper - lower)),
        "error_interval_width_spearman": float(pred_df["absolute_error"].corr(pred_df["interval_width"], method="spearman")),
        "coverage_by_similarity_bin": bin_rows,
        **regression_metrics(y[test_idx], test_pred),
    }
    return metrics, pred_df


def save_figures(metrics_rows: list[dict], predictions: pd.DataFrame) -> None:
    """Save conformal uncertainty figures."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    metrics_df = pd.DataFrame(metrics_rows)
    plt.figure(figsize=(7, 4.5))
    plt.bar(metrics_df["split"], metrics_df["empirical_coverage"], color="#4C78A8")
    plt.axhline(TARGET_COVERAGE, color="black", linestyle="--", linewidth=1)
    plt.ylabel("Empirical coverage")
    plt.title("90% Conformal-Style Coverage")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "conformal_coverage_by_split.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.scatter(predictions["max_tanimoto_to_train"], predictions["interval_width"], s=10, alpha=0.25)
    plt.xlabel("Max Tanimoto to training set")
    plt.ylabel("90% interval width")
    plt.title("Conformal Interval Width vs Similarity")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "conformal_interval_width_vs_similarity.png", dpi=200)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    for split, group in predictions.groupby("split"):
        plt.hist(group["interval_width"], bins=20, alpha=0.5, label=split)
    plt.xlabel("90% interval width")
    plt.ylabel("Molecule count")
    plt.title("Conformal Interval Width Distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "conformal_interval_width_distribution.png", dpi=200)
    plt.close()


def main() -> None:
    """Run random and scaffold conformal-style uncertainty analysis."""
    x_matrix = sparse.load_npz(PROCESSED_DIR / "features_morgan_fingerprints.npz").astype(np.float32)
    index = pd.read_csv(PROCESSED_DIR / "features_morgan_index.csv")
    y = index["median_pIC50"].astype(float).to_numpy()
    all_indices = np.arange(len(index))

    train_cal, random_test = train_test_split(all_indices, test_size=0.20, random_state=RANDOM_STATE)
    random_train, random_cal = train_test_split(train_cal, test_size=0.25, random_state=RANDOM_STATE)
    scaffold_train, scaffold_cal, scaffold_test = scaffold_three_way_split(index)

    metrics_rows: list[dict] = []
    prediction_frames: list[pd.DataFrame] = []
    for split_name, train_idx, cal_idx, test_idx in [
        ("random_split_conformal", random_train, random_cal, random_test),
        ("scaffold_group_conformal", scaffold_train, scaffold_cal, scaffold_test),
    ]:
        metrics, pred_df = evaluate_split(split_name, x_matrix, y, train_idx, cal_idx, test_idx)
        metrics_rows.append(metrics)
        prediction_frames.append(pred_df)

    predictions = pd.concat(prediction_frames, ignore_index=True)
    save_figures(metrics_rows, predictions)

    payload = {
        "status": "completed",
        "target_coverage": TARGET_COVERAGE,
        "method": "split_conformal_absolute_residual_intervals",
        "splits": metrics_rows,
    }
    save_json(METRICS_PATH, payload)

    lines = [
        "# EGFR Conformal-Style Uncertainty Report",
        "",
        "This stage adds conformal-style pIC50 uncertainty checks using absolute calibration residuals. These intervals are retrospective uncertainty summaries.",
        "",
    ]
    for row in metrics_rows:
        lines.extend(
            [
                f"## {row['split']}",
                "",
                f"- Status: {row['status']}",
                f"- Target coverage: {row['target_coverage']:.2f}",
                f"- Empirical coverage: {row['empirical_coverage']:.3f}",
                f"- Mean interval width: {row['mean_interval_width']:.3f}",
                f"- Median interval width: {row['median_interval_width']:.3f}",
                f"- Test MAE: {row['MAE']:.3f}",
                f"- Test RMSE: {row['RMSE']:.3f}",
                f"- Test R2: {row['R2']:.3f}",
                "",
            ]
        )
    lines.append("Intervals quantify retrospective model uncertainty on held-out EGFR records; they are not clinical confidence statements.")
    write_text(REPORT_PATH, "\n".join(lines))

    print("Conformal uncertainty status: completed")
    for row in metrics_rows:
        print(f"{row['split']} coverage: {row['empirical_coverage']:.3f}")


if __name__ == "__main__":
    main()
