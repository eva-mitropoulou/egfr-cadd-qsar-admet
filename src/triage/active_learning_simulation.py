"""Fast retrospective active-learning simulation over existing EGFR records."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import pairwise_distances
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import DESCRIPTOR_COLUMNS, FIGURES_DIR, METRICS_DIR, RANDOM_STATE, REPORTS_DIR, TARGET_COLUMN, load_standardized_or_model_ready, markdown_table, save_figure, save_json, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


METRICS_PATH = METRICS_DIR / "egfr_active_learning_metrics.json"
REPORT_PATH = REPORTS_DIR / "egfr_active_learning_report.md"


def select_diverse(indices: np.ndarray, scores: np.ndarray, scaffolds: np.ndarray, batch_size: int) -> np.ndarray:
    """Select high-scoring molecules while preserving scaffold diversity."""
    order = indices[np.argsort(scores)[::-1]]
    selected = []
    seen = set()
    for idx in order:
        scaffold = scaffolds[idx]
        if scaffold in seen:
            continue
        selected.append(idx)
        seen.add(scaffold)
        if len(selected) >= batch_size:
            break
    if len(selected) < batch_size:
        for idx in order:
            if idx not in selected:
                selected.append(idx)
            if len(selected) >= batch_size:
                break
    return np.array(selected, dtype=int)


def distance_uncertainty(x: np.ndarray, pool_indices: np.ndarray, selected_indices: np.ndarray) -> np.ndarray:
    """Use nearest selected descriptor distance as an uncertainty/diversity proxy."""
    distances = pairwise_distances(x[pool_indices], x[selected_indices], metric="euclidean")
    return distances.min(axis=1)


def run_strategy(strategy: str, x: np.ndarray, y: np.ndarray, scaffolds: np.ndarray, seed_indices: np.ndarray, potent_indices: set[int], rng: np.random.Generator) -> list[dict[str, object]]:
    """Run one budgeted retrospective active-learning strategy."""
    selected = set(int(idx) for idx in seed_indices)
    history: list[dict[str, object]] = []
    batch_size = 500
    rounds = 5

    for round_id in range(rounds + 1):
        selected_array = np.array(sorted(selected), dtype=int)
        recovered = len(set(selected_array).intersection(potent_indices))
        history.append(
            {
                "strategy": strategy,
                "round": round_id,
                "labeled_count": int(len(selected_array)),
                "top_potent_recovered": int(recovered),
                "top_potent_recovery_fraction": float(recovered / max(len(potent_indices), 1)),
                "mean_selected_pIC50": float(np.mean(y[selected_array])),
                "selected_scaffold_count": int(len(set(scaffolds[selected_array]))),
            }
        )
        if round_id == rounds:
            break

        pool_indices = np.array([idx for idx in range(len(y)) if idx not in selected], dtype=int)
        model = make_pipeline(StandardScaler(), Ridge(alpha=1.0))
        model.fit(x[selected_array], y[selected_array])
        pred = model.predict(x[pool_indices])
        uncertainty = distance_uncertainty(x, pool_indices, selected_array)

        if strategy == "random":
            chosen = rng.choice(pool_indices, size=min(batch_size, len(pool_indices)), replace=False)
        elif strategy == "highest_predicted_pIC50":
            chosen = pool_indices[np.argsort(pred)[::-1][:batch_size]]
        elif strategy == "uncertainty_sampling":
            chosen = pool_indices[np.argsort(uncertainty)[::-1][:batch_size]]
        elif strategy == "scaffold_diverse_high_score":
            chosen = select_diverse(pool_indices, pred, scaffolds, batch_size)
        elif strategy == "applicability_domain_aware_high_score":
            scores = pred - 0.05 * uncertainty
            chosen = pool_indices[np.argsort(scores)[::-1][:batch_size]]
        elif strategy == "hybrid_activity_uncertainty_diversity":
            scores = pred + 0.05 * uncertainty
            chosen = select_diverse(pool_indices, scores, scaffolds, batch_size)
        else:
            raise ValueError(strategy)

        selected.update(int(idx) for idx in chosen)

    return history


def main() -> None:
    """Run active-learning simulation and save discovery curves."""
    df = load_standardized_or_model_ready().reset_index(drop=True)
    x = df[DESCRIPTOR_COLUMNS].astype(float).to_numpy()
    y = df[TARGET_COLUMN].astype(float).to_numpy()
    scaffolds = df["scaffold_hash"].to_numpy() if "scaffold_hash" in df.columns else np.arange(len(df)).astype(str)

    rng = np.random.default_rng(RANDOM_STATE)
    seed_size = 400
    seed_indices = rng.choice(np.arange(len(y)), size=seed_size, replace=False)
    potent_cutoff = float(np.quantile(y, 0.90))
    potent_indices = set(int(idx) for idx in np.where(y >= potent_cutoff)[0])
    strategies = [
        "random",
        "highest_predicted_pIC50",
        "uncertainty_sampling",
        "scaffold_diverse_high_score",
        "applicability_domain_aware_high_score",
        "hybrid_activity_uncertainty_diversity",
    ]

    rows = []
    for strategy in strategies:
        rows.extend(run_strategy(strategy, x, y, scaffolds, seed_indices, potent_indices, rng))
    history = pd.DataFrame(rows)
    history.to_csv(REPORTS_DIR / "egfr_active_learning_history.csv", index=False)
    final_rows = history.sort_values("round").groupby("strategy").tail(1).copy()
    best = final_rows.sort_values(["top_potent_recovery_fraction", "selected_scaffold_count"], ascending=False).iloc[0]

    metrics = {
        "simulation_rows": int(len(history)),
        "molecule_count": int(len(df)),
        "seed_size": seed_size,
        "batch_size": 500,
        "rounds": 5,
        "model": "Ridge regression on RDKit descriptor features",
        "potent_cutoff_pIC50": potent_cutoff,
        "potent_pool_count": int(len(potent_indices)),
        "strategies": strategies,
        "best_strategy": str(best["strategy"]),
        "best_final_recovery_fraction": float(best["top_potent_recovery_fraction"]),
        "final_strategy_summary": final_rows.to_dict(orient="records"),
    }
    save_json(METRICS_PATH, metrics)

    plt.figure(figsize=(8, 5))
    for strategy in strategies:
        subset = history[history["strategy"] == strategy]
        plt.plot(subset["labeled_count"], subset["top_potent_recovery_fraction"], marker="o", label=strategy)
    plt.xlabel("Labeled molecule budget")
    plt.ylabel("Top-potent recovery fraction")
    plt.title("Retrospective Active-Learning Discovery Curve")
    plt.legend(fontsize=7)
    plt.grid(alpha=0.25)
    save_figure(FIGURES_DIR / "egfr_active_learning_discovery_curve.png")

    plt.figure(figsize=(8, 5))
    for strategy in strategies:
        subset = history[history["strategy"] == strategy]
        plt.plot(subset["labeled_count"], subset["selected_scaffold_count"], marker="o", label=strategy)
    plt.xlabel("Labeled molecule budget")
    plt.ylabel("Selected scaffold count")
    plt.title("Active-Learning Scaffold Diversity")
    plt.legend(fontsize=7)
    plt.grid(alpha=0.25)
    save_figure(FIGURES_DIR / "egfr_active_learning_scaffold_diversity.png")

    report = [
        "# EGFR Retrospective Active-Learning Report",
        "",
        "This simulation uses existing labeled ChEMBL/project records only. No new molecules are generated.",
        "",
        f"- Molecules available: {len(df):,}",
        f"- Initial labeled seed size: {seed_size:,}",
        "- Batch size: 500",
        "- Rounds: 5",
        f"- Model: {metrics['model']}",
        f"- Potent threshold: top 10 percent, pIC50 >= {potent_cutoff:.3f}",
        f"- Best strategy: {metrics['best_strategy']}",
        "",
        "## Final Strategy Summary",
        "",
        markdown_table(final_rows[["strategy", "labeled_count", "top_potent_recovery_fraction", "selected_scaffold_count", "mean_selected_pIC50"]]),
        "",
    ]
    write_text(REPORT_PATH, "\n".join(report))

    print(f"Active-learning strategies: {len(strategies)}")
    print(f"Best strategy: {metrics['best_strategy']}")


if __name__ == "__main__":
    main()
