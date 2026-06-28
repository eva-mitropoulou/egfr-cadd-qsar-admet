"""Train first descriptor-only baseline QSAR models for EGFR pIC50."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


INPUT_PATH = Path("data/processed/egfr_model_ready.csv")
METRICS_PATH = Path("results/descriptor_baseline_metrics.csv")
PLOT_PATH = Path("figures/descriptor_predicted_vs_actual.png")

RANDOM_STATE = 42
TEST_SIZE = 0.2

FEATURE_COLUMNS = [
    "MolWt",
    "MolLogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "RingCount",
    "HeavyAtomCount",
    "QED",
]
TARGET_COLUMN = "median_pIC50"


def get_models() -> dict[str, object]:
    """Create simple descriptor-only baseline regressors."""
    return {
        "Dummy mean": DummyRegressor(strategy="mean"),
        # Scaling is fitted only on the training set by the pipeline.
        "Ridge regression": Pipeline(
            [
                ("scaler", StandardScaler()),
                ("model", Ridge(alpha=1.0)),
            ]
        ),
        "Random forest": RandomForestRegressor(
            n_estimators=100,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        ),
        "Gradient boosting": GradientBoostingRegressor(random_state=RANDOM_STATE),
    }


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Calculate regression metrics on pIC50 predictions."""
    return {
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "R2": r2_score(y_true, y_pred),
    }


def format_metrics_table(metrics_df: pd.DataFrame) -> str:
    """Format model metrics as a compact Markdown table."""
    display = metrics_df.copy()
    for column in ["MAE", "RMSE", "R2"]:
        display[column] = display[column].map(lambda value: f"{value:.3f}")

    columns = ["model", "train_size", "test_size", "MAE", "RMSE", "R2"]
    header = "| " + " | ".join(columns) + " |"
    separator = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = ["| " + " | ".join(str(row[column]) for column in columns) + " |" for _, row in display.iterrows()]

    return "\n".join([header, separator, *rows])


def save_predicted_vs_actual_plot(y_true: np.ndarray, y_pred: np.ndarray, model_name: str) -> None:
    """Save a predicted-vs-actual plot for the best descriptor model."""
    lower = min(float(np.min(y_true)), float(np.min(y_pred)))
    upper = max(float(np.max(y_true)), float(np.max(y_pred)))

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 6))
    plt.scatter(y_true, y_pred, s=18, alpha=0.55)
    plt.plot([lower, upper], [lower, upper], color="black", linestyle="--", linewidth=1)
    plt.xlabel("Actual median pIC50")
    plt.ylabel("Predicted median pIC50")
    plt.title(f"Descriptor baseline: {model_name}")
    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=200)
    plt.close()


def main() -> None:
    """Train descriptor-only baselines on a random train/test split."""
    df = pd.read_csv(INPUT_PATH)
    print(f"Input rows: {len(df)}")

    # X contains only molecular descriptors. The target and metadata columns are
    # excluded to avoid leakage and non-chemical identifiers.
    x = df[FEATURE_COLUMNS]
    y = df[TARGET_COLUMN].to_numpy()

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
    )
    print(f"Train rows: {len(x_train)}")
    print(f"Test rows: {len(x_test)}")

    results = []
    predictions = {}

    for model_name, model in get_models().items():
        print(f"Training {model_name}")
        model.fit(x_train, y_train)
        y_pred = model.predict(x_test)

        metrics = calculate_metrics(y_test, y_pred)
        results.append(
            {
                "model": model_name,
                "train_size": len(x_train),
                "test_size": len(x_test),
                **metrics,
            }
        )
        predictions[model_name] = y_pred

    metrics_df = pd.DataFrame(results).sort_values("RMSE").reset_index(drop=True)

    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(METRICS_PATH, index=False)

    best_model_name = metrics_df.iloc[0]["model"]
    save_predicted_vs_actual_plot(y_test, predictions[best_model_name], best_model_name)

    print("")
    print(format_metrics_table(metrics_df))
    print("")
    print(f"Saved metrics to {METRICS_PATH}")
    print(f"Saved predicted-vs-actual plot to {PLOT_PATH}")
    print(f"Best descriptor baseline: {best_model_name}")


if __name__ == "__main__":
    main()
