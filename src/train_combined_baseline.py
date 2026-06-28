"""Train combined descriptor + fingerprint QSAR baselines for EGFR pIC50."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


INPUT_PATH = Path("data/processed/egfr_model_ready.csv")
METRICS_PATH = Path("results/combined_baseline_metrics.csv")
DESCRIPTOR_METRICS_PATH = Path("results/descriptor_baseline_metrics.csv")
FINGERPRINT_METRICS_PATH = Path("results/fingerprint_baseline_metrics.csv")
PLOT_PATH = Path("figures/combined_predicted_vs_actual.png")

RANDOM_STATE = 42
TEST_SIZE = 0.2
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048
TARGET_COLUMN = "median_pIC50"

DESCRIPTOR_COLUMNS = [
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


def calculate_morgan_fingerprints(smiles: pd.Series) -> tuple[np.ndarray, int]:
    """Convert canonical SMILES into Morgan fingerprint bit vectors."""
    generator = GetMorganGenerator(radius=FINGERPRINT_RADIUS, fpSize=FINGERPRINT_BITS)
    fingerprints = np.zeros((len(smiles), FINGERPRINT_BITS), dtype=np.uint8)
    invalid_smiles_count = 0

    for row_idx, smiles_value in enumerate(smiles):
        mol = Chem.MolFromSmiles(smiles_value)

        if mol is None:
            invalid_smiles_count += 1
            continue

        bit_vector = generator.GetFingerprint(mol)
        fingerprints[row_idx, list(bit_vector.GetOnBits())] = 1

    return fingerprints, invalid_smiles_count


def build_combined_features(df: pd.DataFrame) -> tuple[np.ndarray, int]:
    """Combine Morgan fingerprint bits with RDKit descriptor columns."""
    fingerprints, invalid_smiles_count = calculate_morgan_fingerprints(df["canonical_smiles"])
    descriptors = df[DESCRIPTOR_COLUMNS].to_numpy(dtype=float)

    return np.hstack([fingerprints, descriptors]), invalid_smiles_count


def get_models() -> dict[str, object]:
    """Create simple combined-feature baseline regressors."""
    return {
        "Dummy mean": DummyRegressor(strategy="mean"),
        # Ridge is sensitive to feature scale, so scaling is fitted inside the
        # pipeline after the train/test split to avoid leakage.
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
    """Save a predicted-vs-actual plot for the best combined model."""
    lower = min(float(np.min(y_true)), float(np.min(y_pred)))
    upper = max(float(np.max(y_true)), float(np.max(y_pred)))

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(6, 6))
    plt.scatter(y_true, y_pred, s=18, alpha=0.55)
    plt.plot([lower, upper], [lower, upper], color="black", linestyle="--", linewidth=1)
    plt.xlabel("Actual median pIC50")
    plt.ylabel("Predicted median pIC50")
    plt.title(f"Combined baseline: {model_name}")
    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=200)
    plt.close()


def print_previous_baseline_comparison(combined_metrics: pd.DataFrame) -> None:
    """Print RF comparison against descriptor-only and fingerprint-only baselines."""
    comparison_rows = []

    if DESCRIPTOR_METRICS_PATH.exists():
        descriptor_metrics = pd.read_csv(DESCRIPTOR_METRICS_PATH)
        descriptor_rf = descriptor_metrics[descriptor_metrics["model"] == "Random forest"]
        if not descriptor_rf.empty:
            row = descriptor_rf.iloc[0]
            comparison_rows.append(
                {
                    "feature_set": "descriptors",
                    "model": "Random forest",
                    "MAE": row["MAE"],
                    "RMSE": row["RMSE"],
                    "R2": row["R2"],
                }
            )

    if FINGERPRINT_METRICS_PATH.exists():
        fingerprint_metrics = pd.read_csv(FINGERPRINT_METRICS_PATH)
        fingerprint_rf = fingerprint_metrics[fingerprint_metrics["model"] == "Random forest"]
        if not fingerprint_rf.empty:
            row = fingerprint_rf.iloc[0]
            comparison_rows.append(
                {
                    "feature_set": "fingerprints",
                    "model": "Random forest",
                    "MAE": row["MAE"],
                    "RMSE": row["RMSE"],
                    "R2": row["R2"],
                }
            )

    combined_rf = combined_metrics[combined_metrics["model"] == "Random forest"]
    if not combined_rf.empty:
        row = combined_rf.iloc[0]
        comparison_rows.append(
            {
                "feature_set": "fingerprints+descriptors",
                "model": "Random forest",
                "MAE": row["MAE"],
                "RMSE": row["RMSE"],
                "R2": row["R2"],
            }
        )

    if not comparison_rows:
        return

    comparison_df = pd.DataFrame(comparison_rows)
    print("")
    print("Random Forest Feature-Set Comparison")
    print("------------------------------------")
    print(comparison_df.to_string(index=False, float_format=lambda value: f"{value:.3f}"))


def main() -> None:
    """Train combined-feature baselines on a random train/test split."""
    df = pd.read_csv(INPUT_PATH)
    print(f"Input rows: {len(df)}")

    x, invalid_smiles_count = build_combined_features(df)
    y = df[TARGET_COLUMN].to_numpy()

    print(f"Fingerprint radius: {FINGERPRINT_RADIUS}")
    print(f"Fingerprint bits: {FINGERPRINT_BITS}")
    print(f"Descriptor columns: {len(DESCRIPTOR_COLUMNS)}")
    print(f"Combined feature count: {x.shape[1]}")
    print(f"Invalid SMILES: {invalid_smiles_count}")

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
    print_previous_baseline_comparison(metrics_df)
    print("")
    print(f"Saved metrics to {METRICS_PATH}")
    print(f"Saved predicted-vs-actual plot to {PLOT_PATH}")
    print(f"Best combined baseline: {best_model_name}")


if __name__ == "__main__":
    main()
