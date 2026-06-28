"""Cross-validate EGFR QSAR Morgan fingerprint baselines."""

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold


INPUT_PATH = Path("data/processed/egfr_model_ready.csv")
OUTPUT_PATH = Path("results/cross_validation_metrics.csv")

RANDOM_STATE = 42
N_SPLITS = 5
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048
TARGET_COLUMN = "median_pIC50"


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


def scaffold_from_smiles(smiles: str, row_index: int) -> str:
    """Generate a Bemis-Murcko scaffold group label for one SMILES."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return f"invalid_smiles_{row_index}"

    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)

    # Empty scaffolds happen for acyclic molecules. Treat each as its own group
    # so unrelated acyclic molecules are not forced into one large bucket.
    if not scaffold:
        return f"acyclic_{row_index}"

    return scaffold


def calculate_scaffold_groups(smiles: pd.Series) -> np.ndarray:
    """Calculate scaffold group labels for GroupKFold validation."""
    return np.array([scaffold_from_smiles(smiles_value, row_idx) for row_idx, smiles_value in smiles.items()])


def get_models() -> dict[str, object]:
    """Create cross-validation baseline models."""
    return {
        "Dummy mean": DummyRegressor(strategy="mean"),
        "Random forest": RandomForestRegressor(
            n_estimators=100,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        ),
    }


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    """Calculate regression metrics on pIC50 predictions."""
    return {
        "MAE": mean_absolute_error(y_true, y_pred),
        "RMSE": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "R2": r2_score(y_true, y_pred),
    }


def evaluate_cv(
    scheme: str,
    splitter: KFold | GroupKFold,
    x: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray | None = None,
) -> list[dict[str, object]]:
    """Run one cross-validation scheme for all baseline models."""
    rows = []
    split_iterator = splitter.split(x, y, groups) if groups is not None else splitter.split(x, y)

    for fold_number, (train_idx, test_idx) in enumerate(split_iterator, start=1):
        print(f"{scheme} fold {fold_number}: train={len(train_idx)}, test={len(test_idx)}")

        if groups is not None:
            train_groups = set(groups[train_idx])
            test_groups = set(groups[test_idx])
            overlap_count = len(train_groups.intersection(test_groups))
            print(f"{scheme} fold {fold_number}: train_scaffolds={len(train_groups)}, test_scaffolds={len(test_groups)}, overlap={overlap_count}")
        else:
            overlap_count = np.nan

        x_train = x[train_idx]
        x_test = x[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        for model_name, model in get_models().items():
            print(f"Training {scheme} fold {fold_number} | {model_name}")
            model.fit(x_train, y_train)
            y_pred = model.predict(x_test)
            metrics = calculate_metrics(y_test, y_pred)

            rows.append(
                {
                    "validation_scheme": scheme,
                    "fold": fold_number,
                    "model": model_name,
                    "train_size": len(train_idx),
                    "test_size": len(test_idx),
                    "scaffold_overlap": overlap_count,
                    **metrics,
                }
            )

    return rows


def print_fold_metrics(metrics_df: pd.DataFrame) -> None:
    """Print fold-level metrics as a compact table."""
    display = metrics_df.copy()
    for column in ["MAE", "RMSE", "R2"]:
        display[column] = display[column].map(lambda value: f"{value:.3f}")

    print("")
    print("Fold-level metrics")
    print("------------------")
    print(display.to_string(index=False))


def print_summary_metrics(metrics_df: pd.DataFrame) -> None:
    """Print mean and standard deviation for each scheme/model."""
    summary = (
        metrics_df.groupby(["validation_scheme", "model"])[["MAE", "RMSE", "R2"]]
        .agg(["mean", "std"])
        .reset_index()
    )

    print("")
    print("Cross-validation summary")
    print("------------------------")
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.3f}"))


def main() -> None:
    """Run random KFold and scaffold GroupKFold validation."""
    df = pd.read_csv(INPUT_PATH).reset_index(drop=True)
    print(f"Input rows: {len(df)}")

    x, invalid_smiles_count = calculate_morgan_fingerprints(df["canonical_smiles"])
    y = df[TARGET_COLUMN].to_numpy()
    scaffold_groups = calculate_scaffold_groups(df["canonical_smiles"])

    print(f"Fingerprint radius: {FINGERPRINT_RADIUS}")
    print(f"Fingerprint bits: {FINGERPRINT_BITS}")
    print(f"Invalid SMILES: {invalid_smiles_count}")
    print(f"Unique scaffold groups: {len(set(scaffold_groups))}")

    random_kfold = KFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    scaffold_group_kfold = GroupKFold(n_splits=N_SPLITS)

    rows = []
    rows.extend(evaluate_cv("random_kfold", random_kfold, x, y))
    rows.extend(evaluate_cv("scaffold_groupkfold", scaffold_group_kfold, x, y, scaffold_groups))

    metrics_df = pd.DataFrame(rows)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(OUTPUT_PATH, index=False)

    print_fold_metrics(metrics_df)
    print_summary_metrics(metrics_df)
    print("")
    print(f"Saved cross-validation metrics to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
