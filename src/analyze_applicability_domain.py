"""Analyze QSAR applicability domain with scaffold CV and Tanimoto similarity."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold


INPUT_PATH = Path("data/processed/egfr_model_ready.csv")
PREDICTIONS_PATH = Path("results/applicability_domain_predictions.csv")
SUMMARY_PATH = Path("results/applicability_domain_summary.csv")
PLOT_PATH = Path("figures/error_vs_max_tanimoto.png")

RANDOM_STATE = 42
N_SPLITS = 5
FINGERPRINT_RADIUS = 2
FINGERPRINT_BITS = 2048
TARGET_COLUMN = "median_pIC50"


def scaffold_from_smiles(smiles: str, row_index: int) -> str:
    """Generate a Bemis-Murcko scaffold group label for one SMILES."""
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return f"invalid_smiles_{row_index}"

    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)

    # Empty scaffolds occur for acyclic molecules. Treat each as its own group
    # so unrelated acyclic molecules are not forced into one artificial scaffold.
    if not scaffold:
        return f"acyclic_{row_index}"

    return scaffold


def calculate_scaffold_groups(smiles: pd.Series) -> np.ndarray:
    """Calculate scaffold labels for GroupKFold validation."""
    return np.array([scaffold_from_smiles(smiles_value, row_idx) for row_idx, smiles_value in smiles.items()])


def calculate_fingerprints(smiles: pd.Series) -> tuple[list[DataStructs.ExplicitBitVect | None], np.ndarray, int]:
    """Calculate RDKit bit vectors and dense arrays from canonical SMILES."""
    generator = GetMorganGenerator(radius=FINGERPRINT_RADIUS, fpSize=FINGERPRINT_BITS)
    bit_vectors: list[DataStructs.ExplicitBitVect | None] = []
    dense = np.zeros((len(smiles), FINGERPRINT_BITS), dtype=np.uint8)
    invalid_smiles_count = 0

    for row_idx, smiles_value in enumerate(smiles):
        mol = Chem.MolFromSmiles(smiles_value)

        if mol is None:
            bit_vectors.append(None)
            invalid_smiles_count += 1
            continue

        bit_vector = generator.GetFingerprint(mol)
        bit_vectors.append(bit_vector)
        dense[row_idx, list(bit_vector.GetOnBits())] = 1

    return bit_vectors, dense, invalid_smiles_count


def max_tanimoto_to_training(
    test_fingerprint: DataStructs.ExplicitBitVect | None,
    train_fingerprints: list[DataStructs.ExplicitBitVect],
) -> float:
    """Return maximum Tanimoto similarity from one test molecule to train molecules."""
    if test_fingerprint is None or not train_fingerprints:
        return np.nan

    similarities = DataStructs.BulkTanimotoSimilarity(test_fingerprint, train_fingerprints)
    return float(max(similarities))


def calculate_metrics(values: pd.Series) -> pd.Series:
    """Calculate count, MAE, and RMSE for one similarity bin."""
    true_values = values["true_pIC50"].to_numpy()
    predicted_values = values["predicted_pIC50"].to_numpy()

    return pd.Series(
        {
            "count": len(values),
            "MAE": mean_absolute_error(true_values, predicted_values),
            "RMSE": float(np.sqrt(mean_squared_error(true_values, predicted_values))),
        }
    )


def summarize_by_similarity_bin(predictions: pd.DataFrame) -> pd.DataFrame:
    """Summarize prediction error by max-Tanimoto similarity bin."""
    bins = [-float("inf"), 0.3, 0.5, 0.7, float("inf")]
    labels = ["<0.3", "0.3-0.5", "0.5-0.7", ">0.7"]

    predictions = predictions.copy()
    predictions["similarity_bin"] = pd.cut(
        predictions["max_tanimoto_to_train"],
        bins=bins,
        labels=labels,
        right=False,
    )

    summary = (
        predictions.groupby("similarity_bin", observed=False)
        .apply(calculate_metrics, include_groups=False)
        .reset_index()
    )

    return summary


def save_error_similarity_plot(predictions: pd.DataFrame) -> None:
    """Save scatter plot of absolute error versus nearest-training similarity."""
    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7, 5))
    plt.scatter(
        predictions["max_tanimoto_to_train"],
        predictions["absolute_error"],
        s=12,
        alpha=0.3,
    )
    plt.xlabel("Max Tanimoto similarity to training set")
    plt.ylabel("Absolute pIC50 error")
    plt.title("Applicability Domain: Error vs Training Similarity")
    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=200)
    plt.close()


def main() -> None:
    """Run scaffold-CV applicability-domain analysis."""
    df = pd.read_csv(INPUT_PATH).reset_index(drop=True)
    print(f"Input rows: {len(df)}")

    scaffold_groups = calculate_scaffold_groups(df["canonical_smiles"])
    bit_vectors, dense_fingerprints, invalid_smiles_count = calculate_fingerprints(df["canonical_smiles"])
    y = df[TARGET_COLUMN].to_numpy()

    print(f"Fingerprint radius: {FINGERPRINT_RADIUS}")
    print(f"Fingerprint bits: {FINGERPRINT_BITS}")
    print(f"Invalid SMILES: {invalid_smiles_count}")
    print(f"Unique scaffold groups: {len(set(scaffold_groups))}")

    group_kfold = GroupKFold(n_splits=N_SPLITS)
    rows = []

    for fold_number, (train_idx, test_idx) in enumerate(
        group_kfold.split(dense_fingerprints, y, scaffold_groups),
        start=1,
    ):
        train_scaffolds = set(scaffold_groups[train_idx])
        test_scaffolds = set(scaffold_groups[test_idx])
        scaffold_overlap = train_scaffolds.intersection(test_scaffolds)

        print(
            f"Fold {fold_number}: train={len(train_idx)}, test={len(test_idx)}, "
            f"train_scaffolds={len(train_scaffolds)}, test_scaffolds={len(test_scaffolds)}, "
            f"overlap={len(scaffold_overlap)}"
        )

        model = RandomForestRegressor(
            n_estimators=100,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        model.fit(dense_fingerprints[train_idx], y[train_idx])
        predictions = model.predict(dense_fingerprints[test_idx])

        train_bit_vectors = [bit_vectors[idx] for idx in train_idx if bit_vectors[idx] is not None]

        for local_idx, row_index in enumerate(test_idx):
            true_value = y[row_index]
            predicted_value = predictions[local_idx]
            max_similarity = max_tanimoto_to_training(bit_vectors[row_index], train_bit_vectors)

            rows.append(
                {
                    "molecule_chembl_id": df.loc[row_index, "molecule_chembl_id"],
                    "canonical_smiles": df.loc[row_index, "canonical_smiles"],
                    "fold": fold_number,
                    "scaffold": scaffold_groups[row_index],
                    "true_pIC50": true_value,
                    "predicted_pIC50": predicted_value,
                    "absolute_error": abs(true_value - predicted_value),
                    "max_tanimoto_to_train": max_similarity,
                }
            )

    predictions_df = pd.DataFrame(rows)
    summary_df = summarize_by_similarity_bin(predictions_df)

    PREDICTIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    predictions_df.to_csv(PREDICTIONS_PATH, index=False)
    summary_df.to_csv(SUMMARY_PATH, index=False)
    save_error_similarity_plot(predictions_df)

    print("")
    print("Applicability-domain summary")
    print("----------------------------")
    print(summary_df.to_string(index=False, float_format=lambda value: f"{value:.3f}"))
    print("")
    print(f"Saved out-of-fold predictions to {PREDICTIONS_PATH}")
    print(f"Saved applicability-domain summary to {SUMMARY_PATH}")
    print(f"Saved error-vs-similarity plot to {PLOT_PATH}")


if __name__ == "__main__":
    main()
