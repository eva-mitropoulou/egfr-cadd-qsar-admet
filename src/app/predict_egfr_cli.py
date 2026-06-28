"""CLI for EGFR QSAR prediction with applicability-domain warnings."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    MODELS_DIR,
    REPORTS_DIR,
    calculate_morgan_matrix,
    load_standardized_or_model_ready,
    max_tanimoto_to_train,
    stable_hash,
    write_text,
)
from triage.admet_risk_scoring import lipinski_violations, model_risk_from_similarity, uncertainty_penalty  # noqa: E402


MODEL_PATH = MODELS_DIR / "egfr_primary_full_model.joblib"


def standardize_smiles(smiles: str) -> str | None:
    """Return canonical isomeric SMILES for CLI input."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def descriptor_row(smiles: str) -> dict[str, float] | None:
    """Calculate basic drug-likeness descriptors for CLI output."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return {
        "QED": float(QED.qed(mol)),
        "MolWt": float(Descriptors.MolWt(mol)),
        "MolLogP": float(Crippen.MolLogP(mol)),
        "TPSA": float(Descriptors.TPSA(mol)),
        "NumHDonors": int(Lipinski.NumHDonors(mol)),
        "NumHAcceptors": int(Lipinski.NumHAcceptors(mol)),
        "NumRotatableBonds": int(Lipinski.NumRotatableBonds(mol)),
    }


def create_example(path: Path, n: int = 5) -> None:
    """Create a small example input from existing project molecules."""
    df = load_standardized_or_model_ready()
    smiles_column = "standardized_smiles" if "standardized_smiles" in df.columns else "canonical_smiles"
    example = df[["molecule_chembl_id", smiles_column]].iloc[:n].copy()
    example = example.rename(columns={"molecule_chembl_id": "molecule_id", smiles_column: "smiles"})
    path.parent.mkdir(parents=True, exist_ok=True)
    example.to_csv(path, index=False)


def tree_prediction_stats(model, matrix) -> tuple[np.ndarray, np.ndarray]:
    """Return mean/std predictions across RF trees if available."""
    if hasattr(model, "estimators_"):
        predictions = np.vstack([tree.predict(matrix) for tree in model.estimators_])
        return predictions.mean(axis=0), predictions.std(axis=0)
    pred = model.predict(matrix)
    return pred, np.zeros(len(pred))


def main() -> None:
    """Run CLI prediction."""
    parser = argparse.ArgumentParser(description="Predict EGFR pIC50 with model-risk warnings.")
    parser.add_argument("--input", type=Path, default=None, help="Input CSV with molecule_id and smiles columns.")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV path.")
    parser.add_argument("--create-example", type=Path, default=None, help="Create a small example input CSV before predicting.")
    args = parser.parse_args()

    if args.create_example is not None:
        create_example(args.create_example)
        if args.input is None:
            args.input = args.create_example

    if args.input is None:
        raise ValueError("--input is required unless --create-example is provided.")
    if not MODEL_PATH.exists():
        raise FileNotFoundError(f"Missing primary model: {MODEL_PATH}")

    model = joblib.load(MODEL_PATH)
    input_df = pd.read_csv(args.input)
    required = {"molecule_id", "smiles"}
    if not required.issubset(input_df.columns):
        raise ValueError("Input must contain molecule_id and smiles columns.")

    standardized = [standardize_smiles(value) for value in input_df["smiles"]]
    valid_mask = pd.Series(standardized).notna()
    valid_smiles = [value for value in standardized if value is not None]
    matrix, input_fps, invalid_count = calculate_morgan_matrix(valid_smiles)
    predicted, std = tree_prediction_stats(model, matrix)

    training = load_standardized_or_model_ready()
    train_smiles_column = "standardized_smiles" if "standardized_smiles" in training.columns else "canonical_smiles"
    _train_matrix, train_fps, _train_invalid = calculate_morgan_matrix(training[train_smiles_column])
    train_fps_valid = [fp for fp in train_fps if fp is not None]

    output_rows: list[dict[str, object]] = []
    valid_counter = 0
    for row_idx, input_row in input_df.iterrows():
        if not valid_mask.iloc[row_idx]:
            output_rows.append(
                {
                    "molecule_id": input_row["molecule_id"],
                    "molecule_hash": stable_hash(input_row["molecule_id"]),
                    "valid_molecule": False,
                    "predicted_pIC50": np.nan,
                    "rf_prediction_std": np.nan,
                    "max_tanimoto_to_train": np.nan,
                    "model_risk_category": "invalid",
                    "applicability_warning": "invalid_molecule",
                }
            )
            continue

        descriptors = descriptor_row(valid_smiles[valid_counter]) or {}
        similarity = max_tanimoto_to_train(input_fps[valid_counter], train_fps_valid)
        temp = pd.DataFrame([descriptors])
        lipinski_count = int(lipinski_violations(temp).iloc[0]) if not temp.empty else None
        risk = model_risk_from_similarity(similarity)
        output_rows.append(
            {
                "molecule_id": input_row["molecule_id"],
                "molecule_hash": stable_hash(valid_smiles[valid_counter]),
                "valid_molecule": True,
                "predicted_pIC50": float(predicted[valid_counter]),
                "rf_prediction_std": float(std[valid_counter]),
                "max_tanimoto_to_train": float(similarity),
                "model_risk_category": risk,
                "applicability_warning": "out_of_domain" if risk == "high" else "within_domain",
                "uncertainty_penalty": uncertainty_penalty(float(std[valid_counter])),
                "QED": descriptors.get("QED"),
                "MolWt": descriptors.get("MolWt"),
                "MolLogP": descriptors.get("MolLogP"),
                "TPSA": descriptors.get("TPSA"),
                "lipinski_violations": lipinski_count,
            }
        )
        valid_counter += 1

    output = pd.DataFrame(output_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)

    report = [
        "# EGFR CLI Demo Report",
        "",
        f"- Input file: `{args.input}`",
        f"- Output file: `{args.output}`",
        f"- Input rows: {len(input_df)}",
        f"- Valid molecules: {int(output['valid_molecule'].sum())}",
        f"- Invalid molecules: {invalid_count}",
        "",
        "The output excludes raw SMILES and includes applicability-domain warnings and uncertainty scores.",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_cli_demo_report.md", "\n".join(report))

    print(f"CLI predictions: {len(output)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
