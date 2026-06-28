"""Build RDKit descriptor feature table for standardized EGFR molecules."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    METRICS_DIR,
    PROCESSED_DIR,
    REPORTS_DIR,
    TARGET_COLUMN,
    load_standardized_or_model_ready,
    save_json,
    write_text,
)


OUTPUT_PATH = PROCESSED_DIR / "features_rdkit_descriptors.csv"


DESCRIPTOR_FUNCTIONS = {
    "MolWt": Descriptors.MolWt,
    "MolLogP": Crippen.MolLogP,
    "TPSA": Descriptors.TPSA,
    "NumHDonors": Lipinski.NumHDonors,
    "NumHAcceptors": Lipinski.NumHAcceptors,
    "NumRotatableBonds": Lipinski.NumRotatableBonds,
    "RingCount": Descriptors.RingCount,
    "HeavyAtomCount": Descriptors.HeavyAtomCount,
    "QED": QED.qed,
}


def descriptors_for_smiles(smiles: str) -> dict[str, float] | None:
    """Calculate descriptor values for one molecule."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return {name: float(func(mol)) for name, func in DESCRIPTOR_FUNCTIONS.items()}


def main() -> None:
    """Create descriptor feature matrix and feature-generation metrics."""
    df = load_standardized_or_model_ready()
    smiles_column = "standardized_smiles" if "standardized_smiles" in df.columns else "canonical_smiles"

    rows = []
    invalid_count = 0
    for _, row in df.iterrows():
        values = descriptors_for_smiles(row[smiles_column])
        if values is None:
            invalid_count += 1
            continue
        rows.append(
            {
                "molecule_chembl_id": row.get("molecule_chembl_id"),
                "molecule_hash": row.get("molecule_hash"),
                "scaffold_hash": row.get("scaffold_hash"),
                TARGET_COLUMN: row[TARGET_COLUMN],
                **values,
            }
        )

    features = pd.DataFrame(rows)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(OUTPUT_PATH, index=False)

    descriptor_columns = list(DESCRIPTOR_FUNCTIONS)
    missing_counts = features[descriptor_columns].isna().sum().to_dict()
    constant_columns = [
        column for column in descriptor_columns if features[column].nunique(dropna=True) <= 1
    ]

    metrics = {
        "input_rows": int(len(df)),
        "descriptor_rows": int(len(features)),
        "invalid_molecule_count": int(invalid_count),
        "descriptor_column_count": len(descriptor_columns),
        "descriptor_columns": descriptor_columns,
        "missing_descriptor_counts": {key: int(value) for key, value in missing_counts.items()},
        "constant_descriptor_columns": constant_columns,
        "feature_label_alignment": int(len(features)) == int(features[TARGET_COLUMN].notna().sum()),
    }

    metrics_path = METRICS_DIR / "feature_generation_metrics.json"
    existing = {}
    if metrics_path.exists():
        import json

        existing = json.loads(metrics_path.read_text(encoding="utf-8"))
    existing["rdkit_descriptors"] = metrics
    save_json(metrics_path, existing)

    lines = [
        "# Feature Generation Report",
        "",
        "## RDKit 2D Descriptors",
        "",
        f"- Input rows: {metrics['input_rows']:,}",
        f"- Descriptor matrix shape: {metrics['descriptor_rows']:,} x {metrics['descriptor_column_count']}",
        f"- Invalid molecule count: {metrics['invalid_molecule_count']:,}",
        f"- Constant descriptor columns: {constant_columns}",
        f"- Feature-label alignment: {metrics['feature_label_alignment']}",
        "",
    ]
    write_text(REPORTS_DIR / "feature_generation_report.md", "\n".join(lines))

    print(f"Descriptor matrix rows: {metrics['descriptor_rows']}")
    print(f"Descriptor columns: {metrics['descriptor_column_count']}")


if __name__ == "__main__":
    main()
