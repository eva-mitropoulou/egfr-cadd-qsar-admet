"""Audit standardized molecular representations without exposing SMILES."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, STANDARDIZED_PATH, read_json, save_json, write_text  # noqa: E402


def formal_charge(smiles: str) -> int | None:
    """Return formal charge for a standardized molecule."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))


def main() -> None:
    """Append molecular representation audit details to standardization metrics/report."""
    if not STANDARDIZED_PATH.exists():
        raise FileNotFoundError(f"Missing standardized molecule table: {STANDARDIZED_PATH}")

    df = pd.read_csv(STANDARDIZED_PATH)
    charges = df["standardized_smiles"].map(formal_charge)
    metrics_path = METRICS_DIR / "molecular_standardization_metrics.json"
    metrics = read_json(metrics_path)
    metrics.update(
        {
            "representation_audit_rows": int(len(df)),
            "missing_standardized_smiles": int(df["standardized_smiles"].isna().sum()),
            "duplicate_molecule_hash_count": int(df.duplicated("molecule_hash").sum()),
            "duplicate_scaffold_hash_count": int(df.duplicated("scaffold_hash").sum()) if "scaffold_hash" in df.columns else None,
            "unique_scaffold_count": int(df["scaffold_hash"].nunique()) if "scaffold_hash" in df.columns else None,
            "nonzero_formal_charge_count": int((charges.fillna(0) != 0).sum()),
            "formal_charge_distribution": {str(key): int(value) for key, value in charges.value_counts(dropna=False).items()},
        }
    )
    save_json(metrics_path, metrics)

    lines = [
        "",
        "## Representation Audit",
        "",
        f"- Rows audited: {metrics['representation_audit_rows']:,}",
        f"- Missing standardized representations: {metrics['missing_standardized_smiles']:,}",
        f"- Duplicate molecule hashes: {metrics['duplicate_molecule_hash_count']:,}",
        f"- Unique scaffold hashes: {metrics['unique_scaffold_count']:,}",
        f"- Nonzero formal charge count: {metrics['nonzero_formal_charge_count']:,}",
        "",
    ]
    with (REPORTS_DIR / "molecular_standardization_report.md").open("a", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    print(f"Representation audit rows: {metrics['representation_audit_rows']}")
    print(f"Unique scaffold hashes: {metrics['unique_scaffold_count']}")


if __name__ == "__main__":
    main()
