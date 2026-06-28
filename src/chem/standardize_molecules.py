"""Standardize EGFR molecule representations with documented RDKit policies."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    METRICS_DIR,
    MODEL_READY_PATH,
    REPORTS_DIR,
    STANDARDIZED_PATH,
    add_scaffold_hashes,
    save_json,
    stable_hash,
    write_text,
)


def load_molstandardize_tools() -> dict[str, object]:
    """Load RDKit MolStandardize tools when available."""
    tools: dict[str, object] = {"available": False}
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize

        tools.update(
            {
                "available": True,
                "fragment_chooser": rdMolStandardize.LargestFragmentChooser(),
                "normalizer": rdMolStandardize.Normalizer(),
                "uncharger": rdMolStandardize.Uncharger(),
                "tautomer_enumerator": rdMolStandardize.TautomerEnumerator(),
            }
        )
    except Exception:
        tools["available"] = False
    return tools


def standardize_one(smiles: str, tools: dict[str, object]) -> tuple[str | None, dict[str, bool]]:
    """Parse, sanitize, standardize, and canonicalize one molecule."""
    status = {
        "parsed": False,
        "fragment_standardized": False,
        "normalized": False,
        "uncharged": False,
        "tautomer_canonicalized": False,
    }
    mol = Chem.MolFromSmiles(str(smiles), sanitize=True)
    if mol is None:
        return None, status
    status["parsed"] = True

    if tools.get("available"):
        try:
            mol = tools["fragment_chooser"].choose(mol)
            status["fragment_standardized"] = True
            mol = tools["normalizer"].normalize(mol)
            status["normalized"] = True
            mol = tools["uncharger"].uncharge(mol)
            status["uncharged"] = True
            # Full tautomer canonicalization across all ChEMBL EGFR records can be
            # slow and fragile for a portfolio-scale automated finish run. The
            # policy is documented as degraded/future rather than applied silently.
            status["tautomer_canonicalized"] = False
        except Exception:
            # Fallback keeps the sanitized molecule and documents degradation via status flags.
            pass

    standardized = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    return standardized, status


def main() -> None:
    """Create standardized EGFR molecule table and policy report."""
    if not MODEL_READY_PATH.exists():
        raise FileNotFoundError(f"Missing model-ready dataset: {MODEL_READY_PATH}")

    df = pd.read_csv(MODEL_READY_PATH)
    tools = load_molstandardize_tools()

    standardized_smiles: list[str | None] = []
    status_rows: list[dict[str, bool]] = []
    for smiles in df["canonical_smiles"]:
        standardized, status = standardize_one(str(smiles), tools)
        standardized_smiles.append(standardized)
        status_rows.append(status)

    status_df = pd.DataFrame(status_rows)
    result = df.copy()
    result["standardized_smiles"] = standardized_smiles
    result["molecule_hash"] = result["standardized_smiles"].map(stable_hash)
    result["standardization_valid"] = result["standardized_smiles"].notna()
    result = add_scaffold_hashes(result[result["standardization_valid"]].reset_index(drop=True), "standardized_smiles")

    duplicate_standardized = int(result.duplicated("standardized_smiles").sum())
    invalid_count = int(len(df) - len(result))

    STANDARDIZED_PATH.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(STANDARDIZED_PATH, index=False)

    metrics = {
        "input_rows": int(len(df)),
        "standardized_rows": int(len(result)),
        "invalid_molecule_count": invalid_count,
        "duplicate_standardized_molecule_count": duplicate_standardized,
        "unique_standardized_molecule_count": int(result["standardized_smiles"].nunique()),
        "molstandardize_available": bool(tools.get("available")),
        "parsed_count": int(status_df["parsed"].sum()),
        "fragment_standardized_count": int(status_df["fragment_standardized"].sum()),
        "normalized_count": int(status_df["normalized"].sum()),
        "uncharged_count": int(status_df["uncharged"].sum()),
        "tautomer_canonicalized_count": int(status_df["tautomer_canonicalized"].sum()),
        "fragment_policy": "largest fragment chosen when RDKit MolStandardize is available; otherwise sanitized canonical molecule retained",
        "salt_policy": "salts/fragments handled by largest-fragment choice when available",
        "charge_policy": "RDKit Uncharger applied when available; otherwise original sanitized charge state retained",
        "tautomer_policy": "Full tautomer canonicalization is documented but skipped in the automated finish run to avoid fragile long runtimes; sanitized canonical isomeric SMILES are retained",
        "stereochemistry_policy": "isomeric canonical SMILES retained; stereochemistry is not collapsed",
        "public_reporting_policy": "reports use molecule IDs, hashes, and aggregate counts rather than raw SMILES",
    }
    save_json(METRICS_DIR / "molecular_standardization_metrics.json", metrics)

    lines = [
        "# Molecular Standardization Report",
        "",
        f"- Input molecules: {metrics['input_rows']:,}",
        f"- Standardized molecules: {metrics['standardized_rows']:,}",
        f"- Invalid molecules: {metrics['invalid_molecule_count']:,}",
        f"- Duplicate standardized molecules: {metrics['duplicate_standardized_molecule_count']:,}",
        f"- Unique standardized molecules: {metrics['unique_standardized_molecule_count']:,}",
        f"- RDKit MolStandardize available: {metrics['molstandardize_available']}",
        "",
        "## Policies",
        "",
        f"- Fragment/salt handling: {metrics['fragment_policy']}",
        f"- Charge handling: {metrics['charge_policy']}",
        f"- Tautomer handling: {metrics['tautomer_policy']}",
        f"- Stereochemistry handling: {metrics['stereochemistry_policy']}",
        f"- Public reporting: {metrics['public_reporting_policy']}",
        "",
        "Detailed standardized molecule representations are saved to `data/processed/egfr_standardized_molecules.csv`.",
        "",
    ]
    write_text(REPORTS_DIR / "molecular_standardization_report.md", "\n".join(lines))

    print(f"Standardized molecules: {metrics['standardized_rows']}")
    print(f"Invalid molecules: {metrics['invalid_molecule_count']}")
    print(f"Metrics: {METRICS_DIR / 'molecular_standardization_metrics.json'}")


if __name__ == "__main__":
    main()
