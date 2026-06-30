"""Calculate basic RDKit descriptors for cleaned EGFR compounds."""

from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED


INPUT_PATH = Path("data/processed/egfr_ic50_clean.csv")
OUTPUT_PATH = Path("data/processed/egfr_descriptors.csv")


def calculate_descriptors(smiles: str) -> dict[str, float] | None:
    """Calculate basic molecular descriptors for one SMILES string.

    Returns None when RDKit cannot parse the SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    return {
        "MolWt": Descriptors.MolWt(mol),
        "MolLogP": Crippen.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Lipinski.NumHDonors(mol),
        "NumHAcceptors": Lipinski.NumHAcceptors(mol),
        "NumRotatableBonds": Lipinski.NumRotatableBonds(mol),
        "RingCount": Descriptors.RingCount(mol),
        "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
        "QED": QED.qed(mol),
    }


def main() -> None:
    """Read cleaned EGFR activity data and save RDKit descriptor features."""
    df = pd.read_csv(INPUT_PATH)
    print(f"Input rows: {len(df)}")

    rows = []
    invalid_smiles_count = 0

    for _, row in df.iterrows():
        descriptors = calculate_descriptors(row["canonical_smiles"])

        if descriptors is None:
            invalid_smiles_count += 1
            continue

        output_row = {
            "molecule_chembl_id": row["molecule_chembl_id"],
            "canonical_smiles": row["canonical_smiles"],
            "median_pIC50": row["median_pIC50"],
            "median_IC50_nM": row["median_IC50_nM"],
            "n_measurements": row["n_measurements"],
        }
        output_row.update(descriptors)
        rows.append(output_row)

    descriptors_df = pd.DataFrame(rows)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    descriptors_df.to_csv(OUTPUT_PATH, index=False)

    print(f"Valid molecules: {len(descriptors_df)}")
    print(f"Invalid SMILES: {invalid_smiles_count}")
    print(f"Saved descriptors to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
