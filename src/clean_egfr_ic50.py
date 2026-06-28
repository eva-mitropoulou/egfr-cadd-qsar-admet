from pathlib import Path

# pyrefly: ignore [missing-import]
import numpy as np
import pandas as pd

RAW_PATH = Path("data/raw/egfr_chembl_ic50_raw.csv")
OUTPUT_PATH = Path("data/processed/egfr_ic50_clean.csv")

df = pd.read_csv(RAW_PATH)
print(f"Raw rows: {len(df)}")

df = df[df["standard_units"] == "nM"].copy()
print(f"After keeping nM units: {len(df)}")

df = df[df["standard_relation"] == "="].copy()
print(f"After keeping exact '=' relations: {len(df)}")

df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
print(f"After numeric conversion: {len(df)}")

df = df[df["standard_value"].notna()].copy()
print(f"After removing missing IC50 values: {len(df)}")

df = df[df["standard_value"] > 0].copy()
print(f"After keeping positive IC50 values: {len(df)}")

df = df[df["canonical_smiles"].notna()].copy()
print(f"After removing missing SMILES: {len(df)}")

df["pIC50"] = 9 - np.log10(df["standard_value"])
df["IC50_nM"] = df["standard_value"]
print(f"After pIC50 calculation: {len(df)}")

clean = (
    df.groupby(["molecule_chembl_id", "canonical_smiles"], as_index=False)
    .agg(
        median_pIC50=("pIC50", "median"),
        median_IC50_nM=("IC50_nM", "median"),
        n_measurements=("IC50_nM", "count"),
    )
)
print(f"After aggregating duplicate molecules: {len(clean)}")

OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
clean.to_csv(OUTPUT_PATH, index=False)
print(f"Saved clean dataset to {OUTPUT_PATH}")
