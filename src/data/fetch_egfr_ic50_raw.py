
from pathlib import Path
import pandas as pd

# pyrefly: ignore[missing-import]
from chembl_webresource_client.new_client import new_client

TARGET_CHEMBL_ID="CHEMBL203"
OUTPUT_PATH=Path('data/raw/egfr_chembl_ic50_raw.csv')

fields = [
      "molecule_chembl_id",
      "canonical_smiles",
      "standard_type",
      "standard_relation",
      "standard_value",
      "standard_units",
      "assay_chembl_id",
      "document_chembl_id",
      "target_chembl_id",
  ]

activity = new_client.activity
records=activity.filter(
    target_chembl_id=TARGET_CHEMBL_ID,
    standard_type='IC50'
)

rows=[]
for record in records:
    row ={field:record.get(field) for field in fields}
    rows.append(row) 

df=pd.DataFrame(rows)
    
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUTPUT_PATH, index=False)

print(f"Saved {len(df)} rows to {OUTPUT_PATH}")

print("\nstandard_units counts:")
print(df["standard_units"].value_counts(dropna=False))

print("\nstandard_relation counts:")
print(df["standard_relation"].value_counts(dropna=False))
