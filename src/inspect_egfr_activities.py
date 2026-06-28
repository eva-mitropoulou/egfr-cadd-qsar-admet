# pyrefly: ignore [missing-import]
from chembl_webresource_client.new_client import new_client

activity=new_client.activity
records=activity.filter(target_chembl_id='CHEMBL203',standard_type='IC50')

fields = [
      "molecule_chembl_id",
      "canonical_smiles",
      "standard_type",
      "standard_relation",
      "standard_value",
      "standard_units",
      "assay_chembl_id",
      "document_chembl_id",
  ]

print(*fields, sep="\t")

for record in records [:10]:
    row=[record.get(field) for field in fields]
    print(*row,sep="\t")
