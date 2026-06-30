# pyrefly: ignore [missing-import]
from chembl_webresource_client.new_client import new_client

target=new_client.target
results=target.search("EGFR")

for record in results:
    print(
        record.get("target_chembl_id"),
        record.get("pref_name"),
        record.get("organism"),
        record.get("target_type"),
        record.get('score'),
        sep="\t",
        
    )

    if (
        record.get('organism')=='Homo sapiens'
        and record.get('target_type')=='SINGLE PROTEIN'
        and record.get('pref_name')=='Epidermal growth factor receptor'
    ):
        print("Likely EGFR target:", record.get("target_chembl_id"), record.get('pref_name'))
    
