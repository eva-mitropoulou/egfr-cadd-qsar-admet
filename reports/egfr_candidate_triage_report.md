# EGFR Candidate Triage Report

This is ADMET-style, drug-likeness, uncertainty, and model-risk triage over existing ChEMBL/project molecules.
It is drug-likeness and model-risk triage over existing molecules.

- Ranked molecules: 10,593
- Diverse top-20 unique scaffolds: 20
- Diverse top-20 low/medium model risk: 20/20
- Diverse top-20 Lipinski-clean: 18/20
- Top-20 with no medicinal-chemistry alert: 20/20
- Diverse top-20 with no medicinal-chemistry alert: 20/20
- PAINS-flagged ranked molecules: 847
- Brenk-flagged ranked molecules: 6,074
- External unwanted-substructure flagged ranked molecules: 0
- Medicinal-chemistry alert status: medchem_alerts_available

PAINS, Brenk, and external unwanted-substructure SMARTS alerts are risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.

## Diverse Top 20 Summary

| molecule_chembl_id | predicted_pIC50 | model_risk_category | triage_risk_bin | QED | lipinski_violations | medchem_alert_flag | medchem_alert_count | final_score | final_triage_category |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CHEMBL5997498 | 10.029 | low | low | 0.595 | 0 | False | 0 | 10.625 | prioritize |
| CHEMBL5790648 | 9.872 | low | low | 0.656 | 0 | False | 0 | 10.529 | prioritize |
| CHEMBL174426 | 9.710 | low | low | 0.591 | 0 | False | 0 | 10.300 | prioritize |
| CHEMBL4749862 | 9.338 | low | low | 0.488 | 0 | False | 0 | 9.826 | prioritize |
| CHEMBL2031299 | 9.131 | low | low | 0.588 | 0 | False | 0 | 9.719 | prioritize |
| CHEMBL4789312 | 9.030 | low | low | 0.523 | 0 | False | 0 | 9.553 | prioritize |
| CHEMBL2031298 | 9.000 | low | low | 0.552 | 0 | False | 0 | 9.551 | prioritize |
| CHEMBL5269375 | 9.113 | low | low | 0.411 | 0 | False | 0 | 9.524 | prioritize |
| CHEMBL53753 | 8.689 | low | low | 0.769 | 0 | False | 0 | 9.458 | prioritize |
| CHEMBL5956184 | 9.042 | low | low | 0.411 | 0 | False | 0 | 9.454 | prioritize |
| CHEMBL2031300 | 8.838 | low | low | 0.574 | 0 | False | 0 | 9.411 | prioritize |
| CHEMBL6026192 | 9.265 | low | medium | 0.388 | 1 | False | 0 | 9.404 | prioritize |
| CHEMBL4164805 | 9.000 | low | low | 0.386 | 0 | False | 0 | 9.386 | prioritize |
| CHEMBL54400 | 8.598 | low | low | 0.769 | 0 | False | 0 | 9.367 | prioritize |
| CHEMBL6032856 | 9.163 | low | medium | 0.406 | 1 | False | 0 | 9.320 | prioritize |
| CHEMBL5087955 | 9.090 | medium | low | 0.765 | 0 | False | 0 | 9.304 | prioritize |
| CHEMBL128987 | 8.634 | low | low | 0.623 | 0 | False | 0 | 9.257 | prioritize |
| CHEMBL162622 | 8.600 | low | low | 0.623 | 0 | False | 0 | 9.223 | prioritize |
| CHEMBL294475 | 8.541 | low | low | 0.617 | 0 | False | 0 | 9.158 | prioritize |
| CHEMBL405398 | 8.503 | low | low | 0.623 | 0 | False | 0 | 9.126 | prioritize |

Detailed ranked output is saved to `reports/egfr_ranked_existing_molecules.csv`.
