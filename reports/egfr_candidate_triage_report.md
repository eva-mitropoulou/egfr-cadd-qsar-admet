# EGFR Candidate Triage Report

This is ADMET-style, drug-likeness, uncertainty, and model-risk triage over existing ChEMBL/project molecules.
It is drug-likeness and model-risk triage over existing molecules.

- Ranked molecules: 10,593
- Diverse top-20 unique scaffolds: 20
- Diverse top-20 low/medium model risk: 20/20
- Diverse top-20 Lipinski-clean: 18/20
- PAINS/Brenk catalog status: available

## Diverse Top 20 Summary

| molecule_chembl_id | predicted_pIC50 | model_risk_category | QED | lipinski_violations | final_score | final_triage_category |
| --- | --- | --- | --- | --- | --- | --- |
| CHEMBL5997498 | 10.029 | low | 0.595 | 0 | 10.625 | prioritize |
| CHEMBL5790648 | 9.872 | low | 0.656 | 0 | 10.529 | prioritize |
| CHEMBL174426 | 9.710 | low | 0.591 | 0 | 10.300 | prioritize |
| CHEMBL4749862 | 9.338 | low | 0.488 | 0 | 9.826 | prioritize |
| CHEMBL2031299 | 9.131 | low | 0.588 | 0 | 9.719 | prioritize |
| CHEMBL4789312 | 9.030 | low | 0.523 | 0 | 9.553 | prioritize |
| CHEMBL2031298 | 9.000 | low | 0.552 | 0 | 9.551 | prioritize |
| CHEMBL5269375 | 9.113 | low | 0.411 | 0 | 9.524 | prioritize |
| CHEMBL53753 | 8.689 | low | 0.769 | 0 | 9.458 | prioritize |
| CHEMBL5956184 | 9.042 | low | 0.411 | 0 | 9.454 | prioritize |
| CHEMBL2031300 | 8.838 | low | 0.574 | 0 | 9.411 | prioritize |
| CHEMBL6026192 | 9.265 | low | 0.388 | 1 | 9.404 | prioritize |
| CHEMBL4164805 | 9.000 | low | 0.386 | 0 | 9.386 | prioritize |
| CHEMBL54400 | 8.598 | low | 0.769 | 0 | 9.367 | prioritize |
| CHEMBL516022 | 8.795 | low | 0.818 | 0 | 9.363 | prioritize |
| CHEMBL6032856 | 9.163 | low | 0.406 | 1 | 9.320 | prioritize |
| CHEMBL5087955 | 9.090 | medium | 0.765 | 0 | 9.304 | prioritize |
| CHEMBL128987 | 8.634 | low | 0.623 | 0 | 9.257 | prioritize |
| CHEMBL3355873 | 9.193 | low | 0.306 | 0 | 9.249 | prioritize |
| CHEMBL162622 | 8.600 | low | 0.623 | 0 | 9.223 | prioritize |

Detailed ranked output is saved to `reports/egfr_ranked_existing_molecules.csv`.
