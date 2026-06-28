# EGFR ChEMBL Data Provenance Audit

- Source: public ChEMBL activity records
- Target IDs found: CHEMBL203
- Raw activity rows: 26,600
- Raw unique molecules: 14,776
- Assays: 2,425
- Documents: 1,298
- IC50 rows: 26,600
- Exact `=` rows: 19,719
- nM rows: 25,244
- Clean molecule-level pIC50 rows: 10,834
- Model-ready molecule rows: 10,593

## Cleaning Policy

- pIC50 = 9 - log10(IC50_nM), using exact nM IC50 values only for the clean regression target.
- Duplicate activity measurements are aggregated to one molecule-level target using median pIC50 and median IC50_nM.

## Distribution Summaries

- Standard type counts: {'IC50': 26600}
- Standard units counts: {'nM': 25244, 'missing': 1259, 'ug.mL-1': 89, '/uM': 6, '10^3 uM': 1, '10^-5 mol/L': 1}
- Standard relation counts: {'=': 19719, '>': 3309, '<': 2067, 'missing': 1398, '<=': 87, '>>': 11, '>=': 8, '~': 1}
