# EGFR Hardening Inventory

This inventory lists artifact paths, shapes, columns, missing-value counts, and metadata availability without exposing raw molecule structures.

## Metadata Availability

- model_ready_found: True
- pIC50_label_found: True
- molecule_identifier_found: True
- molecule_representation_found: True
- assay_metadata_available: True
- document_metadata_available: True
- assay_metadata_source: data/raw/egfr_chembl_ic50_raw.csv
- document_metadata_source: data/raw/egfr_chembl_ic50_raw.csv

## Artifact Summary

- `data/processed/egfr_model_ready.csv`: exists=True, bytes=1817562, shape=[10593, 14]
  columns: molecule_chembl_id, canonical_smiles, median_pIC50, median_IC50_nM, n_measurements, MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds, RingCount, HeavyAtomCount, QED
- `data/processed/egfr_ic50_clean.csv`: exists=True, bytes=1073957, shape=[10834, 5]
  columns: molecule_chembl_id, canonical_smiles, median_pIC50, median_IC50_nM, n_measurements
- `data/processed/egfr_standardized_molecules.csv`: exists=True, bytes=3214570, shape=[10593, 19]
  columns: molecule_chembl_id, canonical_smiles, median_pIC50, median_IC50_nM, n_measurements, MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds, RingCount, HeavyAtomCount, QED, standardized_smiles, molecule_hash, standardization_valid, scaffold, scaffold_hash
- `data/processed/features_morgan_index.csv`: exists=True, bytes=1051307, shape=[10593, 5]
  columns: molecule_chembl_id, median_pIC50, molecule_hash, scaffold_hash, scaffold
- `data/processed/features_morgan_fingerprints.npz`: exists=True, bytes=520422, shape=[10593, 2048]
- `data/processed/features_rdkit_descriptors.csv`: exists=True, bytes=1525073, shape=[10593, 13]
  columns: molecule_chembl_id, molecule_hash, scaffold_hash, median_pIC50, MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds, RingCount, HeavyAtomCount, QED
- `data/raw/egfr_chembl_ic50_raw.csv`: exists=True, bytes=3381778, shape=[26600, 9]
  columns: molecule_chembl_id, canonical_smiles, standard_type, standard_relation, standard_value, standard_units, assay_chembl_id, document_chembl_id, target_chembl_id
- `reports/metrics/qsar_matched_benchmark_metrics.json`: exists=True, bytes=5493, shape=None
- `reports/metrics/applicability_domain_metrics.json`: exists=True, bytes=1060, shape=None
- `reports/metrics/egfr_gnn_benchmark_metrics.json`: exists=True, bytes=1850, shape=None
- `reports/metrics/egfr_redocking_metrics.json`: exists=True, bytes=1041, shape=None
- `reports/egfr_ranked_existing_molecules.csv`: exists=True, bytes=3438171, shape=[10593, 30]
  columns: molecule_chembl_id, molecule_hash, scaffold_hash, predicted_pIC50, max_tanimoto_to_train, similarity_bin, out_of_domain_flag, rf_prediction_std, interval_lower_90, interval_upper_90, model_risk_category, QED, MolWt, MolLogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds, lipinski_violations, pains_alert, brenk_alert, synthetic_accessibility_status, synthetic_accessibility_score, property_penalty, model_risk_penalty, uncertainty_penalty, final_score, final_triage_category, true_pIC50, absolute_error
- `data/structure_prepared/5UG9_receptor.pdbqt`: exists=True, bytes=198589, shape=None
- `data/structure_prepared/5UG9_8AM_ligand.pdbqt`: exists=True, bytes=3387, shape=None
- `data/structure_prepared/5UG9_8AM_redocked_out.pdbqt`: exists=True, bytes=28848, shape=None
