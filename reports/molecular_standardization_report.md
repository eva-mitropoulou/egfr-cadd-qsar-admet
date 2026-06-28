# Molecular Standardization Report

- Input molecules: 10,593
- Standardized molecules: 10,593
- Invalid molecules: 0
- Duplicate standardized molecules: 42
- Unique standardized molecules: 10,551
- RDKit MolStandardize available: True

## Policies

- Fragment/salt handling: largest fragment chosen when RDKit MolStandardize is available; otherwise sanitized canonical molecule retained
- Charge handling: RDKit Uncharger applied when available; otherwise original sanitized charge state retained
- Tautomer/protonation handling: RDKit sanitization, largest-fragment handling, charge normalization where available, and documented tautomer/protonation policy. Full tautomer canonicalization was skipped in this automated finish run.
- Stereochemistry handling: isomeric canonical SMILES retained; stereochemistry is not collapsed
- Public reporting: reports use molecule IDs, hashes, and aggregate counts rather than raw SMILES

Detailed standardized molecule representations are saved to `data/processed/egfr_standardized_molecules.csv`.

## Representation Audit

- Rows audited: 10,593
- Missing standardized representations: 0
- Duplicate molecule hashes: 42
- Unique scaffold hashes: 3,685
- Nonzero formal charge count: 13
