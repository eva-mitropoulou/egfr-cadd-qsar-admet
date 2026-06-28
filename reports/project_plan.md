# Project Plan

## Scope

This project will build a reproducible CADD workflow for prioritizing small-molecule EGFR inhibitor candidates from ChEMBL bioactivity data. It will focus on data curation, RDKit descriptor and fingerprint generation, QSAR modeling, scaffold-split validation, ADMET-style filtering, and candidate ranking.

## Pipeline

1. Fetch EGFR target and bioactivity records from ChEMBL.
2. Filter records to relevant activity endpoints, preferably IC50 or comparable standardized activity measurements.
3. Clean compound records by checking SMILES validity, harmonizing activity units, resolving duplicates, and creating modeling targets.
4. Calculate RDKit descriptors and Morgan fingerprints.
5. Perform EDA on chemical properties and activity distributions.
6. Train baseline QSAR models using scikit-learn.
7. Evaluate models with random splits and scaffold splits to expose model-risk differences.
8. Apply simple ADMET-style filters and drug-likeness rules.
9. Rank compounds for follow-up based on predicted activity, validation risk, and property profile.
10. Summarize results in tables, plots, and short reports.

## Expected Deliverables

- Raw and processed data files.
- Reusable Python modules in `src/`.
- Four workflow notebooks.
- QSAR performance summaries.
- Ranked candidate table.
- Figures supporting EDA, validation, and ranking decisions.

## Current Status

Step 1 creates the repository structure, README, dependency list, notebook placeholders, and documented source-module stubs. Implementation starts in the next step with fetching and cleaning ChEMBL EGFR bioactivity data.
