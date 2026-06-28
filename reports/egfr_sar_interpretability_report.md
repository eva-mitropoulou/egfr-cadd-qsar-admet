# EGFR SAR-Support and Interpretable Error Analysis

This report adds model interpretation and error-analysis evidence for the retrospective EGFR QSAR workflow.
The analysis is SAR-supportive and descriptive; it does not claim causal binding mechanisms or prospective discovery.

## Descriptor Importance

- MolWt: 0.1784
- TPSA: 0.1623
- QED: 0.1388
- MolLogP: 0.1375
- RingCount: 0.0866
- NumRotatableBonds: 0.0804
- HeavyAtomCount: 0.0786
- NumHAcceptors: 0.0742
- NumHDonors: 0.0631

## Morgan Fingerprint Importance

Top Morgan fingerprint bits are reported as bit IDs only because reliable substructure reconstruction was not required for this evidence layer.
- bit 1452: 0.15896
- bit 228: 0.05182
- bit 489: 0.01783
- bit 1152: 0.01172
- bit 1645: 0.01039
- bit 1291: 0.01030
- bit 1435: 0.01006
- bit 687: 0.01000
- bit 935: 0.00962
- bit 378: 0.00792

## Activity Cliff Candidates

- Similarity threshold: >= 0.85
- Activity-difference threshold: >= 1.0 pIC50
- Activity cliff candidate pairs saved: 607
- Detailed table: `reports/egfr_activity_cliffs.csv`

## Scaffold-Level Error

- Count-filtered scaffold rows: 387
- Detailed table: `reports/egfr_scaffold_error_table.csv`

## Nearest-Neighbor Evidence for Top-Ranked Existing Molecules

- Top-ranked molecules summarized: 20
- Reports use molecule IDs/hashes and scaffold hashes only.
