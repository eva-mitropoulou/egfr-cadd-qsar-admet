# EGFR Assay-Aware and Document-Aware Validation

This stage tests whether the Morgan fingerprint Random Forest remains useful when held-out groups are defined by assay or document context.

## Metadata Availability

- Raw activity table available: True
- Assay metadata available: True
- Document metadata available: True
- Metadata mapping status: mapped_from_raw_activity
- Molecules with assay metadata: 10593
- Molecules with document metadata: 10593

## Validation Summary

### assay_group_split

- Status: completed
- Train molecules: 8570
- Test molecules: 2023
- Train unique groups: 867
- Test unique groups: 217
- Group overlap count: 0
- Molecule overlap count: 0
- MAE: 0.7959778908221339
- RMSE: 1.0135992949202917
- R2: 0.4480171002230632
- Pearson: 0.6705173336586066
- Spearman: 0.6870709476123784

### document_group_split

- Status: completed
- Train molecules: 8646
- Test molecules: 1947
- Train unique groups: 724
- Test unique groups: 182
- Group overlap count: 0
- Molecule overlap count: 0
- MAE: 0.8814069616241912
- RMSE: 1.1425119156890893
- R2: 0.21240743631437198
- Pearson: 0.5042280029635546
- Spearman: 0.4893740873666689

No raw molecule structures or raw assay records are shown in this report.