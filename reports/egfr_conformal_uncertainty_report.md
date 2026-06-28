# EGFR Split-Conformal Uncertainty Report

This stage adds formal split-conformal pIC50 prediction intervals using absolute calibration residuals.

## random_split_conformal

- Status: completed
- Target coverage: 0.90
- Empirical coverage: 0.902
- Mean interval width: 2.441
- Median interval width: 2.441
- Test MAE: 0.538
- Test RMSE: 0.729
- Test R2: 0.705

## scaffold_group_conformal

- Status: completed
- Target coverage: 0.90
- Empirical coverage: 0.931
- Mean interval width: 3.094
- Median interval width: 3.094
- Test MAE: 0.642
- Test RMSE: 0.846
- Test R2: 0.579

Intervals quantify retrospective model uncertainty on held-out EGFR records; they are not clinical confidence statements.