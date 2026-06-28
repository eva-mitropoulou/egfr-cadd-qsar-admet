# EGFR Retrospective Active-Learning Report

This simulation uses existing labeled ChEMBL/project records only. No new molecules are generated.

- Molecules available: 10,593
- Initial labeled seed size: 400
- Batch size: 500
- Rounds: 5
- Model: Ridge regression on RDKit descriptor features
- Potent threshold: top 10 percent, pIC50 >= 8.610
- Best strategy: applicability_domain_aware_high_score

## Final Strategy Summary

| strategy | labeled_count | top_potent_recovery_fraction | selected_scaffold_count | mean_selected_pIC50 |
| --- | --- | --- | --- | --- |
| highest_predicted_pIC50 | 2900 | 0.366 | 1304 | 7.408 |
| random | 2900 | 0.246 | 1511 | 6.883 |
| applicability_domain_aware_high_score | 2900 | 0.400 | 1360 | 7.478 |
| scaffold_diverse_high_score | 2900 | 0.328 | 1787 | 7.286 |
| uncertainty_sampling | 2900 | 0.258 | 1438 | 6.723 |
| hybrid_activity_uncertainty_diversity | 2900 | 0.309 | 1827 | 7.217 |
