# EGFR Uncertainty and Calibration Report

Uncertainty was estimated from an applicability-domain reliability proxy using existing out-of-fold scaffold-CV predictions.

| uncertainty_bin | count | mean_uncertainty | MAE | coverage_90 |
| --- | --- | --- | --- | --- |
| low | 2656 | 0.131 | 0.436 | 0.960 |
| medium_low | 2641 | 0.230 | 0.543 | 0.930 |
| medium_high | 2653 | 0.314 | 0.677 | 0.886 |
| high | 2643 | 0.492 | 0.816 | 0.824 |

- Spearman correlation between uncertainty proxy and absolute error: 0.273
- Empirical 90 percent interval coverage: 0.900
- Interval status: degraded_proxy_interval_from_out_of_fold_residual_quantile
