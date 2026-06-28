# EGFR Hardening Fix Report

Branch: `portfolio-hardening-final`

## Summary

- Reconciled public final status to `DONE`.
- Rephrased uncertainty as a conformal-style retrospective uncertainty check where appropriate.
- Clarified molecular standardization policy when tautomer canonicalization is skipped.
- Framed the custom PyTorch GCN as exploratory negative benchmark evidence.
- Framed redocking as a retrospective Vina pose-recovery audit.
- Added reproducibility and CI entry points.
- Added tests for status consistency and conservative public wording.

## Checks Run

- `make reproduce-small PYTHON=/usr/bin/python3`: pass.
- `make test PYTHON=/usr/bin/python3`: pass, 11 tests.
- Public report CSV SMILES-header scan: pass, zero SMILES headers.
- Unsupported-claim keyword scan: pass for public-facing claim terms.
- `git diff --check`: pass.

## Remaining Manual Review

- Full docking and structure workflows depend on local structure-preparation tooling; the lightweight check verifies committed reports and metrics only.
