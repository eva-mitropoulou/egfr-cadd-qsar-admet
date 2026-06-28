# Initial Hardening Audit: egfr-cadd-qsar-admet

- Default branch: `main`
- Hardening branch: `portfolio-hardening-final`
- Start commit: `4b46b38a8a16`
- Tracked files: 202
- Markdown files: 34
- Metrics JSON files: 26
- Reports: 31
- Tests detected: 4
- CI workflows: 0

## Public-Facing Surfaces

- README, docs, projects, portfolio assets, reports, model/data cards where present.
- Public claim scan hits for manual review: 3
- Raw sequence-like public files found: 0

## Immediate Fixes Needed

- Reconcile DONE versus DONE_WITH_WARNINGS status.
- Use conformal-style uncertainty wording for degraded proxy metrics.
- Conservatively word standardization, GNN, and redocking evidence.
- Add Makefile, pyproject/environment if missing, and CI entry point.

## Reproducibility Files

- Makefile present: False
- pyproject present: False
- requirements present: True
- environment.yml present: False
- CI workflows: none

## Notes

This audit records file and claim-scan metadata only. It does not include raw rows, raw sequences, or raw molecule tables.
