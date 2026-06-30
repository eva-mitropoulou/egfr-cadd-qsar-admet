# External Unwanted-Substructure SMARTS Schema

Optional file name: `unwanted_substructures.csv`

Expected columns:

- `name`: human-readable alert name
- `smarts`: SMARTS pattern
- `description`: optional note

Accepted aliases include `alert`, `rule`, or `substructure` for the name field and `SMARTS` or `pattern` for the SMARTS field.

These SMARTS are used as medicinal-chemistry risk annotations and sensitivity-analysis filters, not automatic exclusions from the primary EGFR QSAR benchmark.
