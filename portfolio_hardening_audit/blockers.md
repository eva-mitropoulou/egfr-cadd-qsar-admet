# Hardening Blockers

## EGFR CI Workflow

- Repo: `egfr-cadd-qsar-admet`
- File: `.github/workflows/egfr-cadd.yml`
- Gate: CI workflow
- Reason: GitHub rejected the push because the available token lacks `workflow` scope for creating or updating workflow files.
- Recommended manual fix: add the workflow manually in GitHub or push it with a token that has `workflow` scope. The committed `Makefile`, `pyproject.toml`, and `scripts/verify_public_artifacts.py` contain the local commands the workflow should run.
