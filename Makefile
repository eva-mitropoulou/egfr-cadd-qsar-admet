PYTHON ?= python

.PHONY: reproduce-small test report figures

reproduce-small:
	$(PYTHON) scripts/verify_public_artifacts.py

test:
	$(PYTHON) -m pytest -q

report:
	$(PYTHON) src/analysis/build_final_hardening_status.py

figures:
	$(PYTHON) scripts/verify_public_artifacts.py --figures-only
