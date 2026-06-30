PYTHON ?= python

.PHONY: reproduce-small test report figures

reproduce-small: test

test:
	$(PYTHON) -m pytest -q

report:
	$(PYTHON) src/analysis/build_final_egfr_report.py

figures:
	$(PYTHON) src/analysis/build_readme_figures.py
