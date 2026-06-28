#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

PYTHON_BIN="${PYTHON_BIN:-python}"
if [ ! -x "$PYTHON_BIN" ]; then
  PYTHON_BIN="${PYTHON:-python}"
fi

"$PYTHON_BIN" src/analysis/egfr_hardening_inventory.py
"$PYTHON_BIN" src/structure/harden_redocking_evidence.py
"$PYTHON_BIN" src/analysis/patch_final_egfr_reports.py
"$PYTHON_BIN" src/analysis/build_final_hardening_status.py

echo "Reproduced final EGFR hardening reports from existing artifacts."
