"""Preflight checks for EGFR hardening stages."""

from __future__ import annotations

import importlib
import json
import platform
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPORTS_DIR = PROJECT_ROOT / "reports"
METRICS_DIR = REPORTS_DIR / "metrics"


def module_status(name: str) -> dict:
    """Return import status and version where available."""
    try:
        module = importlib.import_module(name)
        version = getattr(module, "__version__", "available")
        return {"available": True, "version": str(version)}
    except Exception as exc:
        return {"available": False, "error": exc.__class__.__name__}


def main() -> None:
    """Run preflight and save metrics."""
    for relative in [
        "reports",
        "reports/metrics",
        "reports/figures",
        "reports/run_logs",
        "scripts",
        "tests",
        "src/validation",
        "src/analysis",
        "src/structure",
    ]:
        (PROJECT_ROOT / relative).mkdir(parents=True, exist_ok=True)

    packages = {
        name: module_status(name)
        for name in ["pandas", "numpy", "sklearn", "scipy", "matplotlib", "joblib", "rdkit"]
    }
    artifacts = {
        "model_ready_data": "data/processed/egfr_model_ready.csv",
        "morgan_features": "data/processed/features_morgan_fingerprints.npz",
        "morgan_index": "data/processed/features_morgan_index.csv",
        "final_report": "reports/final_egfr_cadd_qsar_report.md",
        "gnn_metrics": "reports/metrics/egfr_gnn_benchmark_metrics.json",
        "redocking_metrics": "reports/metrics/egfr_redocking_metrics.json",
        "redocked_pose": "data/structure_prepared/5UG9_8AM_redocked_out.pdbqt",
    }
    artifact_status = {
        key: {
            "path": path,
            "exists": (PROJECT_ROOT / path).exists(),
            "bytes": (PROJECT_ROOT / path).stat().st_size if (PROJECT_ROOT / path).exists() else 0,
        }
        for key, path in artifacts.items()
    }
    payload = {
        "status": "completed" if artifact_status["model_ready_data"]["exists"] else "blocked_model_ready_missing",
        "project_root": str(PROJECT_ROOT),
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "packages": packages,
        "artifacts": artifact_status,
    }
    (METRICS_DIR / "egfr_hardening_preflight.json").write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Preflight status: {payload['status']}")


if __name__ == "__main__":
    main()
