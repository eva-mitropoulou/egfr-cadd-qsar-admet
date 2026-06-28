"""Preflight checks for the EGFR CADD/QSAR decision workflow."""

from __future__ import annotations

import importlib.util
import json
import platform
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import (  # noqa: E402
    DATA_DIR,
    FIGURES_DIR,
    METRICS_DIR,
    MODELS_DIR,
    PROCESSED_DIR,
    RAW_DIR,
    REPORTS_DIR,
    RUN_LOGS_DIR,
    ensure_project_dirs,
    save_json,
    write_text,
)


ESSENTIAL_PACKAGES = [
    "pandas",
    "numpy",
    "scipy",
    "sklearn",
    "matplotlib",
    "rdkit",
    "joblib",
    "requests",
]

OPTIONAL_PACKAGES = [
    "tqdm",
    "xgboost",
    "torch",
    "torch_geometric",
    "deepchem",
    "chemprop",
    "streamlit",
    "mordred",
    "umap",
    "shap",
]


def package_available(name: str) -> bool:
    """Return whether a package import can be resolved."""
    return importlib.util.find_spec(name) is not None


def maybe_install_missing_essential(packages: list[str]) -> list[str]:
    """Try installing missing essential packages, without using sudo."""
    attempted: list[str] = []
    pip_names = {"sklearn": "scikit-learn", "rdkit": "rdkit-pypi"}
    for package in packages:
        if package_available(package):
            continue
        pip_name = pip_names.get(package, package)
        attempted.append(pip_name)
        subprocess.run(
            [sys.executable, "-m", "pip", "install", pip_name],
            cwd=PROJECT_ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    return attempted


def main() -> None:
    """Run preflight checks and save a machine-readable status file."""
    ensure_project_dirs()

    missing_essential_before = [pkg for pkg in ESSENTIAL_PACKAGES if not package_available(pkg)]
    install_attempts = maybe_install_missing_essential(missing_essential_before)
    missing_essential_after = [pkg for pkg in ESSENTIAL_PACKAGES if not package_available(pkg)]
    optional_status = {pkg: package_available(pkg) for pkg in OPTIONAL_PACKAGES}

    gpu_status = {"cuda_available": False, "gpu_detail": "torch not available"}
    if package_available("torch"):
        try:
            import torch

            gpu_status = {
                "cuda_available": bool(torch.cuda.is_available()),
                "gpu_detail": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "no CUDA device",
            }
        except Exception as exc:  # pragma: no cover - defensive reporting
            gpu_status = {"cuda_available": False, "gpu_detail": f"torch check failed: {type(exc).__name__}"}

    checked_paths = {
        "project_root": PROJECT_ROOT.exists(),
        "raw_data_dir": RAW_DIR.exists(),
        "processed_data_dir": PROCESSED_DIR.exists(),
        "reports_dir": REPORTS_DIR.exists(),
        "figures_dir": FIGURES_DIR.exists(),
        "models_dir": MODELS_DIR.exists(),
        "run_logs_dir": RUN_LOGS_DIR.exists(),
    }

    existing_artifact_counts = {
        "raw_data_files": len(list(RAW_DIR.glob("*"))) if RAW_DIR.exists() else 0,
        "processed_data_files": len(list(PROCESSED_DIR.glob("*"))) if PROCESSED_DIR.exists() else 0,
        "report_files": len(list(REPORTS_DIR.glob("**/*"))) if REPORTS_DIR.exists() else 0,
        "figure_files": len(list(FIGURES_DIR.glob("*"))) if FIGURES_DIR.exists() else 0,
        "model_files": len(list(MODELS_DIR.glob("*"))) if MODELS_DIR.exists() else 0,
    }

    metrics = {
        "python_executable": sys.executable,
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "essential_packages": {pkg: package_available(pkg) for pkg in ESSENTIAL_PACKAGES},
        "missing_essential_before_install": missing_essential_before,
        "missing_essential_after_install": missing_essential_after,
        "install_attempts": install_attempts,
        "optional_packages": optional_status,
        "gpu": gpu_status,
        "checked_paths": checked_paths,
        "existing_artifact_counts": existing_artifact_counts,
        "preflight_pass": len(missing_essential_after) == 0 and checked_paths["project_root"],
    }

    save_json(METRICS_DIR / "preflight_metrics.json", metrics)
    report = [
        "# EGFR Project Preflight",
        "",
        f"- Project root exists: {checked_paths['project_root']}",
        f"- Python executable: `{sys.executable}`",
        f"- Python version: {metrics['python_version']}",
        f"- Missing essential packages after install attempts: {missing_essential_after}",
        f"- CUDA available: {gpu_status['cuda_available']}",
        f"- Raw data files: {existing_artifact_counts['raw_data_files']}",
        f"- Processed data files: {existing_artifact_counts['processed_data_files']}",
        f"- Existing model files: {existing_artifact_counts['model_files']}",
        "",
        "Optional packages are recorded in `reports/metrics/preflight_metrics.json`.",
        "",
    ]
    write_text(REPORTS_DIR / "preflight_report.md", "\n".join(report))

    print(json.dumps({"preflight_pass": metrics["preflight_pass"], "metrics_path": str(METRICS_DIR / "preflight_metrics.json")}))


if __name__ == "__main__":
    main()
