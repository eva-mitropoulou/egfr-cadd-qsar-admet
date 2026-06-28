from __future__ import annotations

import argparse
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED_REPORTS = [
    "README.md",
    "reports/final_egfr_cadd_qsar_report.md",
    "reports/egfr_final_hardening_status.md",
    "reports/egfr_conformal_uncertainty_report.md",
    "reports/egfr_gnn_benchmark_report.md",
    "reports/egfr_redocking_audit_report.md",
]

REQUIRED_METRICS = [
    "reports/metrics/egfr_final_hardening_status.json",
    "reports/metrics/egfr_uncertainty_calibration_metrics.json",
    "reports/metrics/egfr_conformal_uncertainty_metrics.json",
    "reports/metrics/molecular_standardization_metrics.json",
    "reports/metrics/egfr_gnn_benchmark_metrics.json",
    "reports/metrics/egfr_redocking_audit_metrics.json",
]

REQUIRED_FIGURES = [
    "reports/figures/random_vs_scaffold_vs_assay_document_split.png",
    "reports/figures/conformal_coverage_by_split.png",
    "reports/figures/5UG9_8AM_redocking_pose_overlay.png",
]


def assert_nonempty(paths: list[str]) -> None:
    missing = [path for path in paths if not (ROOT / path).is_file() or (ROOT / path).stat().st_size == 0]
    if missing:
        raise SystemExit(f"Missing required artifacts: {missing}")


def parse_json(paths: list[str]) -> None:
    for path in paths:
        payload = json.loads((ROOT / path).read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise SystemExit(f"Expected JSON object: {path}")


def compile_project_python() -> None:
    for folder in ["src", "scripts"]:
        for path in (ROOT / folder).rglob("*.py"):
            py_compile.compile(str(path), doraise=True)


def check_status_consistency() -> None:
    final_report = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    hardening = json.loads((ROOT / "reports/metrics/egfr_final_hardening_status.json").read_text(encoding="utf-8"))
    degraded = hardening.get("degraded_items") or []
    expected = "DONE_WITH_WARNINGS" if degraded else "DONE"
    if f"FINAL_STATUS = {expected}" not in final_report:
        raise SystemExit(f"Final report status does not match hardening status: expected {expected}")


def check_public_wording() -> None:
    report = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    required = [
        "Conformal-Style Uncertainty Check",
        "Exploratory Custom PyTorch GCN Baseline",
        "retrospective Vina redocking pose-recovery audit",
        "FINAL_STATUS = DONE",
    ]
    missing = [text for text in required if text not in report]
    if missing:
        raise SystemExit(f"Missing public wording: {missing}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures-only", action="store_true")
    args = parser.parse_args()
    if args.figures_only:
        assert_nonempty(REQUIRED_FIGURES)
        return 0
    assert_nonempty(REQUIRED_REPORTS + REQUIRED_FIGURES)
    parse_json(REQUIRED_METRICS)
    compile_project_python()
    check_status_consistency()
    check_public_wording()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
