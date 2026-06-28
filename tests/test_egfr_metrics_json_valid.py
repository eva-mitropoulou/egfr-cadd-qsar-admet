import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_metrics_json_files_parse():
    metrics_files = [
        "reports/metrics/egfr_hardening_inventory.json",
        "reports/metrics/egfr_assay_aware_validation_metrics.json",
        "reports/metrics/egfr_conformal_uncertainty_metrics.json",
        "reports/metrics/egfr_sar_interpretability_metrics.json",
        "reports/metrics/egfr_redocking_audit_metrics.json",
        "reports/metrics/egfr_final_hardening_status.json",
        "reports/metrics/agentic_egfr_hardening_state.json",
    ]
    for path in metrics_files:
        payload = json.loads((ROOT / path).read_text(encoding="utf-8"))
        assert isinstance(payload, dict)


def test_key_metrics_present_in_final_report():
    text = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    required_snippets = [
        "26,600",
        "10,593",
        "RMSE 0.871",
        "R2 0.550",
        "Assay/Document-Aware Validation",
        "Conformal-Style Uncertainty Check",
        "Pose recovery RMSD",
    ]
    for snippet in required_snippets:
        assert snippet in text
