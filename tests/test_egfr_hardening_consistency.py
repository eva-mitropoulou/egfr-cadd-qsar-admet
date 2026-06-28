import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read_json(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def test_final_status_matches_hardening_metrics():
    report = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    hardening = read_json("reports/metrics/egfr_final_hardening_status.json")
    expected = "DONE_WITH_WARNINGS" if hardening.get("degraded_items") else "DONE"
    assert f"FINAL_STATUS = {expected}" in report


def test_uncertainty_wording_is_conservative():
    report = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    proxy = read_json("reports/metrics/egfr_uncertainty_calibration_metrics.json")
    if "degraded_proxy" in proxy.get("conformal_status", ""):
        assert "Conformal-Style Uncertainty Check" in report
        assert "retrospective uncertainty proxy" in report


def test_gnn_wording_is_exploratory_negative_benchmark():
    report = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    assert "Exploratory Custom PyTorch GCN Baseline" in report
    assert "negative benchmark evidence" in report


def test_standardization_wording_does_not_imply_full_tautomer_resolution():
    report = (ROOT / "reports/molecular_standardization_report.md").read_text(encoding="utf-8")
    metrics = read_json("reports/metrics/molecular_standardization_metrics.json")
    if metrics.get("tautomer_canonicalized_count") == 0:
        assert "Full tautomer canonicalization was skipped" in report
        assert "documented tautomer/protonation policy" in report


def test_redocking_wording_is_pose_recovery_audit():
    report = (ROOT / "reports/final_egfr_cadd_qsar_report.md").read_text(encoding="utf-8")
    assert "retrospective Vina redocking pose-recovery audit" in report
    assert "not a binding free-energy calculation" in report
