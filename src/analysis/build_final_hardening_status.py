"""Build final EGFR hardening status files."""

from __future__ import annotations

import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REPORTS_DIR = PROJECT_ROOT / "reports"
METRICS_DIR = REPORTS_DIR / "metrics"
STATUS_REPORT = REPORTS_DIR / "egfr_final_hardening_status.md"
STATUS_JSON = METRICS_DIR / "egfr_final_hardening_status.json"


def read_json(path: Path) -> dict:
    """Read JSON with empty fallback."""
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    """Create final hardening status summary."""
    assay = read_json(METRICS_DIR / "egfr_assay_aware_validation_metrics.json")
    conformal = read_json(METRICS_DIR / "egfr_conformal_uncertainty_metrics.json")
    sar = read_json(METRICS_DIR / "egfr_sar_interpretability_metrics.json")
    redocking = read_json(METRICS_DIR / "egfr_redocking_audit_metrics.json")
    smoke = read_json(METRICS_DIR / "egfr_smoke_test_metrics.json")
    inventory = read_json(METRICS_DIR / "egfr_hardening_inventory.json")

    assay_rows = assay.get("validation_rows", [])
    assay_split = next((row for row in assay_rows if row.get("split") == "assay_group_split"), {})
    document_split = next((row for row in assay_rows if row.get("split") == "document_group_split"), {})
    degraded = []
    if assay_split.get("status") != "completed":
        degraded.append({"stage": "assay-aware validation", "reason": assay_split.get("status")})
    if document_split.get("status") != "completed":
        degraded.append({"stage": "document-aware validation", "reason": document_split.get("status")})
    if redocking.get("overlay_artifact_status") != "overlay_figure_created":
        degraded.append({"stage": "redocking overlay", "reason": redocking.get("overlay_artifact_status")})
    if smoke.get("status") != "passed":
        degraded.append({"stage": "smoke tests", "reason": smoke.get("status")})

    blocked = not inventory.get("metadata_availability", {}).get("model_ready_found", True)
    final_status = "BLOCKED" if blocked else ("DONE_WITH_WARNINGS" if degraded else "DONE")
    payload = {
        "FINAL_HARDENING_STATUS": final_status,
        "assay_aware_validation_status": assay_split.get("status"),
        "document_aware_validation_status": document_split.get("status"),
        "conformal_uncertainty_status": conformal.get("status"),
        "SAR_analysis_status": sar.get("status"),
        "redocking_audit_status": redocking.get("status"),
        "overlay_artifact_status": redocking.get("overlay_artifact_status"),
        "smoke_test_status": smoke.get("status"),
        "degraded_items": degraded,
    }
    STATUS_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    lines = [
        "# EGFR Final Hardening Status",
        "",
        f"- FINAL_HARDENING_STATUS: {final_status}",
        f"- Assay-aware validation status: {payload['assay_aware_validation_status']}",
        f"- Document-aware validation status: {payload['document_aware_validation_status']}",
        f"- Conformal uncertainty status: {payload['conformal_uncertainty_status']}",
        f"- SAR analysis status: {payload['SAR_analysis_status']}",
        f"- Redocking audit status: {payload['redocking_audit_status']}",
        f"- Overlay artifact status: {payload['overlay_artifact_status']}",
        f"- Smoke-test status: {payload['smoke_test_status']}",
        "",
    ]
    if degraded:
        lines.append("## Degraded Items")
        lines.append("")
        for item in degraded:
            lines.append(f"- {item['stage']}: {item['reason']}")
    STATUS_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"FINAL_HARDENING_STATUS: {final_status}")


if __name__ == "__main__":
    main()
