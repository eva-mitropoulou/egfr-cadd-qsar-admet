"""Optional EGFR redocking validation status report."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import METRICS_DIR, REPORTS_DIR, read_json, save_json, write_text  # noqa: E402


def main() -> None:
    """Report redocking status without making docking a hard dependency."""
    metrics_path = METRICS_DIR / "egfr_structure_module_metrics.json"
    metrics = read_json(metrics_path)
    vina_available = shutil.which("vina") is not None or shutil.which("autodock_vina") is not None
    structures_available = int(metrics.get("available_structures", 0)) > 0

    if vina_available and structures_available:
        status = "ready_not_run_automatically"
        reason = "Docking executable and at least one structure are available; automated redocking was not launched by this retrospective pipeline."
    else:
        status = "skipped"
        reason = "Docking skipped because a docking executable and/or prepared local structure files were unavailable."

    metrics.update(
        {
            "redocking_status": status,
            "redocking_reason": reason,
            "vina_available": vina_available,
            "pose_recovery_metric": None,
        }
    )
    save_json(metrics_path, metrics)

    report = [
        "# EGFR Redocking Validation Report",
        "",
        f"- Redocking status: {status}",
        f"- Vina-like executable available: {vina_available}",
        f"- Pose recovery metric: {metrics['pose_recovery_metric']}",
        f"- Reason: {reason}",
        "",
        "Docking is treated as optional structure-based evidence and does not affect QSAR benchmark completion.",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_redocking_report.md", "\n".join(report))

    print(f"Redocking status: {status}")


if __name__ == "__main__":
    main()
