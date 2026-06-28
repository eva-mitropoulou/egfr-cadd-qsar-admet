"""Optional EGFR interaction-fingerprint status report."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import FIGURES_DIR, METRICS_DIR, REPORTS_DIR, read_json, save_figure, save_json, setup_matplotlib, write_text  # noqa: E402


setup_matplotlib()
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    """Create interaction-fingerprint report or clear degraded status."""
    metrics_path = METRICS_DIR / "egfr_structure_module_metrics.json"
    metrics = read_json(metrics_path)
    structures_available = int(metrics.get("available_structures", 0)) > 0

    if structures_available:
        status = "metadata_only"
        reason = "Structure files were available, but ligand extraction and PLIP-style interaction generation were not required for this automated run."
    else:
        status = "skipped_no_structures"
        reason = "No prepared local structures were available for interaction fingerprinting."

    metrics.update(
        {
            "interaction_fingerprint_status": status,
            "interaction_fingerprint_reason": reason,
            "interaction_count": 0,
        }
    )
    save_json(metrics_path, metrics)

    plt.figure(figsize=(6, 4))
    plt.bar(["interaction_fingerprint"], [0], color="#BAB0AC")
    plt.ylabel("Interactions counted")
    plt.title("Interaction Fingerprint Status")
    save_figure(FIGURES_DIR / "interaction_frequency.png")

    report = [
        "# EGFR Interaction Fingerprint Report",
        "",
        f"- Interaction-fingerprint status: {status}",
        f"- Interaction count: {metrics['interaction_count']}",
        f"- Reason: {reason}",
        "",
        "This module is documented as optional structure-based analysis and is not used to claim binding efficacy.",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_interaction_fingerprint_report.md", "\n".join(report))

    print(f"Interaction fingerprint status: {status}")


if __name__ == "__main__":
    main()
