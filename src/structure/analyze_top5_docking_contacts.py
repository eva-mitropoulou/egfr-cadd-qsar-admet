"""Deprecated: top-5 contact-overlap analysis was removed from the workflow.

The final top-5 structure layer reports ligand preparation, docking status, and
Vina scores only. The separate 5UG9/8AM redocking validation remains in place.
"""

from __future__ import annotations


def main() -> None:
    """Report the deprecation status without generating contact-overlap outputs."""
    print("Top-5 contact-overlap analysis: removed")
    print("Use src/structure/build_top5_structure_sanity_report.py for Vina-score summaries.")


if __name__ == "__main__":
    main()
