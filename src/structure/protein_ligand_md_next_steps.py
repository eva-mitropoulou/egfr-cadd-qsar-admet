"""Document optional EGFR protein-ligand MD next steps."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from egfr_pipeline_utils import REPORTS_DIR, write_text  # noqa: E402


def main() -> None:
    """Create a non-blocking future-work MD bridge report."""
    report = [
        "# EGFR Protein-Ligand MD Next Steps",
        "",
        "Full protein-ligand MD was not run automatically in this retrospective QSAR and CADD workflow.",
        "",
        "A credible EGFR-ligand MD extension would require:",
        "",
        "- a prepared EGFR kinase-domain structure with validated ligand pose",
        "- ligand protonation/tautomer state selection",
        "- force-field parameters for protein and ligand",
        "- solvated/neutralized simulation system",
        "- minimization, equilibration, and production MD inputs",
        "",
        "Suggested analyses:",
        "",
        "- protein backbone RMSD",
        "- ligand heavy-atom RMSD",
        "- binding-site RMSF",
        "- ligand-contact persistence",
        "- hydrogen-bond occupancy",
        "- interaction-fingerprint persistence",
        "",
        "This is future work for extending the current project evidence.",
        "",
    ]
    write_text(REPORTS_DIR / "egfr_protein_ligand_md_next_steps.md", "\n".join(report))
    print(f"MD next-steps report: {REPORTS_DIR / 'egfr_protein_ligand_md_next_steps.md'}")


if __name__ == "__main__":
    main()
