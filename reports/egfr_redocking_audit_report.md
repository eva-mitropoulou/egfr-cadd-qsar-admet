# EGFR Redocking Audit Report

This audit hardens the completed retrospective Vina redocking evidence for the 5UG9 / 8AM co-crystal case.

## Result

- Redocking status: completed_redocking
- Docking score: -9.471 kcal/mol
- Pose recovery RMSD: 0.968 angstrom
- Pose recovery status: successful

## Tooling

- Vina Python package: 1.2.7
- Vina CLI: unavailable
- Meeko: 0.7.1
- OpenBabel CLI: Open Babel 3.1.0 -- Feb 17 2026 -- 06:47:47
- RDKit: 2026.03.3

## Input/Output Audit

- Receptor PDBQT: `data/structure_prepared/5UG9_receptor.pdbqt` (198589 bytes)
- Ligand PDBQT: `data/structure_prepared/5UG9_8AM_ligand.pdbqt` (3387 bytes)
- Docked pose: `data/structure_prepared/5UG9_8AM_redocked_out.pdbqt` (28848 bytes)
- Overlay artifact status: overlay_figure_created
- Overlay figure: `reports/figures/5UG9_8AM_redocking_pose_overlay.png`
- Overlay script: `reports/structure_visualization/5UG9_8AM_overlay.pml`

## RMSD Policy

- Method: direct fixed-frame coordinate RMSD from prepared ligand PDBQT original coordinates to first Vina pose
- Atom scope: heavy atom
- Hydrogens ignored: True
- Atom mapping: direct prepared-ligand PDBQT atom order

This is retrospective redocking validation of a known co-crystallized ligand, not a binding free-energy calculation or prospective discovery claim.
