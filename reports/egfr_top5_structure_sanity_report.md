# EGFR Top-5 Structure Sanity Report

## Purpose

This module links the final ranked EGFR molecule table to the validated 5UG9 structure workflow.
Docking of top-ranked molecules was used as a structure-aware sanity check, not as proof of binding affinity, therapeutic efficacy, or prospective discovery.

## Selection Criteria

The top-5 molecules were selected from existing ranked EGFR records using high triage score, applicability-domain support, low/acceptable uncertainty, no PAINS/Brenk/unwanted-substructure alert, acceptable drug-likeness, and scaffold diversity where possible.
Selection table: `reports/egfr_top5_structure_selection.csv`

## Reference Redocking Context

- PDB ID: 5UG9
- Reference ligand: 8AM
- Reference redocking score: -9.471 kcal/mol
- Reference redocking RMSD: 0.968 A

## Docking Setup

- Receptor PDBQT: `data/structure_prepared/5UG9_receptor.pdbqt`
- Docking box center: {'x': -13.8704375, 'y': 14.9338125, 'z': -26.729937500000002}
- Docking box size: {'x': 21.697, 'y': 22.433999999999997, 'z': 22.741}
- Backend preference: Python Vina API, then Vina CLI fallback

## Top-5 Results

| molecule_id | rank_before_docking | predicted_pIC50 | docking_status | vina_score_kcal_mol | shared_contact_fraction_with_8AM | distance_to_8AM_centroid | structure_sanity_label |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CHEMBL5997498 | 1 | 10.029 | vina_python_completed | -8.991 | 0.000 | 1.956 | structure_sanity_warning |
| CHEMBL5790648 | 2 | 9.872 | vina_python_completed | -8.386 | 0.000 | 2.150 | structure_sanity_warning |
| CHEMBL174426 | 5 | 9.710 | vina_python_completed | -8.878 | 0.000 | 2.252 | structure_sanity_warning |
| CHEMBL4749862 | 8 | 9.338 | vina_python_completed | -8.671 | 0.000 | 2.488 | structure_sanity_warning |
| CHEMBL2031299 | 9 | 9.131 | vina_python_completed | -8.665 | 0.000 | 0.741 | structure_sanity_warning |

## Contact Overlap With 8AM

- Mean shared contact fraction with 8AM: 0.0
- Reference 8AM contact count: 19
- Primary figure: `reports/figures/egfr_top5_docking_scores_and_contact_overlap.png`
- Overlay helper script: `reports/structure_visualization/top5_pose_overlay.pml`

## Structure-Aware Sanity Labels

- structure_sanity_warning: 5

## Limitations

- Docking scores are Vina scoring-function outputs and should not be interpreted as physical binding energies.
- Contact classes are heuristic residue-contact annotations.
- This is a retrospective sanity check over existing ranked molecules.
- Molecules that pass this check remain computationally prioritized existing records, not experimentally validated binding or inhibition claims.

## Reproducibility

```bash
python scripts/agentic_top5_structure_sanity.py
```
