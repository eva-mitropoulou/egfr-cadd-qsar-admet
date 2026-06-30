# Top-5 EGFR Docking Score Sanity Check

Five existing EGFR molecules were selected from the final ranked and triaged table.
Selection prioritized high predicted activity, applicability-domain support, low uncertainty, medchem-alert cleanliness, and scaffold diversity.

Docked the top 5 clean, scaffold-diverse, high-ranked existing EGFR molecules into the validated 5UG9 binding-site setup and reported Vina scores as structure-aware triage annotations.

The 8AM reference ligand remains used only for the separate 5UG9/8AM redocking validation of the docking setup, not for top-5 scoring.

## Reference Redocking Context

- PDB ID: 5UG9
- Reference ligand: 8AM
- Reference redocking score: -9.471 kcal/mol
- Reference redocking RMSD: 0.968 A

## Top-5 Vina Score Results

| molecule_id | rank_before_docking | scaffold_id | predicted_pIC50 | conformal_interval_width | applicability_domain_bin | medchem_alert_flag | docking_status | vina_score_kcal_mol | docking_note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CHEMBL5997498 | 1 | ebdc73040f55 | 10.029 | 2.725 | high | False | vina_python_completed | -8.991 | Docking succeeded; Vina score retained as structure-aware triage annotation only. |
| CHEMBL5790648 | 2 | a7d5a3b39d24 | 9.872 | 2.725 | high | False | vina_python_completed | -8.386 | Docking succeeded; Vina score retained as structure-aware triage annotation only. |
| CHEMBL174426 | 5 | c28414362dcd | 9.710 | 2.725 | high | False | vina_python_completed | -8.878 | Docking succeeded; Vina score retained as structure-aware triage annotation only. |
| CHEMBL4749862 | 8 | 7a4199cf9767 | 9.338 | 2.725 | high | False | vina_python_completed | -8.671 | Docking succeeded; Vina score retained as structure-aware triage annotation only. |
| CHEMBL2031299 | 9 | 3cefafe781c8 | 9.131 | 2.725 | high | False | vina_python_completed | -8.665 | Docking succeeded; Vina score retained as structure-aware triage annotation only. |

- Successful ligand preparations: 5/5
- Successful dockings: 5/5
- Vina scores ranged from -8.991 to -8.386 kcal/mol.
- Figure: `reports/figures/egfr_top5_vina_scores.png`

## Limitations

- Vina scores are structure-aware triage annotations only.
- The top-5 docking stage does not validate binding, affinity, inhibition, biological activity, or discovery status.
- The top-5 docking stage does not use 8AM as a top-5 comparison target.

## Reproducibility

```bash
python scripts/agentic_top5_structure_sanity.py
```
