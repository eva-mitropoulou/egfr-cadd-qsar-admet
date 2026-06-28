# EGFR Project Interview Talking Points

- Why scaffold split matters: random splits can place close analogs in both train and test, inflating apparent QSAR performance.
- Why applicability domain matters: predictions are more reliable when a molecule is close to training chemistry; the project quantified this with max Tanimoto similarity.
- Why no efficacy claim is made: all outputs are retrospective model benchmarks and existing-molecule triage, not prospective validation.
- Why Morgan fingerprints beat descriptors: fingerprints encode substructure patterns that are more relevant to kinase inhibitor SAR than global properties alone.
- What could come next: docking validation, experimental feedback loops, true ADMET predictors, and protein-ligand MD on selected validated poses.
