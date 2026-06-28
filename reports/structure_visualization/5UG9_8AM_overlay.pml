reinitialize
load data/structure_prepared/5UG9_8AM_ligand.pdbqt, reference_ligand
load data/structure_prepared/5UG9_8AM_redocked_out.pdbqt, redocked_pose
show sticks, reference_ligand
show sticks, redocked_pose
color cyan, reference_ligand
color magenta, redocked_pose
set ray_opaque_background, off
zoom
