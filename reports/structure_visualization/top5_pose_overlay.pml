reinitialize
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/5UG9_protein.pdb, receptor_5UG9
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/5UG9_8AM_ligand.pdb, ligand_8AM_reference
show cartoon, receptor_5UG9
show sticks, ligand_8AM_reference
color cyan, ligand_8AM_reference
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/top5_docked/CHEMBL5997498_vina_out.pdbqt, pose_CHEMBL5997498
show sticks, pose_CHEMBL5997498
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/top5_docked/CHEMBL5790648_vina_out.pdbqt, pose_CHEMBL5790648
show sticks, pose_CHEMBL5790648
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/top5_docked/CHEMBL174426_vina_out.pdbqt, pose_CHEMBL174426
show sticks, pose_CHEMBL174426
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/top5_docked/CHEMBL4749862_vina_out.pdbqt, pose_CHEMBL4749862
show sticks, pose_CHEMBL4749862
load /home/dimit/eva/Computational-Chemistry/pharma_project/egfr-cadd-qsar-admet/data/structure_prepared/top5_docked/CHEMBL2031299_vina_out.pdbqt, pose_CHEMBL2031299
show sticks, pose_CHEMBL2031299
zoom ligand_8AM_reference
set ray_opaque_background, off
