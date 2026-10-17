
python cg2at.py -c 3D4S_CG_lastframe.gro -fg martini_3-2_nextgen_protein_charmm36 martini_3-1_new_lipidome_charmm36 -ff charmm36-jul2022 -w tip3p -a minimized_protein_noh_AARef.pdb

python cg2at_refine.py -cg2at CG2AT* -tpr pro.tpr -xtc pro_center.xtc -step 1 --ntmpi 1 --ntomp 4
