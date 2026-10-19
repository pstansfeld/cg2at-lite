
cg2at_lite -c 6TIQ_CG_lastframe.gro -fg martini_3-2_nextgen_protein_charmm36 martini_3-1_new_lipidome_charmm36 -ff charmm36-cg2at -w tip3p -a minimized_protein_noh_AARef.pdb

cg2at_refine -cg2at CG2AT* -tpr pro.tpr -xtc pro_center.xtc -step 1 --ntmpi 1 --ntomp 4 --dphi 10

The CG trajectory (pro_center.xtc) is not included; supply your own CG trajectory for this system.
