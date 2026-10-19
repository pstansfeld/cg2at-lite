
cg2at_lite -c 4QWP_CG.gro -fg martini_3-0_charmm36 -ff charmm36-cg2at -w tip3p

Backward-compatibility check: 4QWP_CG.gro uses the old Martini 3 protein model (single BB bead per residue), in water with ions. It should convert with the martini_3-0_charmm36 library as before the explicit N/CA/C/O backbone (martini_3-2_nextgen_protein_charmm36) was added. No atomistic reference or trajectory is supplied, so this is a de novo conversion only.
