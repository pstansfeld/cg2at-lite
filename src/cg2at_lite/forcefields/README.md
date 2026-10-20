Force fields shipped with CG2AT-lite:

- charmm36-cg2at.ff = CHARMM36m protein / CHARMM36 lipids (July 2022 toppar), ported to GROMACS with charmm2gmx.
                      See charmm36-cg2at.ff/forcefield.doc for the source files and references.

Select it with -ff charmm36-cg2at. The earlier CHARMM36 releases (charmm36-mar2019,
charmm36-mar2019-updated, charmm36-jul2020, charmm36-jul2020-updated, charmm36-jul2021,
charmm36-jul2022) are accepted as aliases for it; the LJ-PME variants are not.

Adding a force field:

- Place a GROMACS force field folder (name.ff) in this directory; it is then listed by -info and can be chosen with -ff name.
- Alternatively, pass the path to a force field folder directly, e.g. -ff /path/to/name.ff.
- Prefix a folder name with _ to stop CG2AT2 from reading it.


Charmm36 references:

Vanommeslaeghe, K. Hatcher, E. Acharya, C. Kundu, S. Zhong, S. Shim, J. E. Darian, E. Guvench, O. Lopes, P. Vorobyov, I. and MacKerell, Jr. A.D. "CHARMM General Force Field (CGenFF): A force field for drug-like molecules compatible with the CHARMM all-atom additive biological force fields," Journal of Computational Chemistry 31: 671-90, 2010, PMC2888302 

Vanommeslaeghe, K., and MacKerell Jr., A.D., "Automation of the CHARMM General Force Field (CGenFF) I: bond perception and atom typing," Journal of Chemical Informationa and Modeling, 52: 3144-3154, 2012, PMC3528824 

Vanommeslaeghe, K., Raman, E.P., and MacKerell Jr., A.D., "Automation of the CHARMM General Force Field (CGenFF) II: Assignment of bonded parameters and partial atomic charges, Journal of Chemical Informationa and Modeling, 52: 3155-3168, 2012, PMC3528813 

Yu, W., He, X., Vanommeslaeghe, K. and MacKerell, A.D., Jr., "Extension of the CHARMM General Force Field to Sulfonyl-Containing Compounds and Its Utility in Biomolecular Simulations," Journal of Computational Chemistry, 33: 2451-2468, 2012, PMC3477297

Soteras Gutierrez, I., Lin, F.-Y., Vanommeslaeghe, K., Lemkul, J.A., Armacost, K.A., Brooks, Cl., III, and MacKerell, A.D., Jr., "Parametrization of Halogen Bonds in the CHARMM General Force Field: Improved treatment of ligand-protein interactions," Bioorganic & Medicinal Chemistry, In Press, 2016,
