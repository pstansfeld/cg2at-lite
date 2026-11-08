#!/usr/bin/env python3

import os
import copy
import numpy as np
from cg2at_lite.bin import gen, g_var, at_mod
from cg2at_lite.bin.exceptions import CG2ATError



def build_atomistic_system(residue_type):
    system={}
#### for each residue type covert to atomistic except protein
    print('Converting residue type: ' +residue_type)
#### creates folder for residue type
    gen.mkdir_directory(g_var.working_dir+residue_type)
    if not os.path.exists(g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb'):
        atomistic_fragments, system[residue_type] = at_np_solvent(residue_type, g_var.cg_residues[residue_type])
        write_solvent(atomistic_fragments, residue_type)
    elif residue_type in g_var.sol_residues:
        atomistic_fragments, system[residue_type] = read_solvent_conversion(residue_type, g_var.cg_residues[residue_type])
    else:
        system[residue_type]=len(g_var.cg_residues[residue_type])
            
    print('{:<100}'.format(''), end='\r')
    print('Finished converting: '+residue_type)
    return system 

def write_solvent(atomistic_fragments, residue_type):
    merge_coord = [line['coord'] for line in atomistic_fragments]
    if residue_type in g_var.np_residues:
        merge_coord = fix_periodic_overlaps(merge_coord)
    coord, index_conversion = at_mod.index_conversion_generate(atomistic_fragments, merge_coord)
    at_mod.write_pdb(atomistic_fragments, coord, index_conversion, g_var.working_dir+residue_type+'/'+residue_type+'_all.pdb')

def fix_periodic_overlaps(coordinates):
#### jitters atoms closer than g_var.args.ov to another atom, looking across the periodic boundary (a
#### molecule that crosses it has atoms outside the box, which a plain distance check would miss; the
#### groups re-placed by correct_group_chirality are not fitted to beads and can land on a neighbour).
#### Orthorhombic boxes only; others are returned unchanged
    box = at_mod.pdb_box([g_var.box_vec]) if isinstance(getattr(g_var, 'box_vec', None), str) else None
    if box is None or len(coordinates) < 2 or np.abs(box - np.diag(np.diag(box))).max() > 1e-6:
        return coordinates
    from scipy.spatial import cKDTree
    lengths = np.diag(box)
    xyz = np.array(coordinates, dtype=float)
    for attempt in range(50):
        pairs = cKDTree(np.mod(xyz, lengths), boxsize=lengths).query_pairs(g_var.args.ov + 0.05)
        if not pairs:
            break
        for i in {min(pair) for pair in pairs}:
            xyz[i] += np.random.uniform(-0.2, 0.2, 3)
    return [x for x in xyz]

def read_solvent_conversion(cg_residue_type,cg_residues):
    residue_type={}
    residue_type_mass={}
    for cg_resid, cg_residue in enumerate(cg_residues):
        frag_location=gen.fragment_location(cg_residue_type) ### get fragment location from database
        residue_type[cg_residue_type], residue_type_mass[cg_residue_type] = at_mod.get_atomistic(frag_location, cg_residue_type)

        for group in residue_type[cg_residue_type]:
            sol_p_bead = 0
            for frag in residue_type[cg_residue_type][group].values():
                for atom in frag.values():
                    if int(atom['resid_ori']) > sol_p_bead:
                        sol_p_bead = int(atom['resid_ori'])
            return sol_p_bead, sol_p_bead*len(cg_residues)
    
    raise CG2ATError('There is an issue with the solvent recalculation')

def at_np_solvent(cg_residue_type,cg_residues):   
    atomistic_fragments={}  #### residue dictionary
#### run through every residue in a particular residue type
    atomistic_fragments_list = []
    # Fragment location and template are constant for all copies of the same
    # residue type — look them up once before the loop.
    frag_location = gen.fragment_location(cg_residue_type)
    frag_template, frag_mass_template = at_mod.get_atomistic(frag_location, cg_residue_type)
    template = at_mod.fragment_template(cg_residue_type, frag_template)

    total = len(cg_residues)
    for cg_resid, cg_residue in enumerate(cg_residues):
        gen.print_progress('Converting ' + cg_residue_type, cg_resid + 1, total)
        atomistic_fragments[cg_resid] = fit_np_residue(cg_residue_type, cg_residue, cg_residues[cg_residue],
                                                       frag_template, frag_mass_template, template, cg_resid)
        atomistic_fragments_list, sol_p_bead =  sort_np_dictionary(atomistic_fragments[cg_resid], atomistic_fragments_list)
    gen.finish_progress('Converting ' + cg_residue_type, total)
    if cg_residue_type in g_var.sol_residues:
        return atomistic_fragments_list, sol_p_bead*len(cg_residues)
    else:
        return atomistic_fragments_list, len(atomistic_fragments)

def fit_np_residue(cg_residue_type, cg_residue, cg_beads, frag_template, frag_mass_template, template, cg_resid=0):
#### fits one copy of a non-protein fragment onto its CG beads; returns {fragment atom number: atom}
    fitted = {cg_resid: {}}
    # Deep-copy the template so each residue gets an independent, mutable copy.
    residue = copy.deepcopy(frag_template)
    residue_mass = copy.deepcopy(frag_mass_template)
    for group in residue:
        center, at_frag_centers, cg_frag_centers, group_fit = at_mod.rigid_fit(residue[group], residue_mass,
                                                                                cg_residue, cg_beads)
        at_connect, cg_connect = at_mod.connectivity(cg_beads, at_frag_centers, cg_frag_centers, group_fit, group)
        xyz_rot_apply = at_mod.get_rotation(cg_connect, at_connect, center, cg_residue_type, group, cg_resid)
        fitted = at_mod.apply_rotations(fitted, cg_resid, group_fit, center, xyz_rot_apply)
        fitted[cg_resid] = at_mod.correct_group_spin(fitted[cg_resid], group_fit, template)
    if cg_residue_type in g_var.np_residues:
#### the mending below works on the molecule made whole (a molecule split across the periodic boundary
#### would have its groups compared across it); each atom then goes back to its own periodic image, as
#### the overlap check that follows does not look across the boundary
        box = at_mod.pdb_box([g_var.box_vec]) if isinstance(getattr(g_var, 'box_vec', None), str) else None
        shifts = at_mod.make_whole(fitted[cg_resid], cg_residue_type, frag_template, box)
#### stereocentres between groups take the fragment's handedness (CG beads can be in the mirror-image
#### arrangement around them)
        if len(residue) > 1:
            groups = {group: [atom for bead in residue[group].values() for atom in bead] for group in residue}
            fitted[cg_resid], _ = at_mod.correct_group_chirality(fitted[cg_resid], cg_residue_type, frag_template, groups)
        fitted[cg_resid] = at_mod.check_hydrogens(fitted[cg_resid])
        fitted[cg_resid] = at_mod.correct_planar_links(fitted[cg_resid], cg_residue_type, frag_template)
        fitted[cg_resid] = at_mod.correct_stereocentres(fitted[cg_resid], cg_residue_type, frag_template)
        fitted[cg_resid] = at_mod.place_planar_hydrogens(fitted[cg_resid], cg_residue_type, frag_template)
        for atom, shift in shifts.items():
            if atom in fitted[cg_resid]:
                fitted[cg_resid][atom]['coord'] = np.array(fitted[cg_resid][atom]['coord'], dtype=float) - shift
    if len(g_var.res_top[cg_residue_type]['CHIRAL']) > 1:
        chiral = at_mod.get_chiral_non_carbonyl(fitted[cg_resid])
        fitted[cg_resid] = at_mod.correct_chiral_atoms(fitted[cg_resid], chiral)
    return fitted[cg_resid]

def sort_np_dictionary(atomistic_fragments, atomistic_fragments_list):
    sol_p_bead = 0
    dict_start = np.min(list(atomistic_fragments.keys()))
    for atom_val in range(dict_start, len(atomistic_fragments)+dict_start):
        atom = atomistic_fragments[atom_val]
        atom['atom_name'] = atom.pop('atom')
        atom['residue_name'] = atom.pop('res_type')
        atom['residue_id'] = atom.pop('resid')
        atomistic_fragments_list.append(atom)    
        if atom['resid_ori'] > sol_p_bead:
            sol_p_bead = atom['resid_ori']
    return atomistic_fragments_list, sol_p_bead 
