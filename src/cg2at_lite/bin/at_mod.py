#!/usr/bin/env python3

import os, sys
import copy
import numpy as np
import itertools
import math
from functools import lru_cache
from scipy.spatial import cKDTree
from cg2at_lite.bin import gen, g_var, at_mod_p, read_in, gro
from cg2at_lite.bin.exceptions import CG2ATError, InputError, FragmentNotFoundError, TopologyError


### sanity checking

def sanity_check_fragments(res, cg, sin_bead):
#### fetches bead and atom info from fragment database
    location = gen.fragment_location(res)
    residue, fragment_mass = get_atomistic(location, res)
    atom_list = []
    bead_list = []
    for group in residue.values():
        for bead in group:
            bead_list.append(bead)
            if not sin_bead:
                for atom in group[bead]:
                    atom_list.append(atom)
            else:
                if bead == sin_bead:
                    for atom in group[bead]:
                        atom_list.append(atom)
    return sorted(bead_list), atom_list

def sanity_check_atoms(atom_list, res):
#### checks atom order
    for at_num in range(1, len(atom_list)+1):
        if at_num not in atom_list:
            raise FragmentNotFoundError('atom number '+str(at_num)+' is missing from fragment library: '+res)

def sanity_check_beads(bead_list, cg, res):
#### checks if bead is in fragment library
    bead_list_new=[]
    for bead in cg:
        if bead not in bead_list:
            new_bead = bead[1:]+bead[0]
            if new_bead in bead_list and new_bead not in cg:
                pass
            else:
                raise FragmentNotFoundError('The bead '+bead+' is missing from the fragment library: '+res)   
        bead_list_new.append(bead)
    return bead_list_new

def sanity_check_protein_other(res_type, test=False):
    bead_list, atom_list = {},{}
    for residue in g_var.cg_residues[res_type]:
        for bead in g_var.cg_residues[res_type][residue]:
            resname = g_var.cg_residues[res_type][residue][bead]['residue_name']
            break
        if resname not in bead_list or resname not in atom_list:
            bead_list[resname], atom_list[resname] = sanity_check_fragments(resname, g_var.cg_residues[res_type][residue], False)
            sanity_check_atoms(atom_list[resname], resname)
        bead_list_cg = sanity_check_beads(bead_list[resname], g_var.cg_residues[res_type][residue], resname)  
        if bead_list[resname] != sorted(bead_list_cg):
            if len(bead_list[resname]) == len(bead_list_cg):
                fix_atom_wrap(bead_list[resname], bead_list_cg, res_type, residue)
            else:
                if not test:
                    print('There is a issue with residue: '+resname+' '+str(residue+1)+'. If expected ignore this message.' )
                if res_type != 'OTHER':
                    raise TopologyError('number of atomistic fragments: '+str(len(bead_list[resname]))+' does not equal number of CG beads: '+str(len(bead_list_cg)))

def sanity_check_solvent(res_type):
    bead_list, atom_list = {},{}
    for residue in g_var.cg_residues[res_type]:
        for bead in g_var.cg_residues[res_type][residue]:
            sin_bead=bead
            break
        if sin_bead not in bead_list or sin_bead not in atom_list:
            bead_list[sin_bead], atom_list[sin_bead] = sanity_check_fragments(res_type, g_var.cg_residues[res_type][residue], sin_bead)
            sanity_check_atoms(atom_list[sin_bead], res_type)
        sanity_check_beads(bead_list[sin_bead], g_var.cg_residues[res_type][residue], res_type)

def sanity_check_non_protein(res_type):
    bead_list, atom_list = {},{}
    if res_type not in bead_list or res_type not in atom_list:
        bead_list[res_type], atom_list[res_type] = sanity_check_fragments(res_type, g_var.cg_residues[res_type], False)
    sanity_check_atoms(atom_list[res_type], res_type)
    for residue in g_var.cg_residues[res_type]:
        bead_list_cg = sanity_check_beads(bead_list[res_type], g_var.cg_residues[res_type][residue], res_type)
        if bead_list[res_type] != sorted(bead_list_cg):
            if len(bead_list[res_type]) == len(bead_list_cg):
                fix_atom_wrap(bead_list[res_type], bead_list_cg, res_type, residue)
            else:
                print('There is a issue with residue: '+res_type+' '+str(residue+1))
                raise TopologyError('number of atomistic fragments: '+str(len(bead_list[res_type]))+' does not equal number of CG beads: '+str(len(bead_list_cg)))

def sanity_check():
#### runs through every bead and checks whether it exists
    for res_type in g_var.cg_residues:
        if res_type in ['PROTEIN', 'OTHER']:
            sanity_check_protein_other(res_type)
        elif res_type in g_var.sol_residues or res_type in g_var.ion_residues:
            sanity_check_solvent(res_type)
        else:
            sanity_check_non_protein(res_type)

def fix_atom_wrap(bead_list_frag, bead_list_cg, section, resid):
    for bead in bead_list_cg:
        if bead not in bead_list_frag:
            new_bead = bead[1:]+bead[0]
            if new_bead in bead_list_frag and new_bead not in bead_list_cg:
                g_var.cg_residues[section][resid][new_bead] = g_var.cg_residues[section][resid][bead]
                del g_var.cg_residues[section][resid][bead]
            else:
                print('There is a issue with residue: '+section+' '+str(resid+1))
                print('input file list:\n',sorted(bead_list_cg))
                print('\ncannot find: '+bead+' or '+new_bead+' in fragment list:')
                raise FragmentNotFoundError('Fragment bead mismatch. Available: '+str(sorted(bead_list_frag)))

#####  Sanity check end

#### rotations section

def rotate_atom(coord, center,xyz_rot_apply):
#### rotates atom around center
    coord =  coord-center  #### centers COM coordinates to 0,0,0
    if np.any(xyz_rot_apply): 
        coord=np.array(coord).dot(xyz_rot_apply)
    coord =  coord+center #### translates coord back by original offset
    return coord

def kabsch_rotate(at_connections,cg_connections):
    ## A solution for the best rotation to relate two sets of vectors. By WOLFGANG KABSCH
    ## http://en.wikipedia.org/wiki/Kabsch_algorithm
    Cov_mat = np.dot(np.transpose(at_connections), cg_connections)
    V, S, W = np.linalg.svd(Cov_mat)
    reflection = np.linalg.det(V) * np.linalg.det(W)
    if reflection < 0:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    rot_mat = np.dot(V, W)
    return rot_mat

## rotations end
## alignments 

def find_cross_vector(ca):
#### finds cross vector of the atoms CA, CA+1, CA+2 
    AB = ca[0]-ca[1]
    AC = ca[0]-ca[2] 
    cross_vector = np.cross(AB, AC)
    cross_vector = cross_vector/np.linalg.norm(cross_vector)
    return cross_vector

def noramlised_vector(c1, c2):
    initial_vector= c1-c2
    initial_vector = initial_vector/np.linalg.norm(initial_vector)
    return initial_vector


def align_to_vector(v1, v2):
#### returns the rotation matrix to rotate v1 to v2
    v = np.cross(v1,v2)
    c = np.dot(v1,v2)
    s = np.linalg.norm(v)

    rotation=np.array([[   0, -v[2],  v[1]],
                       [v[2],     0, -v[0]],
                       [-v[1], v[0],    0], 
                       ])

    r = np.identity(3) - rotation + np.matmul(rotation,rotation) * ((1 - c)/(s**2))
    return r

def align_at_frag_to_CG_frag(at_com, cg_com, group):
#### aligns atomistic fragment to cg bead
    COM_vector=at_com-cg_com ### gets vector between COM of atoms in fragment and cg bead 
    for bead in group:
        for atom in group[bead]: ### runs through atoms in fragments and centers on the cg bead 
            group[bead][atom]['coord']=group[bead][atom]['coord']-COM_vector
    return group, COM_vector

def COM(mass, fragment):
#### returns center of mass of fragment
    try:
        if np.any(np.array(mass)[:,3]):   
            return np.average(np.array(mass)[:,:3], axis=0, weights=np.array(mass)[:,3])
    except:
        for bead in fragment:
            print(bead, fragment[bead], '\n')
            print(mass)
            raise TopologyError('missing the mass one of the atoms in '+fragment[bead][1]['res_type'])      

def rigid_fit(group, frag_mass, resid, cg):
#### rigid fits group to CG beads
    rigid_mass_at = []
    rigid_mass_cg = []
    for bead in group:
        rigid_mass_at+=frag_mass[bead]
        rigid_mass_cg.append(cg[bead]['coord'])
    rigid_mass_at = COM(rigid_mass_at, group)
    rigid_mass_cg = np.mean(rigid_mass_cg, axis=0)
    
    group, COM_vector = align_at_frag_to_CG_frag(rigid_mass_at, rigid_mass_cg, group)
    at_frag_centers = {}
    cg_frag_centers = {}
    for bead in group:
        cg_frag_centers[bead] = cg[bead]['coord']
        at_frag_centers[bead] = COM(frag_mass[bead], bead)
        at_frag_centers[bead] -= COM_vector
    rigid_mass_cg, at_frag_centers, cg_frag_centers, group
    return rigid_mass_cg, at_frag_centers, cg_frag_centers, group

## alignment end

## overlap checker

def overlapping_atoms(tree):
    overlapped_ndx = tree.query_ball_tree(tree, r=g_var.args.ov)
    overlapped_cut = [ndx for ndx in overlapped_ndx if len(ndx) >1]      
    overlapped_cut.sort()
    overlapped=list(overlapped_cut for overlapped_cut,_ in itertools.groupby(overlapped_cut))
    return overlapped

def check_atom_overlap(coordinates, fixed=None):
    """Jitter any overlapping atoms until all pairwise distances exceed g_var.args.ov.

    fixed: coordinates of atoms that are never moved but must not be overlapped either (the
    massive dummies of rigid CH3/NH3 groups in virtual-site topologies); an atom on top of one
    is moved off it. The cKDTree is rebuilt once per pass over all overlaps rather than once
    per atom move, which gives a significant speedup for large systems.
    """
    n = len(coordinates)
    fixed = [] if fixed is None else list(fixed)
    def build():
        return cKDTree(np.array(list(coordinates) + fixed)) if fixed else cKDTree(coordinates)
    def movable(groups):
        out = []
        for group in groups:
            mobile = [i for i in group if i < n]
            if mobile:
                out.append(mobile + [i for i in group if i >= n])
        return out
    tree = build()
    overlapped = movable(overlapping_atoms(tree))
    while len(overlapped) > 0:
        if len(overlapped) > 30:
            gen.print_progress('Fixing atom overlaps', 0, len(overlapped))
        for ndx_val, ndx in enumerate(overlapped):
            if len(overlapped) > 30:
                gen.print_progress('Fixing atom overlaps', ndx_val + 1, len(overlapped))
            # Jitter the atom using the current tree for the inner acceptance check.
            # The tree is intentionally not rebuilt here; we rebuild once per pass below.
            xyz_check = coordinates[ndx[0]] + np.random.uniform(-0.2, 0.2, 3)
            while len(tree.query_ball_point(xyz_check, r=g_var.args.ov)) > 1:
                xyz_check = coordinates[ndx[0]] + np.random.uniform(-0.2, 0.2, 3)
            coordinates[ndx[0]] = xyz_check
        # Rebuild tree once per pass now that all atoms in this batch have moved.
        tree = build()
        if len(overlapped) > 30:
            gen.finish_progress('Fixing atom overlaps', len(overlapped))
        overlapped = movable(overlapping_atoms(tree))
    return coordinates

## overlap end
## get fragment information

def split_fragment_names(line, residue, resname):
    bead = gen.strip_header(line)
    group = g_var.res_top[resname]['GROUPS'][bead]
    if group not in residue:
        residue[group] = {}
    if bead not in residue[group]:
        residue[group][bead]={} 
    return residue, group, bead   

@lru_cache(maxsize=None)
def _get_atomistic_cached(frag_location: str, resname: str):
    """Parse a fragment PDB once and cache the result.

    Keyed on (frag_location, resname) so each unique residue type is read from
    disk exactly once, regardless of how many copies exist in the system.
    """
    residue = {}
    fragment_mass = {}
    with open(frag_location, 'r') as pdb_input:
        for line in pdb_input.readlines():
            if line.startswith('['):
                residue, group, bead = split_fragment_names(line, residue, resname)
                fragment_mass[bead] = []
            if line.startswith('ATOM'):
                line_sep = gen.pdbatom(line)
                residue[group][bead][line_sep['atom_number']] = {
                    'coord': np.array([line_sep['x'] * g_var.sf,
                                       line_sep['y'] * g_var.sf,
                                       line_sep['z'] * g_var.sf]),
                    'atom': line_sep['atom_name'], 'resid': 1,
                    'resid_ori': line_sep['residue_id'],
                    'res_type': line_sep['residue_name'],
                    'frag_mass': 1,
                }
                if not gen.is_hydrogen(line_sep['atom_name']):
                    if line_sep['atom_name'] in g_var.res_top[resname]['atom_masses']:
                        mass = g_var.res_top[resname]['atom_masses'][line_sep['atom_name']]
                        residue[group][bead][line_sep['atom_number']]['frag_mass'] = mass
                        fragment_mass[bead].append([line_sep['x'] * g_var.sf,
                                                    line_sep['y'] * g_var.sf,
                                                    line_sep['z'] * g_var.sf, mass])
                else:
                    fragment_mass[bead].append([line_sep['x'] * g_var.sf,
                                                line_sep['y'] * g_var.sf,
                                                line_sep['z'] * g_var.sf, 1])
    return residue, fragment_mass


def get_atomistic(frag_location: str, resname: str = False):
    """Return a per-call deep copy of the cached fragment data.

    Callers mutate the returned dicts (e.g. deleting unused beads in
    at_mod_p.py), so each call gets its own independent copy while the
    underlying parse result is only computed once per residue type.
    """
    if not resname:
        resname = frag_location.split('/')[-1][:-4]
    residue, fragment_mass = _get_atomistic_cached(frag_location, resname)
    return copy.deepcopy(residue), copy.deepcopy(fragment_mass)

def connectivity(cg, at_frag_centers, cg_frag_centers, group, group_number):
#### returns the connections between the atomistic and coarse grain beads
    at_connection, cg_connection=[],[]
    for bead in cg:
        resname = cg[bead]['residue_name']
        break
    if resname in g_var.sorted_connect:
        if len(g_var.sorted_connect[resname]) > 0:
            for group_bead in group: 
                for bead_atom in group[group_bead]:
                    if bead_atom in g_var.sorted_connect[resname][int(group_number)]:
                        skip = False
                        for bead_connect in g_var.sorted_connect[resname][int(group_number)][bead_atom]:
                            if 'cg_connect' not in locals():
                                if bead_connect in cg:
                                    cg_connect = [cg[bead_connect]['coord']]
                                else:
                                    skip = True
                            else:
                                if bead_connect in cg:
                                    cg_connect.append(cg[bead_connect]['coord'])
                        if not skip:
                            at_connection.append(group[group_bead][bead_atom]['coord'])
                            cg_connection.append(np.mean(np.array(cg_connect), axis=0))    
    if len(at_frag_centers)  > 1:         
        for bead in at_frag_centers:
            cg_connection.append(cg_frag_centers[bead])
            at_connection.append(at_frag_centers[bead])   
    return at_connection, cg_connection    

def get_rotation(cg_connect, at_connect, center, resname, group, cg_resid):
    if len(at_connect) == len(cg_connect):
        if len(cg_connect) == 0:
            return gen.AnglesToRotMat([np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2), np.random.uniform(0, math.pi*2)])
        else:
            return kabsch_rotate(np.array(at_connect)-center, np.array(cg_connect)-center)
    else:
        print('atom connections: '+str(len(at_connect))+' does not match CG connections: '+str(len(cg_connect)))
        raise CG2ATError('residue number: '+str(cg_resid)+', residue type: '+str(resname)+', group: '+group)

def apply_rotations(atomistic_fragments,cg_resid, group_fit, center, xyz_rot_apply):
    for bead in group_fit:
        for atom in group_fit[bead]:
            group_fit[bead][atom]['coord'] = rotate_atom(group_fit[bead][atom]['coord'], center, xyz_rot_apply)   
            atomistic_fragments[cg_resid][atom] = group_fit[bead][atom].copy()
    return atomistic_fragments

## Spin correction for fragment groups.
## A group joined to the rest of the residue by single bonds (e.g. a lipid tail bead)
## is only fixed by its connection points, so the rotational fit leaves its spin about
## the axis through them undetermined. Left arbitrary, neighbouring beads are twisted
## against each other and double bonds near bead boundaries come out trans. The spin
## is chosen instead to reproduce the fragment's own dihedrals across the join.

_template_cache = {}

def fragment_template(resname, frag_template):
#### heavy-atom coordinates and bonds of the whole-molecule fragment (cached per residue)
    if resname not in _template_cache:
        coords = {}
        for group in frag_template.values():
            for bead in group.values():
                for atom_number, atom in bead.items():
                    if not gen.is_hydrogen(atom['atom']):
                        coords[atom_number] = np.array(atom['coord'], dtype=float)
        numbers = list(coords)
        bonded = {number: set() for number in numbers}
        if len(numbers) > 1:
            xyz = np.array([coords[number] for number in numbers])
            for i, j in cKDTree(xyz).query_pairs(1.75 * g_var.sf):
                bonded[numbers[i]].add(numbers[j])
                bonded[numbers[j]].add(numbers[i])
        _template_cache[resname] = (coords, bonded)
    return _template_cache[resname]

def six_membered_rings(bonded):
#### six-membered rings of the template's heavy-atom bond graph, as ordered atom tuples
    rings, seen = [], set()
    def extend(path):
        if len(path) == 6:
            if path[0] in bonded[path[-1]] and frozenset(path) not in seen:
                seen.add(frozenset(path))
                rings.append(tuple(path))
            return
        for nxt in bonded[path[-1]]:
            if nxt > path[0] and nxt not in path:
                extend(path + [nxt])
    for start in bonded:
        extend([start])
#### drop cycles with a chord (fused-ring outlines), keeping only true rings
    return [ring for ring in rings if all(len(bonded[atom] & set(ring)) == 2 for atom in ring)]

def ring_pucker_theta(coords):
#### Cremer-Pople theta (degrees) of a six-membered ring: ~0 or ~180 chair, ~90 boat/twist-boat
    r = np.array(coords) - np.mean(coords, axis=0)
    j = np.arange(6)
    normal = np.cross((r * np.sin(2*np.pi*j/6)[:, None]).sum(0), (r * np.cos(2*np.pi*j/6)[:, None]).sum(0))
    z = r.dot(normal / np.linalg.norm(normal))
    q2 = np.hypot(np.sqrt(2/6) * (z * np.cos(4*np.pi*j/6)).sum(), np.sqrt(2/6) * (z * np.sin(4*np.pi*j/6)).sum())
    q3 = np.sqrt(1/6) * (z * (-1)**j).sum()
    return math.degrees(math.atan2(q2, q3)), math.hypot(q2, q3)

def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w))

_stereo_cache = {}

def stereo_template(resname, frag_template):
#### sp3 stereocentres of the fragment (a carbon with three or four heavy neighbours and at most one
#### hydrogen): (centre, three reference neighbours, fourth substituent, its hydrogens, handedness).
#### The fourth substituent is the hydrogen, or for a quaternary carbon its least-connected heavy
#### neighbour (e.g. a methyl), whose own hydrogens move with it.
    if resname not in _stereo_cache:
        coords, bonded = fragment_template(resname, frag_template)
        hydrogens = {}
        for group in frag_template.values():
            for bead in group.values():
                for number, atom in bead.items():
                    if gen.is_hydrogen(atom['atom']):
                        hydrogens[number] = np.array(atom['coord'], dtype=float)
        h_numbers = list(hydrogens)
        h_tree = cKDTree(np.array([hydrogens[n] for n in h_numbers])) if h_numbers else None
        def attached(number):
            if h_tree is None:
                return []
            return [h_numbers[i] for i in h_tree.query_ball_point(coords[number], 1.25 * g_var.sf)]
        centres = []
        for c in coords:
            if not frag_atom_name(frag_template, c).startswith('C') or len(bonded[c]) < 3:
                continue
            hs = attached(c)
            if len(bonded[c]) + len(hs) != 4 or len(hs) > 1:
                continue
            if hs:
                fourth, carried, three = hs[0], [], sorted(bonded[c])
            else:
                fourth = min(bonded[c], key=lambda x: (len(bonded[x]), x))
                if len(bonded[fourth]) > 1:
                    continue
                carried, three = attached(fourth), sorted(bonded[c] - {fourth})
            v = [coords[x] - coords[c] for x in three]
            volume = np.dot(v[0], np.cross(v[1], v[2])) / np.prod([np.linalg.norm(x) for x in v])
            if abs(volume) > 0.1:
                centres.append((c, three, fourth, carried, np.sign(volume)))
        _stereo_cache[resname] = centres
    return _stereo_cache[resname]

def frag_atom_name(frag_template, number):
    for group in frag_template.values():
        for bead in group.values():
            if number in bead:
                return bead[number]['atom']
    return ''

def correct_stereocentres(residue, resname, frag_template):
#### a placed residue whose stereocentre has the opposite handedness to the fragment (its substituents
#### were fitted in different beads) is mended locally: the centre is reflected through the plane of
#### its three reference neighbours, which keeps those three bonds, and the fourth substituent (with
#### any hydrogens it carries) is put back opposite them at its bond length. CISRES then holds it.
    for c, three, fourth, carried, sign in stereo_template(resname, frag_template):
        if any(x not in residue for x in [c, fourth] + three + carried):
            continue
        X = {x: np.array(residue[x]['coord'], dtype=float) for x in [c, fourth] + three + carried}
        v = [X[x] - X[c] for x in three]
        if np.sign(np.dot(v[0], np.cross(v[1], v[2]))) == sign:
            continue
#### a ring atom is mended by swapping its two substituents outside the ring instead: reflecting the
#### centre itself flips the ring's pucker there, and the ring restraints then stall the minimisation
#### (GlcN I C3 of lipid A, dragged through by its ester in the soft pass)
        mended = swap_ring_substituents(residue, resname, frag_template, c)
        if mended is not None:
            continue
        normal = np.cross(X[three[1]] - X[three[0]], X[three[2]] - X[three[0]])
        normal /= np.linalg.norm(normal)
        new_c = X[c] - 2 * np.dot(X[c] - X[three[0]], normal) * normal
        outward = new_c - np.mean([X[x] for x in three], axis=0)
        new_fourth = new_c + np.linalg.norm(X[fourth] - X[c]) * outward / np.linalg.norm(outward)
        shift = new_fourth - X[fourth]
        residue[c]['coord'] = new_c
        residue[fourth]['coord'] = new_fourth
        for h in carried:
            residue[h]['coord'] = X[h] + shift
    return residue

_ring_atom_cache = {}

def swap_ring_substituents(residue, resname, frag_template, c):
#### inverts a ring stereocentre without moving the ring: its two substituents outside the ring (a
#### heavy atom and a hydrogen) are put in the two tetrahedral positions outside the ring, the way round that
#### gives the fragment's handedness, at their fragment bond lengths. A substituent's hydrogens move
#### with it, and so does a small group hanging off it (up to 4 heavy atoms, e.g. a carboxylate or a
#### hydroxymethyl); a larger one (an acyl chain, the rest of the molecule) stays, leaving one stretched
#### bond for the minimisation. Returns None when c is not in a six-membered ring of the fragment (it is
#### then mended by reflection), False for a ring centre left as it is, True when mended
    coords, bonded = fragment_template(resname, frag_template)
    if resname not in _ring_atom_cache:
        _ring_atom_cache[resname] = {atom for ring in six_membered_rings(bonded) for atom in ring}
    ring_atoms = _ring_atom_cache[resname]
    if c not in ring_atoms:
        return None
    parents = hydrogen_parents(resname, frag_template)
    in_ring = [n for n in bonded[c] if n in ring_atoms]
    outside = [n for n in bonded[c] if n not in ring_atoms] + [h for h, p in parents.items() if p == c]
#### only a centre with a hydrogen outside the ring: with two heavy substituents (e.g. the anomeric
#### carbon of Kdo, between its carboxylate and the rest of lipid A) one of them would have to move a
#### whole branch's length, so it is left to the minimisation, as is any other ring centre
    if (len(in_ring) != 2 or len(outside) != 2 or any(x not in residue for x in [c] + outside)
            or not any(gen.is_hydrogen(frag_atom_name(frag_template, x)) for x in outside)):
        return False
    def group(root):
        if gen.is_hydrogen(frag_atom_name(frag_template, root)):
            return [root]
        heavy, stack = {root}, [root]
        while stack:
            for n in bonded[stack.pop()]:
                if n != c and n not in heavy:
                    heavy.add(n)
                    stack.append(n)
        if len(heavy) > 4:
            heavy = {root}
        return list(heavy) + [h for h, p in parents.items() if p in heavy]
    unit = lambda d: d / np.linalg.norm(d)
    P = {x: np.array(residue[x]['coord'], dtype=float) for x in [c] + in_ring + outside}
    C = P[c]
#### the two tetrahedral positions outside the ring: about the bisector pointing away from the ring,
#### either side of the ring atoms' plane (a centre dragged flat has both substituents on one side, so
#### swapping them would not mend it)
    away = -unit(unit(P[in_ring[0]] - C) + unit(P[in_ring[1]] - C))
    normal = unit(np.cross(P[in_ring[0]] - C, P[in_ring[1]] - C))
    half = math.radians(54.75)
    slots = [unit(away * math.cos(half) + normal * math.sin(half)), unit(away * math.cos(half) - normal * math.sin(half))]
    sign = next(s for centre, three, fourth, carried, s in stereo_template(resname, frag_template) if centre == c)
    three = next(t for centre, t, fourth, carried, s in stereo_template(resname, frag_template) if centre == c)
#### bond lengths from the fragment: a centre dragged through can leave a bond squashed or stretched
    frag_xyz = {number: np.array(atom['coord'], dtype=float) for group in frag_template.values()
                for bead in group.values() for number, atom in bead.items()}
    length = {x: np.linalg.norm(frag_xyz[x] - frag_xyz[c]) / g_var.sf for x in outside}
    for order in (slots, slots[::-1]):
        new = {x: C + length[x] * d for x, d in zip(outside, order)}
        v = [(new[x] if x in new else P[x]) - C for x in three]
        if np.sign(np.dot(v[0], np.cross(v[1], v[2]))) == sign:
            break
    shifts = {x: new[x] - P[x] for x in outside}
    for root, shift in shifts.items():
        for atom in group(root):
            if atom in residue:
                residue[atom]['coord'] = np.array(residue[atom]['coord'], dtype=float) + shift
    return True

_planar_h_cache = {}

def planar_hydrogens(resname, frag_template):
#### hydrogens on planar centres of the fragment (an amide N-H, C=C-H or aromatic C-H: a parent with
#### two heavy neighbours, this one hydrogen, and all three substituents in the parent's plane):
#### (hydrogen, parent, the two heavy neighbours, bond length)
    if resname not in _planar_h_cache:
        coords, bonded = fragment_template(resname, frag_template)
        hydrogens = {}
        for group in frag_template.values():
            for bead in group.values():
                for number, atom in bead.items():
                    if gen.is_hydrogen(atom['atom']):
                        hydrogens[number] = np.array(atom['coord'], dtype=float)
        h_numbers = list(hydrogens)
        h_tree = cKDTree(np.array([hydrogens[n] for n in h_numbers])) if h_numbers else None
        out = []
        for p in coords:
            if len(bonded[p]) != 2 or h_tree is None:
                continue
            hs = [h_numbers[i] for i in h_tree.query_ball_point(coords[p], 1.25 * g_var.sf)]
            if len(hs) != 1:
                continue
            a, b = sorted(bonded[p])
            v = [coords[a] - coords[p], coords[b] - coords[p], hydrogens[hs[0]] - coords[p]]
            if abs(np.dot(v[0], np.cross(v[1], v[2]))) / np.prod([np.linalg.norm(x) for x in v]) < 0.1:
                out.append((hs[0], p, a, b, np.linalg.norm(v[2])))
        _planar_h_cache[resname] = out
    return _planar_h_cache[resname]

def place_planar_hydrogens(residue, resname, frag_template, scale=1.0):
#### such a hydrogen is rebuilt on the outward bisector of its parent's two heavy bonds: fitted in a
#### different bead from a neighbour (e.g. the previous residue's carbonyl C of an amide) it can land
#### across that bond, where the minimisation cannot pull it out
    for h, p, a, b, length in planar_hydrogens(resname, frag_template):
        if any(x not in residue for x in (h, p, a, b)):
            continue
        P, A, B = [np.array(residue[x]['coord'], dtype=float) for x in (p, a, b)]
        u = (P - A) / np.linalg.norm(P - A) + (P - B) / np.linalg.norm(P - B)
        if np.linalg.norm(u) > 1e-6:
            residue[h]['coord'] = P + scale * length * u / np.linalg.norm(u)
    return residue

def correct_minimised_residue(residue_type, pdb_in, pdb_out):
#### applies correct_stereocentres, correct_planar_links and place_planar_hydrogens to every molecule
#### of a per-residue structure. Stereocentres come first: mending a ring centre moves a substituent
#### (e.g. an ester oxygen), and the link through it is then turned back to its fragment dihedral
    def mend(residue, frag_template):
        residue = correct_stereocentres(residue, residue_type, frag_template)
        residue = correct_planar_links(residue, residue_type, frag_template)
        return place_planar_hydrogens(residue, residue_type, frag_template, 1.0 / g_var.sf)
    apply_per_molecule(residue_type, pdb_in, pdb_out, mend)

def distorted(residue, resname, frag_template):
#### a molecule the minimisation could not mend: a heavy bond stretched more than 10% (a ring held open
#### by its own restraints), or a stereocentre inverted or splayed nearly flat (two of its heavy bonds
#### close to linear, where dihedral restraints through it are singular)
    coords, bonded = fragment_template(resname, frag_template)
    for a in coords:
        for b in bonded[a]:
            if b > a and a in residue and b in residue:
                length = np.linalg.norm(np.array(residue[a]['coord'], dtype=float) - np.array(residue[b]['coord'], dtype=float))
                if length > 1.1 * np.linalg.norm(coords[a] - coords[b]) / g_var.sf:
                    return True
    for c, three, fourth, carried, sign in stereo_template(resname, frag_template):
        if any(x not in residue for x in [c] + three):
            continue
        v = [np.array(residue[x]['coord'], dtype=float) - np.array(residue[c]['coord'], dtype=float) for x in three]
        if np.sign(np.dot(v[0], np.cross(v[1], v[2]))) != sign:
            return True
        for a, b in itertools.combinations(v, 2):
            if np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)) < math.cos(math.radians(150)):
                return True
    return False

def replace_distorted_molecules(residue_type, pdb_in, pdb_out):
#### every distorted molecule of a per-residue structure is replaced by the whole fragment, rigidly fitted
#### to its heavy atoms: it loses the detail of its fit to the CG beads but keeps the fragment's rings,
#### stereocentres and links, and minimises cleanly from there. Returns the number replaced
    replaced = []
    def replace(residue, frag_template):
        if not distorted(residue, residue_type, frag_template):
            return residue
        whole = {number: np.array(atom['coord'], dtype=float) / g_var.sf
                 for group in frag_template.values() for bead in group.values() for number, atom in bead.items()}
        heavy = [n for n, atom in ((n, a) for group in frag_template.values() for bead in group.values() for n, a in bead.items())
                 if not gen.is_hydrogen(atom['atom']) and n in residue]
        P = np.array([whole[n] for n in heavy]); Q = np.array([residue[n]['coord'] for n in heavy], dtype=float)
        pc, qc = P.mean(axis=0), Q.mean(axis=0)
        rot = kabsch_rotate(P - pc, Q - qc)
        for number in residue:
            if number in whole:
                residue[number]['coord'] = (whole[number] - pc).dot(rot) + qc
        replaced.append(1)
        return residue
    apply_per_molecule(residue_type, pdb_in, pdb_out, replace)
    return len(replaced)

def pdb_box(lines):
#### box vectors (rows, Angstrom) from a pdb's CRYST1 record, or None without one
    for line in lines:
        if line.startswith('CRYST1'):
            try:
                a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
                alpha, beta, gamma = [math.radians(float(line[k:k + 7])) for k in (33, 40, 47)]
            except ValueError:
                return None
            if min(a, b, c) <= 1:
                return None
            cx = c * math.cos(beta)
            cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
            return np.array([[a, 0, 0], [b * math.cos(gamma), b * math.sin(gamma), 0],
                             [cx, cy, math.sqrt(max(c**2 - cx**2 - cy**2, 0))]])
        if line.startswith('ATOM'):
            return None
    return None

def minimum_image(d, box):
#### the shortest periodic image of a separation vector
    if box is None:
        return d
    frac = np.linalg.solve(box.T, d)
    return box.T.dot(frac - np.round(frac))

_h_parent_cache = {}

def hydrogen_parents(resname, frag_template):
#### each hydrogen of the fragment -> its nearest heavy atom (cached per residue)
    if resname not in _h_parent_cache:
        coords, _ = fragment_template(resname, frag_template)
        heavy = list(coords)
        parents = {}
        if heavy:
            tree = cKDTree(np.array([coords[n] for n in heavy]))
            for group in frag_template.values():
                for bead in group.values():
                    for number, atom in bead.items():
                        if gen.is_hydrogen(atom['atom']):
                            parents[number] = heavy[tree.query(np.array(atom['coord'], dtype=float))[1]]
        _h_parent_cache[resname] = parents
    return _h_parent_cache[resname]

def make_whole(residue, resname, frag_template, box):
#### moves the atoms of a molecule into one periodic image, along the fragment's heavy-atom bonds
#### (each hydrogen next to its nearest heavy atom in the fragment); returns the shift given to each
#### atom moved. Leaves the molecule as it is without a box
    if box is None:
        return {}
    coords, bonded = fragment_template(resname, frag_template)
    shifts, seen = {}, set()
    def move(number, reference):
        X = np.array(residue[number]['coord'], dtype=float)
        R = np.array(residue[reference]['coord'], dtype=float)
        new = R + minimum_image(X - R, box)
        if np.linalg.norm(new - X) > 1e-6:
            residue[number]['coord'] = new
            shifts[number] = new - X
    for root in coords:
        if root in seen or root not in residue:
            continue
        seen.add(root)
        queue = [root]
        while queue:
            current = queue.pop()
            for nxt in bonded[current]:
                if nxt not in seen and nxt in residue:
                    seen.add(nxt)
                    move(nxt, current)
                    queue.append(nxt)
    for h, parent in hydrogen_parents(resname, frag_template).items():
        if h in residue and parent in residue:
            move(h, parent)
    return shifts

_cis_cache = {}

def cis_double_bonds(resname, frag_template):
#### C=C bonds of the fragment (between two carbons with three neighbours each, hydrogens included,
#### outside rings) whose heavy substituents are cis: (i, a, b, l) for each cis heavy-atom dihedral
    if resname not in _cis_cache:
        coords, bonded = fragment_template(resname, frag_template)
        n_h = {number: 0 for number in coords}
        for h, parent in hydrogen_parents(resname, frag_template).items():
            n_h[parent] += 1
        names = {number: frag_atom_name(frag_template, number) for number in coords}
        ring_atoms = {atom for ring in six_membered_rings(bonded) for atom in ring}
        sp2 = {n for n in coords if names[n].startswith('C') and len(bonded[n]) + n_h[n] == 3 and n not in ring_atoms}
        out = []
        for a in sp2:
            for b in bonded[a]:
                if b <= a or b not in sp2:
                    continue
                for i in bonded[a] - {b}:
                    for l in bonded[b] - {a}:
                        if abs(math.degrees(dihedral(*[coords[x] for x in (i, a, b, l)]))) < 30:
                            out.append((i, a, b, l))
        _cis_cache[resname] = out
    return _cis_cache[resname]

def trans_double_bonds(residue_type, pdb_in):
#### counts the fragment's cis C=C bonds that came out trans (|dihedral| > 90 deg) in a per-residue
#### structure: (bonds turned, bonds checked, molecules with one or more)
    counts = [0, 0, 0]
    def check(residue, frag_template):
        turned = 0
        for i, a, b, l in cis_double_bonds(residue_type, frag_template):
            if all(x in residue for x in (i, a, b, l)):
                counts[1] += 1
                if abs(math.degrees(dihedral(*[np.array(residue[x]['coord'], dtype=float) for x in (i, a, b, l)]))) > 90:
                    turned += 1
        counts[0] += turned
        counts[2] += turned > 0
        return residue
    apply_per_molecule(residue_type, pdb_in, os.devnull, check)
    return tuple(counts)

def apply_per_molecule(residue_type, pdb_in, pdb_out, fix):
#### applies fix(residue, frag_template) to every molecule of a per-residue structure (molecules in itp
#### atom order, one after another); fragment atom -> itp position is by name when the itp names are
#### unique, otherwise the fragment atom number is the itp index
    location = gen.fragment_location(residue_type)
    frag_template, _ = get_atomistic(location, residue_type)
    names, in_atoms, molecules = {}, False, 0
    with open(location[:-4]+'.itp') as itp_input:
        for line in itp_input:
            line_sep = line.split(';')[0].split()
            if not line_sep:
                continue
            if line_sep[0] == '[':
                in_atoms = line_sep[1] == 'atoms'
                molecules += line_sep[1] == 'moleculetype'
                continue
            if in_atoms and molecules == 1 and line_sep[0].isdigit():
                names.setdefault(line_sep[4], []).append(int(line_sep[0]))
    n_atoms = sum(len(index) for index in names.values())
    unique = all(len(index) == 1 for index in names.values())
    position = {}
    for group in frag_template.values():
        for bead in group.values():
            for number, atom in bead.items():
                position[number] = (names[atom['atom']][0] if atom['atom'] in names else None) if unique else number
    with open(pdb_in) as pdb_input:
        lines = pdb_input.readlines()
    box = pdb_box(lines)
    atom_lines = [k for k, line in enumerate(lines) if line.startswith('ATOM')]
    for start in range(0, len(atom_lines) - n_atoms + 1, n_atoms):
        rows = atom_lines[start:start + n_atoms]
        residue = {number: {'coord': np.array([float(lines[rows[p - 1]][30:38]), float(lines[rows[p - 1]][38:46]),
                                                float(lines[rows[p - 1]][46:54])])}
                   for number, p in position.items() if p is not None and p <= n_atoms}
#### a molecule split across the periodic boundary is made whole first (its bonds would otherwise
#### measure a box length, and every one crossing the boundary was taken for a distortion); each atom
#### is put back in its own periodic image afterwards, so the structure is wrapped as before
        shifts = make_whole(residue, residue_type, frag_template, box)
        residue = fix(residue, frag_template)
        for number, shift in shifts.items():
            if number in residue:
                residue[number]['coord'] = np.array(residue[number]['coord'], dtype=float) - shift
        for number, p in position.items():
            if p is not None and number in residue:
                k = rows[p - 1]
                lines[k] = lines[k][:30] + '{:8.3f}{:8.3f}{:8.3f}'.format(*residue[number]['coord']) + lines[k][54:]
    with open(pdb_out, 'w') as pdb_output:
        pdb_output.writelines(lines)

_planar_link_cache = {}

def planar_links(resname, frag_template):
#### amide and ester links of the fragment (a carbonyl carbon bonded to an N or O that has a second heavy
#### neighbour): (carbonyl C, N/O, its other heavy neighbour, the carbonyl's other two heavy neighbours
#### with their hydrogens, fragment dihedral about the link for each of those two)
    if resname not in _planar_link_cache:
        coords, bonded = fragment_template(resname, frag_template)
        hydrogens = {}
        for group in frag_template.values():
            for bead in group.values():
                for number, atom in bead.items():
                    if gen.is_hydrogen(atom['atom']):
                        hydrogens[number] = np.array(atom['coord'], dtype=float)
        h_numbers = list(hydrogens)
        h_tree = cKDTree(np.array([hydrogens[n] for n in h_numbers])) if h_numbers else None
        attached = lambda n: [h_numbers[i] for i in h_tree.query_ball_point(coords[n], 1.25 * g_var.sf)] if h_tree is not None else []
        out = []
        for a in coords:
            if not frag_atom_name(frag_template, a).startswith('C') or len(bonded[a]) + len(attached(a)) != 3:
                continue
            for b in bonded[a]:
                if frag_atom_name(frag_template, b)[0] not in 'NO' or len(bonded[b]) != 2:
                    continue
                l = next(iter(bonded[b] - {a}))
                movers = sorted(bonded[a] - {b})
                targets = [dihedral(coords[i], coords[a], coords[b], coords[l]) for i in movers]
                out.append((a, b, l, movers, [h for i in movers for h in attached(i)], targets))
        _planar_link_cache[resname] = out
    return _planar_link_cache[resname]

def correct_planar_links(residue, resname, frag_template):
#### an amide or ester whose carbonyl and N/O substituent were fitted in different beads can be placed
#### cis or twisted, and minimisation cannot turn the attached chain back through the barrier. The
#### carbonyl's two other substituents (the C=O oxygen and the alpha carbon, with its hydrogens) are
#### turned about the link to the fragment's (trans) dihedral; only the alpha carbon's onward bond is
#### stretched, which the minimisation then relaxes. CISRES holds the link afterwards
    for a, b, l, movers, carried, targets in planar_links(resname, frag_template):
        if any(x not in residue for x in [a, b, l] + movers + carried):
            continue
        X = {x: np.array(residue[x]['coord'], dtype=float) for x in [a, b, l] + movers + carried}
        current = dihedral(X[movers[0]], X[a], X[b], X[l])
        delta = math.atan2(math.sin(targets[0] - current), math.cos(targets[0] - current))
        if abs(delta) < math.radians(30):
            continue
        axis = X[a] - X[b]
        best = None
        for angle in (delta, -delta):
            rot = rotation_about_axis(axis, angle)
            moved = {x: (X[x] - X[a]).dot(rot) + X[a] for x in movers + carried}
            error = abs(math.atan2(math.sin(targets[0] - dihedral(moved[movers[0]], X[a], X[b], X[l])),
                                   math.cos(targets[0] - dihedral(moved[movers[0]], X[a], X[b], X[l]))))
            if best is None or error < best[0]:
                best = (error, moved)
        for x, coord in best[1].items():
            residue[x]['coord'] = coord
    return residue

def rotation_about_axis(axis, angle):
#### row-vector rotation matrix (used as coord.dot(matrix)) for a right-handed turn about axis
    x, y, z = axis / np.linalg.norm(axis)
    c, s = math.cos(angle), math.sin(angle)
    C = 1 - c
    return np.array([[c + x*x*C, x*y*C - z*s, x*z*C + y*s],
                     [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
                     [z*x*C - y*s, z*y*C + x*s, c + z*z*C]]).T

def correct_group_spin(residue, group_fit, template):
#### spins a newly placed group about the axis through its connection atoms so the
#### dihedrals across its bond to the already placed atoms match the fragment template
    t_coords, bonded = template
    group_atoms = {atom for bead in group_fit.values() for atom in bead if atom in t_coords}
    placed = {atom for atom in residue if atom in t_coords and atom not in group_atoms}
    if not group_atoms or not placed:
        return residue
    joins = [(a, b) for a in group_atoms for b in bonded[a] if b in placed]
    if len(joins) != 1:      #### only the single-attachment case has a free spin
        return residue
    a, b = joins[0]
    onward = [x for x in group_atoms for y in bonded[x] if y not in group_atoms and y not in placed]
    pivot = residue[a]['coord']
    if onward:
        axis = residue[onward[0]]['coord'] - pivot
    else:
        axis = np.mean([residue[x]['coord'] for x in group_atoms], axis=0) - pivot
    if np.linalg.norm(axis) < 1e-3 * g_var.sf:
        return residue
    quartets = [(i, a, b, l) for i in bonded[a] if i in group_atoms
                for l in bonded[b] if l in placed]
    if not quartets:
        return residue
    targets = np.array([dihedral(*[t_coords[x] for x in q]) for q in quartets])
    best_angle, best_score = 0.0, None
    for angle in np.radians(np.arange(0, 360, 5)):
        rot = rotation_about_axis(axis, angle)
        current = np.array([dihedral((residue[i]['coord'] - pivot).dot(rot) + pivot,
                                     residue[j]['coord'], residue[k]['coord'], residue[l]['coord'])
                            for i, j, k, l in quartets])
        score = np.sum(1 - np.cos(current - targets))
        if best_score is None or score < best_score - 1e-9:
            best_angle, best_score = angle, score
    if best_angle:
        rot = rotation_about_axis(axis, best_angle)
        for bead in group_fit.values():
            for atom in bead:
                residue[atom]['coord'] = (residue[atom]['coord'] - pivot).dot(rot) + pivot
    return residue

## Chirality across group junctions.
## Each rigid group keeps the fragment's handedness, but a stereocentre whose substituents sit in
## different groups takes its handedness from where the groups land. CG beads carry no chirality, so
## in a CG simulation the beads around such a junction can drift into the mirror-image arrangement
## (Martini KLA: about half the molecules at the Kdo II anomeric carbon), and every group then fits
## its own beads well while the centre between them comes out inverted - beyond what restraints and
## minimisation can undo. Such a junction is mended by re-placing the group on the side away from the
## middle of the molecule from the fragment's own geometry relative to its neighbour: its bead fit is
## given up, the chirality is kept.

_junction_cache = {}

def group_junctions(resname, frag_template, groups):
#### groups: {group: [atom numbers]} as placed. Returns (order, checks): the groups in breadth-first
#### order from the middle of the molecule (the group whose removal leaves the smallest largest part),
#### each with its parent and descendants, and for each child group the stereocentres at its bonds to
#### the parent: (centre, three neighbours on the two groups' side, the fragment's handedness)
    if resname in _junction_cache:
        return _junction_cache[resname]
    coords, bonded = fragment_template(resname, frag_template)
    parents_h = hydrogen_parents(resname, frag_template)
    frag_xyz = {number: np.array(atom['coord'], dtype=float) for group in frag_template.values()
                for bead in group.values() for number, atom in bead.items()}
    of_group = {atom: g for g, atoms in groups.items() for atom in atoms}
    links = {g: set() for g in groups}
    bonds = {}
    for a in coords:
        for b in bonded[a]:
            if a in of_group and b in of_group and of_group[a] != of_group[b]:
                links[of_group[a]].add(of_group[b])
                bonds.setdefault(frozenset((of_group[a], of_group[b])), set()).add((a, b))
    size = {g: len(atoms) for g, atoms in groups.items()}
    def largest_part_without(g):
        seen, worst = {g}, 0
        for start in links[g]:
            if start in seen:
                continue
            part, stack = 0, [start]
            seen.add(start)
            while stack:
                x = stack.pop()
                part += size[x]
                for y in links[x]:
                    if y not in seen:
                        seen.add(y)
                        stack.append(y)
            worst = max(worst, part)
        return worst
    root = min(groups, key=lambda g: (largest_part_without(g), -size[g], str(g)))
    parent, order, queue = {root: None}, [], [root]
    while queue:
        g = queue.pop(0)
        order.append(g)
        for h in sorted(links[g], key=str):
            if h not in parent:
                parent[h] = g
                queue.append(h)
    descendants = {g: [] for g in parent}
    for g in reversed(order):
        if parent[g] is not None:
            descendants[parent[g]] += [g] + descendants[g]
    hydrogens_of = {}
    for h, p in parents_h.items():
        hydrogens_of.setdefault(p, []).append(h)
    checks = {}
    for g in order:
        p = parent[g]
        if p is None:
            continue
        side = set(groups[g]) | set(groups[p])
        found = []
        for a, b in bonds.get(frozenset((g, p)), ()):
            for c in (a, b):
                hs = hydrogens_of.get(c, [])
                if not frag_atom_name(frag_template, c).startswith('C') or len(bonded[c]) + len(hs) != 4 or len(hs) > 1:
                    continue
                near = sorted(n for n in bonded[c] if n in side) + [h for h in hs if h in side]
                if len(near) < 3:
                    continue
                near = near[:3]
                v = [frag_xyz[x] - frag_xyz[c] for x in near]
                volume = np.dot(v[0], np.cross(v[1], v[2])) / np.prod([np.linalg.norm(x) for x in v])
                if abs(volume) > 0.1 and (c, tuple(near)) not in [(f[0], f[1]) for f in found]:
                    found.append((c, tuple(near), np.sign(volume)))
        checks[g] = found
    _junction_cache[resname] = ({g: (parent[g], descendants[g]) for g in order}, order, checks, frag_xyz)
    return _junction_cache[resname]

def correct_group_chirality(residue, resname, frag_template, groups):
#### residue: {atom number: atom} with every group placed. Walks out from the middle of the molecule;
#### where a stereocentre at a group's bond to its parent has the wrong handedness, the group is put
#### where the fragment has it relative to its parent (superposing the fragment on the parent's placed
#### atoms) and everything beyond it follows rigidly
    family, order, checks, frag_xyz = group_junctions(resname, frag_template, groups)
    mended = 0
    for g in order:
        p, beyond = family[g]
        if p is None or not checks.get(g):
            continue
        wrong = False
        for c, near, sign in checks[g]:
            if any(x not in residue for x in (c,) + near):
                continue
            v = [np.array(residue[x]['coord'], dtype=float) - np.array(residue[c]['coord'], dtype=float) for x in near]
            if np.sign(np.dot(v[0], np.cross(v[1], v[2]))) != sign:
                wrong = True
                break
        if not wrong:
            continue
        anchor = [x for x in groups[p] if x in residue and x in frag_xyz]
        moving = [x for x in groups[g] if x in residue and x in frag_xyz]
        if len(anchor) < 3 or not moving:
            continue
        A = np.array([frag_xyz[x] for x in anchor])
        B = np.array([residue[x]['coord'] for x in anchor], dtype=float)
        ac, bc = A.mean(axis=0), B.mean(axis=0)
        rot = kabsch_rotate(A - ac, B - bc)
        old = np.array([residue[x]['coord'] for x in moving], dtype=float)
        new = np.array([(frag_xyz[x] - ac).dot(rot) + bc for x in moving])
#### the groups beyond follow the same rigid move that took this group from its fit to its new place
        oc, nc = old.mean(axis=0), new.mean(axis=0)
        follow = kabsch_rotate(old - oc, new - nc)
        for x, xyz in zip(moving, new):
            residue[x]['coord'] = xyz
        for h in beyond:
            for x in groups[h]:
                if x in residue:
                    residue[x]['coord'] = (np.array(residue[x]['coord'], dtype=float) - oc).dot(follow) + nc
        mended += 1
    return residue, mended

### end fragment information
### linked residue connectivity start

def BB_connectivity(at_connections,cg_connections, cg_residues, at_residues, residue_number, BB_bead, resname):
    con_atoms = {}
    for atom in at_residues:
        # resname = at_residues[atom]['res_type']
        if at_residues[atom]['atom'] in g_var.res_top[resname]['CONNECT'][BB_bead]['atom']:
            con_atoms[g_var.res_top[resname]['CONNECT'][BB_bead]['atom'].index(at_residues[atom]['atom'])]=atom
    new_chain=False
    for con in con_atoms:
        con_resid = residue_number+g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con]
        if con_resid in cg_residues:
            xyz_cur = cg_residues[residue_number][BB_bead]['coord']
            if g_var.res_top[resname]['CONNECT'][BB_bead]['Con_Bd'][con] in cg_residues[con_resid]:
                xyz_con = cg_residues[con_resid][g_var.res_top[resname]['CONNECT'][BB_bead]['Con_Bd'][con]]['coord']
                if gen.calculate_distance(xyz_con, xyz_cur) < 6:
                    cg_connections.append(xyz_con)
                    at_connections.append(at_residues[con_atoms[con]]['coord'])
                else:
                    if g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con] > 0:
                        new_chain=True
            else:
                    if g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con] > 0:
                        new_chain=True
        else:
            if g_var.res_top[resname]['CONNECT'][BB_bead]['dir'][con] > 0:
                new_chain=True
    return at_connections,cg_connections, new_chain

################################################################### Merged system

def merge_individual_chain_pdbs(file, end, res_type):
#### reads in each chain into merge list
    merge, merged_coords = [],[]
    count = 0
    restraint_count = -1
    for chain in range(0,g_var.system[res_type]):
        merge_temp = []
        if os.path.exists(file+'_'+str(chain)+end):
            with open(file+'_'+str(chain)+end, 'r') as pdb_input:
                merge_temp += read_in.filter_input(pdb_input.readlines(), False)
        else:
            raise InputError('cannot find chain: '+file+'_'+str(chain)+end)
        if res_type+'_aligned' in file:  
            count, restraint_count = at_mod_p.create_disres(merge_temp, chain, file, count, restraint_count)
        merge, merge_coords = fix_chirality(merge,merge_temp,merged_coords, res_type)   
    if res_type+'_aligned' not in file:
        coords, index_conversion = index_conversion_generate(merge, merge_coords)
    if 'aligned' in file.split('/')[-1]:
        write_pdb(merge, merge_coords, {}, g_var.working_dir+'PROTEIN/PROTEIN_aligned_merged.pdb')
    else:
        write_pdb(merge, coords,index_conversion, g_var.working_dir+res_type+'/'+res_type+'_de_novo_merged.pdb')

def index_conversion_generate(merge, merge_coords):
    index_conversion = {}
    coords, dummies = [], []
    for at_val, atom in enumerate(merge):
        if not atom['atom_name'].startswith('M'):
            index_conversion[at_val] = len(coords)
            coords.append(merge_coords[at_val])
        else:
            dummies.append(merge_coords[at_val])
    coords = check_atom_overlap(coords, dummies)
    return coords, index_conversion

def write_pdb(merge, merge_coords, index_conversion, write_file):
#### creates merged pdb and writes chains to it
    if not os.path.exists(write_file):
        pdb_output=gen.create_pdb(write_file)
        atom_counter=0
        for line_val, line in enumerate(merge):
            if line_val in index_conversion:
                x, y, z = gen.trunc_coord(merge_coords[index_conversion[line_val]])
            else:
                if 'coord' in line:
                    x, y, z = gen.trunc_coord(line['coord'])
                else:
                    x, y, z = gen.trunc_coord([line['x'],line['y'],line['z']])
            if atom_counter >= 99_999:
                atom_counter=1
            else:
                atom_counter+=1      
            pdb_output.write(g_var.pdbline%((atom_counter, line['atom_name'], line['residue_name'],' ',line['residue_id'],\
            x,y,z,1,0))+'\n')
        pdb_output.close()

def merge_system_pdbs(protein):
    os.chdir(g_var.merged_directory)
#### create merged pdb 
    if not os.path.exists(g_var.merged_directory+'merged_cg2at'+protein+'.pdb'):
        merge=[]
        merge_coords=[]
    #### run through every residue type in cg_residues
        done = []
        if 'novo' in protein:
            print('checking for atom overlap in : '+protein[1:])
        for segment, residue_type in enumerate(g_var.system):
            if residue_type not in done:
                if residue_type not in ['PROTEIN', 'OTHER']:
                        start_num = len(merge_coords)
                        merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, 
                                                                g_var.working_dir+residue_type+'/'+residue_type+'_merged.pdb')
                        g_var.np_blocks[residue_type]=[start_num,len(merge_coords)]
                elif residue_type == 'OTHER':
                    merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, 
                                                            g_var.working_dir+residue_type+'/'+residue_type+'_de_novo_merged.pdb')
                elif residue_type == 'PROTEIN':
                    merge, merge_coords = read_in_merged_pdbs(merge, merge_coords, 
                                                            g_var.working_dir+residue_type+'/'+residue_type+protein+'_merged.pdb')     
                done.append(residue_type)
        index_conversion = {}
        if 'novo' in protein:
            coords, index_conversion = index_conversion_generate(merge, merge_coords)
            write_pdb(merge, coords, index_conversion, g_var.merged_directory+'merged_cg2at'+protein+'.pdb')
        else:
            write_pdb(merge, merge_coords, {}, g_var.merged_directory+'merged_cg2at'+protein+'.pdb')

def read_in_merged_pdbs(merge, merge_coords, location):
    if os.path.exists(location):
    #### opens pdb files and writes straight to merged_cg2at pdb
        with open(location, 'r') as pdb_input:
            read_in_atoms = read_in.filter_input(pdb_input.readlines(), False)
            merge_coords += [[line_sep['x'],line_sep['y'],line_sep['z']] for line_sep in read_in_atoms ]
            return merge+read_in_atoms, merge_coords
    else:
        raise InputError('cannot find minimised residue: '+location) 

def check_overlap_chain(chain, input,res_type):
    if not os.path.exists(g_var.working_dir+res_type+'/'+res_type+'_'+input+str(chain)+'_gmx_checked.pdb'):
        lines, coords = read_in_merged_pdbs([], [], res_type+'_'+input+str(chain)+'_gmx.pdb')
        if res_type == 'PROTEIN':
            coords = at_mod_p.correct_amide_h(lines, coords)

        coords, index_conversion = index_conversion_generate(lines, coords)
        write_pdb(lines, coords, index_conversion, g_var.working_dir+res_type+'/'+res_type+'_'+input+str(chain)+'_gmx_checked.pdb')

## fix chirality errors

def fetch_chiral_coord(merge_temp, residue_type):
    chiral_atoms={}
    coord=[]
    for atom in range(len(merge_temp)):
        if residue_type in ['PROTEIN', 'OTHER']:
            resname = gen.check_alternate_resname(merge_temp[atom]['residue_name'])
        else:
            resname = gen.check_alternate_resname(residue_type)
        if len(g_var.res_top[resname]['CHIRAL']['atoms']) > 0:
            if merge_temp[atom]['atom_name'] in g_var.res_top[resname]['CHIRAL']['atoms']:
                if merge_temp[atom]['residue_id'] not in chiral_atoms:
                    chiral_atoms[merge_temp[atom]['residue_id']]={}
                chiral_atoms[merge_temp[atom]['residue_id']][merge_temp[atom]['atom_name']]=atom
        coord.append(np.array([merge_temp[atom]['x'],merge_temp[atom]['y'],merge_temp[atom]['z']]))
    return chiral_atoms, coord


def get_atom_move(merge_temp, resname, residue, chiral_group, chiral_atoms):
    stat = merge_temp[chiral_atoms[residue][chiral_group]].copy()
    atom_move = {'stat':np.array([stat['x'],stat['y'],stat['z']]), 'm':'', 'c1':'', 'c2':'', 'c3':''}
    for chir_atom in atom_move:
        if chir_atom != 'stat':
            test = merge_temp[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group][chir_atom]]].copy()
            atom_move[chir_atom]= np.array([test['x'],test['y'],test['z']])
            if gen.calculate_distance(atom_move['stat'], atom_move[chir_atom]) > 10:
                atom_move[chir_atom] = np.array(read_in.brute_mic(atom_move['stat'],atom_move[chir_atom]))
    return atom_move


def fix_chirality(merge, merge_temp, merged_coords, residue_type):
#### fixes chiral groups
    chiral_atoms, coord = fetch_chiral_coord(merge_temp, residue_type)
    for residue in chiral_atoms:
        if residue_type in ['PROTEIN', 'OTHER']:
            for atom in chiral_atoms[residue]:
                resname = merge_temp[chiral_atoms[residue][atom]]['residue_name']
                resname = gen.check_alternate_resname(resname)
                break
        else:
            resname = gen.check_alternate_resname(residue_type)

        for chiral_group in g_var.res_top[resname]['CHIRAL']:
            if chiral_group != 'atoms':
                atom_move = get_atom_move(merge_temp, resname, residue, chiral_group, chiral_atoms)

                S_M = atom_move['m'] - atom_move['stat']
                rotation = align_to_vector(S_M, [0,0,1])
                c1_coord = (atom_move['c1'] - atom_move['stat']).dot(rotation)
                c2_coord = (atom_move['c2'] - atom_move['stat']).dot(rotation)
                c3_coord = (atom_move['c3'] - atom_move['stat']).dot(rotation)

                if gen.angle_clockwise(c1_coord[0:2], c2_coord[0:2]) > gen.angle_clockwise(c1_coord[0:2], c3_coord[0:2]):
                    for ax_val, axis in enumerate(['x', 'y', 'z']):
                        merge_temp[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group]['m']]][axis] = merge_temp[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group]['m']]][axis] - (3*S_M[ax_val])   
                        merge_temp[chiral_atoms[residue][chiral_group]][axis] = merge_temp[chiral_atoms[residue][chiral_group]][axis] - (S_M[ax_val])        
                    coord[chiral_atoms[residue][g_var.res_top[resname]['CHIRAL'][chiral_group]['m']]] -=  (2*S_M) #move_coord -
                    coord[chiral_atoms[residue][chiral_group]] -=  (0.25*S_M) 
    merge+=merge_temp
    merged_coords+=coord
    return merge, merged_coords

def get_chiral_non_carbonyl(at):
    chiral = {}
    for atom in at:
        if at[atom]['atom'] in g_var.res_top[at[atom]['res_type']]['CHIRAL']['atoms']:
            chiral[at[atom]['atom']] = atom
    return chiral

def correct_chiral_atoms(at, chiral):
    for chiral_group in g_var.res_top[at[1]['res_type']]['CHIRAL']:
        if chiral_group != 'atoms':            
            atom_list = [chiral[chiral_group]]
            for chiral_atoms in ['c1','c2','c3']:
                atom_list.append(chiral[g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group][chiral_atoms]])
            cross_vector_chiral = find_cross_vector( [at[atom_list[3]]['coord'], at[atom_list[2]]['coord'], at[atom_list[1]]['coord']])
            norm_vector_chiral = noramlised_vector(at[chiral[g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m']]]['coord'], at[atom_list[0]]['coord'])
            if not np.dot(cross_vector_chiral, norm_vector_chiral) > 0.5:
                if g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m'] in chiral:
                    at[atom_list[0]]['coord'] = at[atom_list[0]]['coord'] + cross_vector_chiral*0.5
                    m = chiral[g_var.res_top[at[1]['res_type']]['CHIRAL'][chiral_group]['m']]
                    at[m]['coord'] = at[m]['coord'] + cross_vector_chiral*2
    return at

## chirality end
## hydrogen orientaion checker

def check_hydrogens(residue):
#### finds the connecting carbons and their associated carbons [carbon atom, hydrogen ref number, connecting ref number]    
    for atom_num, atom in enumerate(residue):
        resname=residue[atom]['res_type']
        break
    for atom in g_var.hydrogen[resname]:
        h_coord = []
        for group in g_var.sorted_connect[resname]:
            if atom in g_var.sorted_connect[resname][group]:
                for hydrogen in g_var.hydrogen[resname][atom]:
                    h_coord.append(residue[hydrogen]['coord'])
                h_com=np.mean(np.array(h_coord), axis=0)
                for heavy_bond in g_var.heavy_bond[resname][atom]:
                    for group_check in g_var.sorted_connect[resname]:
                        if heavy_bond in g_var.sorted_connect[resname][group_check] and group_check != group:
                            skip = False
                            if heavy_bond in residue:
                                con_heavy_atom = heavy_bond
                                con_heavy_atom_co =  residue[con_heavy_atom]['coord']
                            else:
                                skip = True
                if not skip:
                #### vector between H COM and bonded carbon 
                    vector=np.array([h_com[0]-residue[atom]['coord'][0],h_com[1]-residue[atom]['coord'][1],h_com[2]-residue[atom]['coord'][2]])
                    h_com_f=h_com+vector*2
                    d1 = gen.calculate_distance(h_com, con_heavy_atom_co)    
                    d2 = gen.calculate_distance(h_com_f, con_heavy_atom_co)   
                    if d2 < d1:
                        for h_at in g_var.hydrogen[resname][atom]:
                            residue[h_at]['coord']=residue[h_at]['coord']-vector*2
    return residue

## check for threaded lipids
def check_ringed_lipids(protein):
    print('Checking for ringed lipids')
    if not os.path.exists(g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb'):
        if not os.path.exists(g_var.merged_directory+'merged_cg2at_threaded.pdb'):
            os.chdir(g_var.merged_directory)
            merge, merge_coords = read_in_merged_pdbs([], [], protein)
            ringed=False
            lipid_atoms = []
            with open(g_var.merged_directory+'threaded_lipids.dat', 'w') as ring_ouput:
                for at_val, atom in enumerate(merge): 
                    resname = get_np_resname(at_val)
                    if resname in g_var.np_residues:
                        offset = fetch_start_of_residue_np(at_val, resname)
                        if atom['atom_number']-offset in g_var.heavy_bond[resname]:
                            for at_bond in g_var.heavy_bond[resname][atom['atom_number']-offset]:
                                at_bond -=1
                                if merge[at_bond+offset]['atom_number'] > merge[at_val]['atom_number']:
                                    merge[at_bond+offset]['x'], merge[at_bond+offset]['y'], merge[at_bond+offset]['z'] = np.array(read_in.brute_mic(merge_coords[at_val],merge_coords[at_bond+offset]))
                                    merge_coords[at_bond+offset] = merge[at_bond+offset]['x'], merge[at_bond+offset]['y'], merge[at_bond+offset]['z']
                                    dist = gen.calculate_distance(merge_coords[at_val], merge_coords[at_bond+offset])
                                    if 2 < dist < 6:
                                        lipid_atoms.append([at_val, at_bond+offset, (np.array(merge_coords[at_val])+np.array(merge_coords[at_bond+offset]))/2])
                                        ring_ouput.write('{0:6}{1:6}{2:2}{3:4}{4:2}{5:7}{6:5}{7:5}{8:5}{9:5}{10:5}{11:5}\n'.format(
                                                            'distance: ',str(np.round(dist,2)),'residue: ', merge[at_val]['residue_name'], merge[at_val]['residue_id'],
                                                            ' atom_1: ', merge[at_val]['atom_name'],
                                                            'atom_2: ', merge[at_bond+offset]['atom_name'], 'rough line num: ', at_val, at_bond+offset))
                                        ringed = True
        if ringed or os.path.exists(g_var.merged_directory+'merged_cg2at_threaded.pdb'):
            print('Found '+str(len(lipid_atoms))+' abnormal bonds, now attempting to fix.')
            print('See this file for a complete list: '+g_var.merged_directory+'threaded_lipids.dat')
            fix_threaded_lipids(lipid_atoms, merge, merge_coords)
        else:
            gen.file_copy_and_check(g_var.merged_directory+'MIN/merged_cg2at_de_novo_minimised.pdb', g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb')

def fetch_start_of_residue(at, merge):
    count = at
    while True:
        if merge[count]['residue_id'] != merge[at]['residue_id']:
            count+=1
            return count
        else:
            count-=1

def fetch_start_of_residue_np(atom, resname):
            prev=g_var.np_blocks[resname][0]
            range_np_block = g_var.np_blocks[resname][1]-g_var.np_blocks[resname][0]
            for range_np_block_ind in range(g_var.np_blocks[resname][0],g_var.np_blocks[resname][1]+1, int(range_np_block/g_var.system[resname])):
                if prev <= atom < range_np_block_ind:
                    return prev
                prev = range_np_block_ind

def get_np_resname(atom):
    for res_check in g_var.np_blocks:
        if g_var.np_blocks[res_check][0] <= atom < g_var.np_blocks[res_check][1]:
            return res_check 

def fix_threaded_lipids(lipid_atoms, merge, merge_coords):
    if not os.path.exists(g_var.merged_directory+'merged_cg2at_threaded.pdb'):
        tree = cKDTree(merge_coords)         
        for threaded in lipid_atoms:  
            resname = get_np_resname(threaded[0]) 
            
            atoms = tree.query_ball_point(threaded[2], r=3)
            P_count = None
            for at in atoms:
                if merge[at]['residue_id'] != merge[threaded[0]]['residue_id']:
                    P_count = fetch_start_of_residue(at, merge)
                    break
#### no other residue next to the stretched bond: it is not threaded through anything, just stretched,
#### and is left to the NVT
            if P_count is None:
                continue

            NP_count = fetch_start_of_residue_np(threaded[0], resname)
            bb = []
            for at in merge[P_count:]:
                if at['residue_id'] != merge[P_count]['residue_id']:
                    break
                if at['atom_name'] in g_var.res_top[at['residue_name']]['ATOMS']:
                    bb.append([at['x'], at['y'], at['z']])
#### only a lipid threaded through a residue with backbone atoms (a protein) can be pushed out; a
#### stretched bond next to anything else (e.g. inside a glycolipid) is left to the minimisation
            if not bb:
                continue
            bb = np.mean(np.array(bb), axis=0)
            if np.linalg.norm(threaded[2]-bb) < 1e-6:
                continue
            BB_M3 = (threaded[2]-bb)/np.linalg.norm((threaded[2]-bb))      
            
            for heavy_atom in threaded[:2]:
                merge_coords[heavy_atom] += BB_M3*3 
                merge[heavy_atom]['x'], merge[heavy_atom]['y'], merge[heavy_atom]['z'] = merge_coords[heavy_atom]
                for hydrogen in g_var.hydrogen[resname].get(heavy_atom-NP_count+1, []):
                    merge_coords[NP_count+hydrogen-1] += BB_M3*3 
                    merge[NP_count+hydrogen-1]['x'], merge[NP_count+hydrogen-1]['y'], merge[NP_count+hydrogen-1]['z'] = merge_coords[NP_count+hydrogen-1]

        coords, index_conversion = index_conversion_generate(merge, merge_coords)
        write_pdb(merge, coords, index_conversion, g_var.merged_directory+'merged_cg2at_threaded.pdb')

    if not os.path.exists(g_var.merged_directory+'MIN/merged_cg2at_threaded_minimised.pdb'):
        gro.minimise_merged_pdbs('_threaded')
    gen.file_copy_and_check(g_var.merged_directory+'MIN/merged_cg2at_threaded_minimised.pdb', g_var.merged_directory+'checked_ringed_lipid_de_novo.pdb')
