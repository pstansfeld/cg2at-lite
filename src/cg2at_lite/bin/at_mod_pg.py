#!/usr/bin/env python3

#### Peptidoglycan (PG) layers built by MemPrO/Insane (Martini 3): residues NAG, NAM and the peptide
#### stems SPEP (pentapeptide), UPEP (tetrapeptide) and UUPEP (tripeptide), all in one PGL molecule
#### that is bonded across the periodic box.
####
#### Each residue is fitted from its fragment (the free, unlinked state) like any non-protein residue.
#### The links are then applied as patches (as CHARMM PRES): beta1,4 glycosidic bonds (14bb charges),
#### the MurNAc lactyl amide to the stem, and 3-4 (D-Ala4 C=O) or 3-3 (mDAP L-centre C=O) cross-links
#### to the acceptor's mDAP D-centre amine. One PGL moleculetype is written with the bonded terms
#### generated from the bond graph (angles, proper dihedrals, 1-4 pairs, amide/carboxylate impropers
#### and CMAP), so GROMACS takes every parameter from the CHARMM36 atom types.

import os
import copy
import numpy as np
from cg2at_lite.bin import gen, g_var, at_mod, at_mod_np, gro
from cg2at_lite.bin.exceptions import InputError, TopologyError

PG_RESIDUES = ['NAG', 'NAM', 'SPEP', 'UPEP', 'UUPEP']
STEMS = ['SPEP', 'UPEP', 'UUPEP']
#### donor carbonyl of each stem type: (bead in the CG model, carbonyl C, free C-terminal O kept, O deleted)
DONOR = {'UPEP': ('B8', 'C4', 'OT41', 'OT42', 'O4'), 'UUPEP': ('B5', 'C3', 'OT31', 'OT32', 'O3')}
#### largest CG bead distances (A) accepted when links are inferred from coordinates
LINK_CUTOFF = {'glycosidic': 4.2, 'stem': 4.2, 'cross': 4.6}

#### patches: atoms deleted, renamed and retyped (type, charge) when a link is made
PATCH = {
    'c1':       {'delete': ['O1', 'HO1'], 'rename': {}, 'set': {'C1': ('CC3162', 0.29)}},
    'o4':       {'delete': ['HO4'], 'rename': {}, 'set': {'C4': ('CC3161', 0.09), 'O4': ('OC301', -0.36)}},
    'acceptor': {'delete': ['HZ32', 'HZ33'], 'rename': {'HZ31': 'HZ3'},
                 'set': {'CE3': ('CTD1', 0.07), 'HE3': ('HB1', 0.09), 'NZ3': ('NH1', -0.47), 'HZ3': ('H', 0.31)}},
    'donor_UPEP':  {'delete': ['OT42'], 'rename': {'OT41': 'O4'}, 'set': {'C4': ('C', 0.51), 'O4': ('O', -0.51)}},
    'donor_UUPEP': {'delete': ['OT32'], 'rename': {'OT31': 'O3'}, 'set': {'C3': ('C', 0.51), 'O3': ('O', -0.51)}},
}


def pg_present():
    return any(residue in g_var.cg_residues for residue in PG_RESIDUES)


def collect_pg_residues():
#### moves the PG residues out of cg_residues into one ordered list (CG file order = PGL residue order)
    residues = []
    for resname in PG_RESIDUES:
        for count, beads in g_var.cg_residues.pop(resname, {}).items():
            residues.append((count, resname, beads))
    g_var.pg_residues = [(resname, beads) for count, resname, beads in sorted(residues, key=lambda x: x[0])]
    print('\nPeptidoglycan layer: '+str(len(g_var.pg_residues))+' residues ('+
          ', '.join(str(sum(r == res for r, b in g_var.pg_residues))+' '+res for res in PG_RESIDUES)+')')
    check_relaxed()


def min_image(vector):
#### shortest periodic image of a vector (A), triclinic boxes included
    frac = np.dot(g_var.r_b_inv, vector)
    return vector - np.dot(g_var.r_b_vec, np.round(frac))


def bead_distance(a, b):
    return np.linalg.norm(min_image(np.array(a) - np.array(b)))


def check_relaxed():
#### a freshly built Insane layer has its sugar beads at about half their bonded distances (e.g. NAG
#### B1-B3 ~1.5 A for b0 3.0 A); fitted fragments then overlap and the minimisation cannot recover the
#### rings. Such input is refused: relax the CG system first (minimisation and a short MD).
    distances = [bead_distance(beads['B1']['coord'], beads['B3']['coord'])
                 for resname, beads in g_var.pg_residues if resname in ['NAG', 'NAM'] and 'B1' in beads and 'B3' in beads]
    if distances and np.median(distances) < 2.2:
        raise InputError('the peptidoglycan layer looks unrelaxed: the median sugar B1-B3 bead distance is '+
                         str(np.round(np.median(distances), 2))+' A (Martini b0 ~2.9-3.0 A).\n'
                         'Minimise the CG system and run a short CG MD (periodic-molecules = yes) before back-mapping.')

######################################################################## links

def read_cg_topology(file_name):
#### PGL moleculetype of a CG topology: the itp itself, or a .top whose #include files contain it.
#### Returns (atoms {index: (resnr, resname, bead, charge)}, bonds [(i, j)])
    if file_name.endswith('.top'):
        directory = os.path.dirname(os.path.abspath(file_name))
        with open(file_name) as top:
            includes = [line.split('"')[1] for line in top if line.strip().startswith('#include') and '"' in line]
        for include in includes:
            path = include if os.path.isabs(include) else os.path.join(directory, include)
#### collect_input copies a relative include next to the .top under its own name
            if not os.path.exists(path):
                path = os.path.join(directory, os.path.basename(include))
            if os.path.exists(path):
                found = read_cg_topology(path)
                if found[0]:
                    return found
        raise InputError('cannot find the PGL moleculetype in the files included by: '+file_name)
    atoms, bonds, section, molecule = {}, [], None, None
    with open(file_name) as itp:
        for line in itp:
            line_sep = line.split(';')[0].split()
            if not line_sep or line_sep[0].startswith('#'):
                continue
            if line_sep[0].startswith('['):
                section = gen.strip_header(line).lower()
                continue
            if section == 'moleculetype':
                molecule = line_sep[0]
            elif molecule == 'PGL' and section == 'atoms':
                atoms[int(line_sep[0])] = (int(line_sep[2]), line_sep[3], line_sep[4], float(line_sep[6]))
            elif molecule == 'PGL' and section == 'bonds':
                bonds.append((int(line_sep[0]), int(line_sep[1])))
    return atoms, bonds


def classify_link(res_a, bead_a, res_b, bead_b):
#### (donor residue, atom, acceptor residue, atom) of a CG bond between two PG residues, or None
    names = g_var.pg_residues
    ra, rb = names[res_a][0], names[res_b][0]
    for (r1, b1, n1, r2, b2, n2) in ((res_a, bead_a, ra, res_b, bead_b, rb), (res_b, bead_b, rb, res_a, bead_a, ra)):
        if n1 in ['NAG', 'NAM'] and b1 == 'B1' and n2 in ['NAG', 'NAM'] and b2 == 'B4' and n1 != n2:
            return (r1, 'C1', r2, 'O4', 'glycosidic')
        if n1 == 'NAM' and b1 == 'B5' and n2 in STEMS and b2 == 'B1':
            return (r1, 'CE3', r2, 'N1', 'stem')
        if n1 in DONOR and b1 == DONOR[n1][0] and n2 in STEMS and b2 == 'B7':
            return (r1, DONOR[n1][1], r2, 'NZ3', 'cross')
#### MemPrO before df7a03b joined 3-3 links bead 7 to bead 7: the tripeptide is the donor
    if ra in STEMS and rb in STEMS and bead_a == 'B7' and bead_b == 'B7' and 'UUPEP' in [ra, rb]:
        donor, acceptor = (res_a, res_b) if ra == 'UUPEP' else (res_b, res_a)
        print('Note: 3-3 cross-link given as bead 7 - bead 7 (older MemPrO); the tripeptide is taken as the donor')
        return (donor, 'C3', acceptor, 'NZ3', 'cross')
    return None


def links_from_topology(file_name):
    atoms, bonds = read_cg_topology(file_name)
    if not atoms:
        raise InputError('no PGL moleculetype found in: '+file_name)
    resnrs = sorted(set(atom[0] for atom in atoms.values()))
    if len(resnrs) != len(g_var.pg_residues):
        raise InputError('the CG topology has '+str(len(resnrs))+' PGL residues, the CG structure '+str(len(g_var.pg_residues)))
    order = {resnr: index for index, resnr in enumerate(resnrs)}
    cg_charge = {}
    for index, (resnr, resname, bead, charge) in atoms.items():
        if resname != g_var.pg_residues[order[resnr]][0]:
            raise InputError('PGL residue '+str(resnr)+' is '+resname+' in the CG topology but '+
                             g_var.pg_residues[order[resnr]][0]+' in the CG structure')
        cg_charge[order[resnr]] = cg_charge.get(order[resnr], 0) + charge
    links = []
    for i, j in bonds:
        res_a, res_b = order[atoms[i][0]], order[atoms[j][0]]
        if res_a == res_b:
            continue
        link = classify_link(res_a, atoms[i][2], res_b, atoms[j][2])
        if link is None:
            raise TopologyError('unrecognised bond between PGL residues '+atoms[i][1]+str(atoms[i][0])+' '+atoms[i][2]+
                                ' and '+atoms[j][1]+str(atoms[j][0])+' '+atoms[j][2])
        links.append(link)
    return fix_strand_closures(links), cg_charge


def fix_strand_closures(links):
#### MemPrO closes a periodic strand with a bond from its first NAG B1 to its last NAM B4, whose C1 and
#### O4 are both linked already; the closure that makes chemical sense is last NAM C1 -> first NAG O4
    used = {}
    for i, a1, j, a2, kind in links:
        if kind == 'glycosidic':
            used[(i, a1)] = used.get((i, a1), 0) + 1
            used[(j, a2)] = used.get((j, a2), 0) + 1
    fixed, out = 0, []
    for i, a1, j, a2, kind in links:
        if (kind == 'glycosidic' and used[(i, a1)] > 1 and used[(j, a2)] > 1
                and (j, 'C1') not in used and (i, 'O4') not in used):
            used[(i, a1)] -= 1
            used[(j, a2)] -= 1
            used[(j, 'C1')] = used[(i, 'O4')] = 1
            out.append((j, 'C1', i, 'O4', kind))
            fixed += 1
        else:
            out.append((i, a1, j, a2, kind))
    if fixed:
        print('Note: '+str(fixed)+' periodic strand closures given as first NAG B1 - last NAM B4 were built as '
              'last NAM C1 - first NAG O4 (the beta1,4 link that closes the strand)')
    return out


def links_from_distances():
#### links inferred from CG bead distances under the minimum image: each pair type is matched greedily,
#### shortest first, each bead taking at most one partner
    print('No CG topology given for the PG layer (-pgl): inferring its links from bead distances')
    res = g_var.pg_residues
    def candidates(donors, acceptors, cutoff):
        pairs = []
        for i, b1, a1 in donors:
            for j, b2, a2 in acceptors:
                if i != j:
                    d = bead_distance(res[i][1][b1]['coord'], res[j][1][b2]['coord'])
                    if d < cutoff:
                        pairs.append((d, i, a1, j, a2))
        return sorted(pairs)
    links, used = [], set()
    def take(pairs, kind):
        for d, i, a1, j, a2 in pairs:
            if (i, a1) not in used and (j, a2) not in used:
                used.update([(i, a1), (j, a2)])
                links.append((i, a1, j, a2, kind))
    of = lambda names, bead, atom: [(k, bead, atom) for k, (n, b) in enumerate(res) if n in names]
    glycosidic = sorted(candidates(of(['NAG'], 'B1', 'C1'), of(['NAM'], 'B4', 'O4'), LINK_CUTOFF['glycosidic']) +
                        candidates(of(['NAM'], 'B1', 'C1'), of(['NAG'], 'B4', 'O4'), LINK_CUTOFF['glycosidic']))
#### Insane writes each strand in order, so sequence neighbours are linked first; MemPrO's periodic closure
#### bond (first NAG B1 - last NAM B4) holds a wrong pair close together and would otherwise win
    take([pair for pair in glycosidic if pair[3] == pair[1] + 1], 'glycosidic')
    take(glycosidic, 'glycosidic')
#### strand closures: MemPrO joins a periodic strand's first NAG B1 to its last NAM B4 (which holds them
#### together); the free last NAM C1 and first NAG O4 of such a pair are linked (beads up to 9 A apart)
    closures = [pair for pair in candidates(of(['NAM'], 'B1', 'C1'), of(['NAG'], 'B4', 'O4'), 9.0)
                if bead_distance(res[pair[3]][1]['B1']['coord'], res[pair[1]][1]['B4']['coord']) < LINK_CUTOFF['glycosidic']]
    take(closures, 'glycosidic')
    take(candidates([(k, 'B5', 'CE3') for k, (n, b) in enumerate(res) if n == 'NAM'],
                    [(k, 'B1', 'N1') for k, (n, b) in enumerate(res) if n in STEMS], LINK_CUTOFF['stem']), 'stem')
    take(candidates([(k, DONOR[n][0], DONOR[n][1]) for k, (n, b) in enumerate(res) if n in DONOR],
                    [(k, 'B7', 'NZ3') for k, (n, b) in enumerate(res) if n in STEMS], LINK_CUTOFF['cross']), 'cross')
    return links, None


def residue_states(links):
    states = [{'c1': None, 'o4': None, 'donor': None, 'acceptor': None, 'stem_of': None} for r in g_var.pg_residues]
    for i, a1, j, a2, kind in links:
        if kind == 'glycosidic':
            states[i]['c1'], states[j]['o4'] = j, i
        elif kind == 'stem':
            states[j]['stem_of'] = i
        else:
            if states[j]['acceptor'] is not None or states[i]['donor'] is not None:
                raise TopologyError('PGL residue '+str(j+1)+' or '+str(i+1)+' takes part in two cross-links of the same kind')
            states[i]['donor'], states[j]['acceptor'] = j, i
    stems_of = {state['stem_of'] for state in states}
    for k, (resname, beads) in enumerate(g_var.pg_residues):
        if resname in STEMS and states[k]['stem_of'] is None:
            print('Warning: stem '+resname+' '+str(k+1)+' is not bonded to a NAM')
        if resname == 'NAM' and k not in stems_of:
            print('Warning: NAM '+str(k+1)+' carries no stem; its lactyl C=O is left without the amide N')
        if resname == 'UUPEP' and states[k]['donor'] is None:
            print('Warning: tripeptide '+str(k+1)+' is not a cross-link donor; it is built with a free mDAP L-carboxylate')
        if resname == 'SPEP' and states[k]['donor'] is not None:
            raise TopologyError('a pentapeptide (SPEP '+str(k+1)+') cannot be a cross-link donor')
    return states

######################################################################## templates

_base_cache = {}

def base_template(resname):
#### atoms (name, type, charge, mass, bead) and bonds of the free residue, from its fragment itp and pdb
    if resname not in _base_cache:
        location = gen.fragment_location(resname)
        frag_template, _ = at_mod.get_atomistic(location, resname)
        bead_of = {atom['atom']: bead for group in frag_template.values() for bead, atoms in group.items() for atom in atoms.values()}
        atoms, bonds, section, names = [], [], None, {}
        with open(location[:-4]+'.itp') as itp:
            for line in itp:
                line_sep = line.split(';')[0].split()
                if not line_sep or line_sep[0].startswith('#'):
                    continue
                if line_sep[0].startswith('['):
                    section = gen.strip_header(line)
                elif section == 'atoms':
                    names[int(line_sep[0])] = line_sep[4]
                    atoms.append({'name': line_sep[4], 'type': line_sep[1], 'q': float(line_sep[6]),
                                  'mass': float(line_sep[7]), 'bead': bead_of[line_sep[4]], 'index': int(line_sep[0])})
                elif section == 'bonds':
                    bonds.append((names[int(line_sep[0])], names[int(line_sep[1])]))
        _base_cache[resname] = (atoms, bonds)
    return _base_cache[resname]


def state_patches(resname, state):
    patches = []
    if resname in ['NAG', 'NAM']:
        patches += ['c1'] if state['c1'] is not None else []
        patches += ['o4'] if state['o4'] is not None else []
    else:
        patches += ['donor_'+resname] if state['donor'] is not None else []
        patches += ['acceptor'] if state['acceptor'] is not None else []
    return patches


def residue_template(resname, state):
#### atoms of a residue in its link state: [(fragment name, name, type, charge, mass, bead)], bonds
    atoms, bonds = base_template(resname)
    delete, rename, retype = set(), {}, {}
    for patch in state_patches(resname, state):
        delete.update(PATCH[patch]['delete'])
        rename.update(PATCH[patch]['rename'])
        retype.update(PATCH[patch]['set'])
    out = []
    for atom in atoms:
        if atom['name'] in delete:
            continue
        name = rename.get(atom['name'], atom['name'])
        atom_type, charge = retype.get(name, (atom['type'], atom['q']))
        out.append((atom['name'], name, atom_type, charge, atom['mass'], atom['bead']))
    kept = {atom[0]: atom[1] for atom in out}
    return out, [(kept[a], kept[b]) for a, b in bonds if a in kept and b in kept]

######################################################################## build

def build_pg_layer():
#### fits, links, writes and minimises the PG layer in <loc>/PGL; returns {'PGL': 1}
    print('\nConverting the peptidoglycan layer (PGL)')
    gen.mkdir_directory(g_var.working_dir+'PGL')
    if g_var.args.pgl is not None:
        links, cg_charge = links_from_topology(g_var.input_directory+gen.path_leaf(g_var.args.pgl)[1])
    else:
        links, cg_charge = links_from_distances()
    states = residue_states(links)
    counts = {kind: sum(link[4] == kind for link in links) for kind in ['glycosidic', 'stem', 'cross']}
    three_three = sum(1 for link in links if link[4] == 'cross' and g_var.pg_residues[link[0]][0] == 'UUPEP')
    print('PG links: '+str(counts['glycosidic'])+' glycosidic, '+str(counts['stem'])+' stem, '+
          str(counts['cross']-three_three)+' 3-4 and '+str(three_three)+' 3-3 cross-links')
    os.chdir(g_var.working_dir+'PGL')
    if not os.path.exists('PGL_merged.pdb'):
        g_var.pg_atoms, coords = fit_pg_layer(states)
        write_pg_itp(states, links, cg_charge)
        write_pg_cisres(states, links)
#### pdb residue names have four characters (UUPEP is written UUPE, its alternative name in the library)
        merge = [{'atom_name': atom['name'], 'residue_name': atom['resname'][:4], 'residue_id': atom['resid'] % 10000}
                 for atom in g_var.pg_atoms]
        coords, index_conversion = at_mod.index_conversion_generate(merge, coords)
        at_mod.write_pdb(merge, coords, index_conversion, g_var.working_dir+'PGL/PGL_all.pdb')
        minimise_pg_layer()
    else:
        g_var.pg_atoms = pg_atom_list(states)
    os.chdir(g_var.working_dir)
    return {'PGL': 1}


def pg_atom_list(states):
#### the PGL atoms in topology order, with the residue, link partners and fragment name of each
    pg_atoms = []
    for k, (resname, beads) in enumerate(g_var.pg_residues):
        atoms, bonds = residue_template(resname, states[k])
        for frag_name, name, atom_type, charge, mass, bead in atoms:
            pg_atoms.append({'name': name, 'frag_name': frag_name, 'type': atom_type, 'q': charge, 'mass': mass,
                             'bead': bead, 'resname': resname, 'resid': k + 1, 'res': k})
    g_var.pg_index = {(atom['res'], atom['name']): number for number, atom in enumerate(pg_atoms)}
    return pg_atoms


def fit_pg_layer(states):
    pg_atoms = pg_atom_list(states)
    coords = [None] * len(pg_atoms)
    templates = {}
    total = len(g_var.pg_residues)
    for k, (resname, beads) in enumerate(g_var.pg_residues):
        gen.print_progress('Converting PGL', k + 1, total)
        if resname not in templates:
            frag_template, frag_mass = at_mod.get_atomistic(gen.fragment_location(resname), resname)
            templates[resname] = (frag_template, frag_mass, at_mod.fragment_template(resname, frag_template))
        fitted = at_mod_np.fit_np_residue(resname, k, beads, *templates[resname], cg_resid=k)
        by_name = {atom['atom']: atom['coord'] for atom in fitted.values()}
        for number, atom in enumerate(pg_atoms):
            if atom['res'] == k:
                coords[number] = np.array(by_name[atom['frag_name']], dtype=float)
    gen.finish_progress('Converting PGL', total)
    return pg_atoms, coords

######################################################################## topology

def bond_graph(states, links):
    bonds = set()
    for k, (resname, beads) in enumerate(g_var.pg_residues):
        atoms, residue_bonds = residue_template(resname, states[k])
        for a, b in residue_bonds:
            bonds.add(tuple(sorted((g_var.pg_index[(k, a)], g_var.pg_index[(k, b)]))))
    for i, a1, j, a2, kind in links:
        bonds.add(tuple(sorted((g_var.pg_index[(i, a1)], g_var.pg_index[(j, a2)]))))
    neighbours = {number: set() for number in range(len(g_var.pg_atoms))}
    for a, b in bonds:
        neighbours[a].add(b)
        neighbours[b].add(a)
    return sorted(bonds), neighbours


def write_pg_itp(states, links, cg_charge):
    bonds, nb = bond_graph(states, links)
    atom_type = [atom['type'] for atom in g_var.pg_atoms]
    angles = sorted({(i, j, k) if i < k else (k, j, i) for j in nb for i in nb[j] for k in nb[j] if i != k})
    dihedrals = set()
    for j, k in bonds:
        for i in nb[j] - {k}:
            for l in nb[k] - {j}:
                if i != l:
                    dihedrals.add(min((i, j, k, l), (l, k, j, i)))
    dihedrals = sorted(dihedrals)
    excluded = set(bonds) | {(min(a[0], a[2]), max(a[0], a[2])) for a in angles}
    pairs = sorted({(min(d[0], d[3]), max(d[0], d[3])) for d in dihedrals} - excluded)
#### impropers as in the CHARMM36 residues: amide C (C, Cx, N, O), carboxylate (C, Cx, O2, O1), amide N (N, C, Cx, H)
    impropers = []
    for a in nb:
        n = sorted(nb[a])
        if atom_type[a] in ['C', 'CC2O1'] and len(n) == 3:
            o = [x for x in n if atom_type[x] in ['O', 'OC2D1']]
            nn = [x for x in n if atom_type[x] in ['NH1', 'NC2D1']]
            if o and nn:
                impropers.append((a, [x for x in n if x not in (o[0], nn[0])][0], nn[0], o[0]))
        elif atom_type[a] == 'CC' and len(n) == 3:
            ox = [x for x in n if atom_type[x] == 'OC']
            if len(ox) == 2:
                impropers.append((a, [x for x in n if x not in ox][0], ox[1], ox[0]))
        elif atom_type[a] in ['NH1', 'NC2D1'] and len(n) == 3:
            h = [x for x in n if atom_type[x] in ['H', 'HCP1']]
            c = [x for x in n if atom_type[x] in ['C', 'CC2O1']]
            if h and c:
                impropers.append((a, c[0], [x for x in n if x not in (h[0], c[0])][0], h[0]))
#### CMAP on every C-NH1-CT1/CTD1-C-NH1 backbone (L-Ala1, mDAP L-centre, D-Ala4)
    cmaps = []
    for n_atom in nb:
        if atom_type[n_atom] != 'NH1':
            continue
        for ca in nb[n_atom]:
            if atom_type[ca] not in ['CT1', 'CTD1']:
                continue
            for c_prev in nb[n_atom]:
                if atom_type[c_prev] != 'C':
                    continue
                for c in nb[ca]:
                    if atom_type[c] != 'C':
                        continue
                    for n_next in nb[c]:
                        if atom_type[n_next] == 'NH1' and n_next != n_atom:
                            cmaps.append((c_prev, n_atom, ca, c, n_next))
    charge = sum(atom['q'] for atom in g_var.pg_atoms)
    with open('PGL.itp', 'w') as itp:
        itp.write('; peptidoglycan layer written by CG2AT (at_mod_pg): CHARMM36 carbohydrate and protein types\n')
        itp.write('; '+str(len(g_var.pg_residues))+' residues, '+str(len(links))+' inter-residue links, net charge '+
                  str(round(charge, 3))+'\n; needs periodic-molecules = yes (the layer is bonded across the box)\n\n')
        itp.write('[ moleculetype ]\n; name  nrexcl\nPGL       3\n\n[ atoms ]\n; nr  type  resnr residu  atom  cgnr  charge  mass\n')
        for number, atom in enumerate(g_var.pg_atoms, 1):
            itp.write('{:7d} {:>8s} {:6d} {:>6s} {:>6s} {:7d} {:8.3f} {:9.4f}   ; {}\n'.format(
                number, atom['type'], atom['resid'], atom['resname'], atom['name'], number, atom['q'], atom['mass'], atom['bead']))
        itp.write('\n[ bonds ]\n'+''.join('{:7d}{:7d}  1\n'.format(a + 1, b + 1) for a, b in bonds))
        itp.write('\n[ pairs ]\n'+''.join('{:7d}{:7d}  1\n'.format(a + 1, b + 1) for a, b in pairs))
        itp.write('\n[ angles ]\n'+''.join('{:7d}{:7d}{:7d}  5\n'.format(a + 1, b + 1, c + 1) for a, b, c in angles))
        itp.write('\n[ dihedrals ]\n'+''.join('{:7d}{:7d}{:7d}{:7d}  9\n'.format(*[x + 1 for x in d]) for d in dihedrals))
        itp.write('\n[ dihedrals ]\n'+''.join('{:7d}{:7d}{:7d}{:7d}  2\n'.format(*[x + 1 for x in d]) for d in impropers))
        if cmaps:
            itp.write('\n[ cmap ]\n'+''.join('{:7d}{:7d}{:7d}{:7d}{:7d}  1\n'.format(*[x + 1 for x in c]) for c in cmaps))
    print('PGL topology: '+str(len(g_var.pg_atoms))+' atoms, net charge '+str(round(charge, 3))+', '+str(len(cmaps))+' CMAP terms')
    if cg_charge:
        report_charge(cg_charge)


def report_charge(cg_charge):
#### the all-atom charge of each residue against its CG charge (from the CG topology)
    aa_charge = {}
    for atom in g_var.pg_atoms:
        aa_charge[atom['res']] = aa_charge.get(atom['res'], 0) + atom['q']
#### sugars differ by +-0.18 at strand ends (CHARMM 14bb charges); only whole charges are reported
    differ = [k for k in aa_charge if abs(aa_charge[k] - cg_charge.get(k, 0)) > 0.5]
    total = sum(aa_charge.values()) - sum(cg_charge.values())
    if differ:
        kinds = {}
        for k in differ:
            key = g_var.pg_residues[k][0]+' (CG '+str(round(cg_charge[k], 2))+', all-atom '+str(round(aa_charge[k], 2))+')'
            kinds[key] = kinds.get(key, 0) + 1
        print('Warning: '+str(len(differ))+' PG residues carry a different charge all-atom than in the CG model: '+
              ', '.join(str(n)+' x '+key for key, n in kinds.items()))
        print('         the atomistic system is '+str(round(total, 2))+' e away from neutral; add counter-ions before production MD')
    g_var.pg_charge_offset = round(total, 3)

######################################################################## restraints and minimisation

def residue_index_map(states, k):
#### fragment atom number of residue k -> PGL index (0-based); an atom deleted by a link is mapped to the
#### atom that takes its place (C1's O1 -> partner O4, a donor's second carboxylate O -> acceptor NZ3)
    resname = g_var.pg_residues[k][0]
    atoms, bonds = base_template(resname)
    rename = {}
    for patch in state_patches(resname, states[k]):
        rename.update(PATCH[patch]['rename'])
    mapping = {}
    for atom in atoms:
        name = rename.get(atom['name'], atom['name'])
        if (k, name) in g_var.pg_index:
            mapping[atom['index']] = g_var.pg_index[(k, name)]
        elif atom['name'] == 'O1' and states[k]['c1'] is not None:
            mapping[atom['index']] = g_var.pg_index[(states[k]['c1'], 'O4')]
        elif resname in DONOR and atom['name'] == DONOR[resname][3] and states[k]['donor'] is not None:
            mapping[atom['index']] = g_var.pg_index[(states[k]['donor'], 'NZ3')]
    return mapping


def write_pg_cisres(states, links):
#### the CISRES restraints of each residue's fragment (stereocentres, chairs, amides) on the PGL indices,
#### plus the trans amides of the links (lactyl-Ala and cross-links)
    restraints, cache = [], {}
    for k, (resname, beads) in enumerate(g_var.pg_residues):
        if resname not in cache:
            cache[resname] = gro.cisres_restraints(resname, gen.fragment_location(resname)[:-4]+'.itp')
        mapping = residue_index_map(states, k)
        for ai, aj, ak, al, phi, dphi, kfac in cache[resname]:
            if all(x in mapping for x in (ai, aj, ak, al)):
                restraints.append([mapping[ai] + 1, mapping[aj] + 1, mapping[ak] + 1, mapping[al] + 1, phi, dphi, kfac])
    for i, a1, j, a2, kind in links:
        if kind == 'glycosidic':
            continue
        carbonyl_o = {'CE3': 'O31', 'C4': 'O4', 'C3': 'O3'}[a1]
        amide_h = 'HN1' if a2 == 'N1' else 'HZ3'
        quartet = [g_var.pg_index[(i, carbonyl_o)], g_var.pg_index[(i, a1)], g_var.pg_index[(j, a2)], g_var.pg_index[(j, amide_h)]]
        restraints.append([x + 1 for x in quartet] + [180.0, 15.0, 1000.0])
    gro.write_cisres_file('PGL_cisres.itp', restraints)


def minimise_pg_layer():
#### PGL on its own (periodic-molecules): a soft pass (stereocentre restraints at 1000) to relax the links
#### and the bonds between beads, then the full restraints, as for other non-protein residues
    gro.make_min('PGL')
    with open('topol_PGL.top', 'w') as topol:
        topol.write('; Include force field parameters\n#include \"../FINAL/'+g_var.forcefield+'/forcefield.itp\"\n')
        topol.write('#include \"PGL.itp\"\n'+gro.cisres_include('PGL')+'\n[ system ]\n; Name\npeptidoglycan\n\n[ molecules ]\n; Compound        #mols\nPGL    1\n')
    soft_pdb, converged = gro.em_pass('PGL', 'em_PGL_soft.mdp', '-DCISRES -DCHIRAL_K=1000.0', 'PGL_all.pdb', 'soft')
    success, converged = gro.strong_minimise('PGL', soft_pdb if soft_pdb else 'PGL_all.pdb', 'min')
    if not os.path.exists('PGL_merged.pdb'):
        raise TopologyError('minimisation of the peptidoglycan layer failed, see: '+g_var.working_dir+'PGL/MIN')
    report_geometry('PGL_merged.pdb', 'PGL after minimisation')

######################################################################## checks

def read_pg_coordinates(pdb_file):
#### PGL coordinates from a system pdb: the block of PG residues starting with the first PGL atom
    with open(pdb_file) as pdb_input:
        lines = [line for line in pdb_input if line.startswith('ATOM')]
    first = g_var.pg_atoms[0]
    start = next(k for k, line in enumerate(lines) if line[17:21].strip() == first['resname'][:4]
                 and line[12:16].strip() == first['name'])
    rows = lines[start:start + len(g_var.pg_atoms)]
    return np.array([[float(line[30:38]), float(line[38:46]), float(line[46:54])] for line in rows])


def report_geometry(pdb_file, label, start=None):
#### inverted stereocentres and rings out of their fragment chair, and the RMSD of the bead centres of
#### mass to the CG beads
    if start is None:
        with open(pdb_file) as pdb_input:
            X = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in pdb_input if l.startswith('ATOM')])
    else:
        X = start
    inverted, open_rings, rings, centres = [], 0, 0, 0
    by_res = {}
    for number, atom in enumerate(g_var.pg_atoms):
        by_res.setdefault(atom['res'], {})[atom['frag_name']] = number
    for k, (resname, beads) in enumerate(g_var.pg_residues):
        frag_template, _ = at_mod.get_atomistic(gen.fragment_location(resname), resname)
        coords, bonded = at_mod.fragment_template(resname, frag_template)
        names = {number: at_mod.frag_atom_name(frag_template, number) for number in coords}
        index = {number: by_res[k].get(names[number]) for number in coords}
        def unwrapped(numbers, centre):
            return [X[centre] + min_image(X[index[n]] - X[centre]) for n in numbers]
        for c, three, fourth, carried, sign in at_mod.stereo_template(resname, frag_template):
            if index[c] is None or any(index[x] is None for x in three):
                continue
            centres += 1
            v = np.array(unwrapped(three, index[c])) - X[index[c]]
            if np.sign(np.dot(v[0], np.cross(v[1], v[2]))) != sign:
                inverted.append(resname+str(k+1)+':'+names[c])
        for ring in at_mod.six_membered_rings(bonded):
            if any(index[x] is None for x in ring):
                continue
            theta_ref, amp_ref = at_mod.ring_pucker_theta([coords[x] for x in ring])
            theta, amplitude = at_mod.ring_pucker_theta(unwrapped(ring, index[ring[0]]))
            rings += 1
            if abs(theta - theta_ref) > 30:
                open_rings += 1
    rmsd = bead_rmsd(X)
    print(label+': '+str(len(inverted))+' of '+str(centres)+' stereocentres inverted'+
          (' ('+', '.join(inverted[:10])+(' ...' if len(inverted) > 10 else '')+')' if inverted else '')+
          ', '+str(open_rings)+' of '+str(rings)+' sugar rings out of their 4C1 chair, bead RMSD to CG '+str(np.round(rmsd, 2))+' A')
    return inverted, open_rings, rmsd


def bead_rmsd(X):
#### RMSD (A) of the all-atom bead centres of mass (fragment mapping) to the CG beads
    sums = {}
    for number, atom in enumerate(g_var.pg_atoms):
        key = (atom['res'], atom['bead'])
        ref = g_var.pg_residues[atom['res']][1][atom['bead']]['coord']
        x = ref + min_image(X[number] - ref)
        total, mass = sums.get(key, (np.zeros(3), 0.0))
        sums[key] = (total + x * atom['mass'], mass + atom['mass'])
    d2 = [np.sum((total / mass - g_var.pg_residues[k][1][bead]['coord'])**2) for (k, bead), (total, mass) in sums.items()]
    return float(np.sqrt(np.mean(d2)))
