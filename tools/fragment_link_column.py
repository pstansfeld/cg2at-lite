#!/usr/bin/env python3
"""Check (or with --write, regenerate) the link column of CG2AT fragment PDBs.

The B-factor column (characters 61-66) of a fragment PDB marks the atoms that join beads: a number
is shared by the heavy atoms at the two ends of each bond (from the fragment's .itp) that crosses
from one bead to another, atoms joined through several crossings share one number, and every other
atom is 0. CG2AT-lite does not read the column (it takes the links from the .itp bonds); it
documents the mapping, and this keeps it in step when beads are remapped.

usage: python tools/fragment_link_column.py [FRAGMENT_DIR ...] [--write]
       (default: every non_protein folder of the bundled fragment libraries)
"""
import os, sys
from collections import defaultdict

def read_itp(path):
    names, bonds, sec, mols = [], [], None, 0
    for l in open(path):
        s = l.split(';')[0].strip()
        if not s or s.startswith('#'):
            continue
        if s.startswith('['):
            sec = s.strip('[] ').strip()
            mols += sec == 'moleculetype'
            continue
        if mols != 1:
            continue
        if sec == 'atoms':
            names.append(s.split()[4])
        elif sec == 'bonds':
            a, b = s.split()[:2]
            bonds.append((int(a), int(b)))
    return names, bonds

def analyse(frag_dir):
    res = os.path.basename(frag_dir.rstrip('/'))
    pdb, itp = os.path.join(frag_dir, res + '.pdb'), os.path.join(frag_dir, res + '.itp')
    if not os.path.exists(pdb) or not os.path.exists(itp):
        return res, 'skip (no pdb/itp)', None
    lines = open(pdb).readlines()
    names, bonds = read_itp(itp)
    unique = len(set(names)) == len(names)
    atoms, bead = [], None          # (line index, itp index, name, bead)
    for k, l in enumerate(lines):
        if l.startswith('['):
            bead = l.strip('[] \n')
        elif l.startswith('ATOM'):
            name, number = l[12:16].strip(), int(l[6:11])
            idx = names.index(name) + 1 if unique and name in names else (number if not unique else None)
            atoms.append((k, idx, name, bead))
    if len(atoms) < 2 or len({a[3] for a in atoms}) < 2:
        return res, 'skip (single bead)', None
    if any(a[1] is None or a[1] > len(names) for a in atoms):
        return res, 'skip (atoms do not map to itp)', None
    by_idx = {a[1]: a for a in atoms}
    heavy = lambda a: not a[2].startswith('H')
    cross = [(by_idx[i], by_idx[j]) for i, j in bonds if i in by_idx and j in by_idx
             and heavy(by_idx[i]) and heavy(by_idx[j]) and by_idx[i][3] != by_idx[j][3]]
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in cross:
        parent[find(a[1])] = find(b[1])
    label, nxt, new = {}, 1, {}
    for a in atoms:
        if a[1] in parent:
            r = find(a[1])
            if r not in label:
                label[r] = nxt
                nxt += 1
            new[a[0]] = float(label[r])
        else:
            new[a[0]] = 0.0
    def old(k):
        try:
            return float(lines[k][60:66])
        except ValueError:
            return None
    # compare the partition (which atoms share a number), not the numbers themselves
    def groups(values):
        g = defaultdict(set)
        for k, v in values.items():
            if v:
                g[v].add(k)
        return {frozenset(s) for s in g.values()}
    olds = {a[0]: old(a[0]) for a in atoms}
    short = any(len(lines[a[0]].rstrip('\n')) < 66 for a in atoms) or any(v is None for v in olds.values())
    status = 'ok' if not short and groups(olds) == groups(new) else ('fix (short lines)' if short else 'fix')
    return res, status + ' (%d crossings)' % len(cross), (lines, new)

def default_roots():
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'cg2at_lite', 'fragments')
    return sorted(os.path.join(base, lib, 'non_protein') for lib in os.listdir(base)
                  if os.path.isdir(os.path.join(base, lib, 'non_protein')))

def check_all(roots):
    """[(fragment folder, status)] for every fragment under roots."""
    out = []
    for root in roots:
        for res in sorted(os.listdir(root)):
            d = os.path.join(root, res)
            if os.path.isdir(d) and not res.startswith(('_', '.')):
                out.append((d, analyse(d)[1]))
    return out

def main(argv):
    write = '--write' in argv
    counts = defaultdict(int)
    for root in [a for a in argv if a != '--write'] or default_roots():
        for res in sorted(os.listdir(root)):
            d = os.path.join(root, res)
            if not os.path.isdir(d) or res.startswith(('_', '.')):
                continue
            name, status, data = analyse(d)
            counts[status.split(' (')[0]] += 1
            if status.startswith('fix') or status.startswith('skip (atoms'):
                print('%-55s %s' % (os.path.relpath(d, os.path.dirname(os.path.dirname(root))), status))
            if write and status.startswith('fix'):
                lines, new = data
                for k, v in new.items():
                    l = lines[k].rstrip('\n').ljust(66)
                    lines[k] = l[:60] + '{:6.2f}'.format(v) + l[66:] + '\n'
                with open(os.path.join(d, name + '.pdb'), 'w') as out:
                    out.writelines(lines)
    print(dict(counts))
    return 0 if write or not counts.get('fix') else 1

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
