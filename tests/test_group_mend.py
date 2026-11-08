"""
Molecules trapped in the per-residue minimisation (a squashed or inverted centre, a ring held open) go back to
their own fit to the CG beads; nothing is replaced by a whole-fragment fit, and molecules that are merely
strained (an acyl chain turned, bonds stretched) keep their minimised coordinates. A chain bond placed through
a sugar ring of a neighbour is moved out before minimisation. KLA, no GROMACS needed.

Run from the repository root with:  python -m unittest discover tests
"""

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = tempfile.mkdtemp(prefix='cg2at_mend_')
RES = 'KLA'


def tearDownModule():
    shutil.rmtree(TMP, ignore_errors=True)


class TestGroupMend(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # g_var reads the command line when it is first imported
        os.chdir(TMP)
        sys.argv = ['cg2at', '-c', 'x.gro', '-fg', 'martini_3-0_charmm36']
        sys.path.insert(0, os.path.join(ROOT, 'src'))
        from cg2at_lite.bin import gen, at_mod, g_var
        cls.at_mod, cls.g_var = at_mod, g_var
        g_var.np_directories = [[os.path.join(ROOT, 'src', 'cg2at_lite', 'fragments', 'martini_3-0_charmm36',
                                              'non_protein') + '/']]
        location = gen.fragment_location(RES)
        gen.get_fragment_topology(RES, location, ring_groups=True)
        g_var.res_top[RES].setdefault('atom_masses', {})
        cls.frag, _ = at_mod.get_atomistic(location, RES)
        names = []
        with open(location[:-4] + '.itp') as itp:
            section, molecules = None, 0
            for line in itp:
                fields = line.split(';')[0].split()
                if not fields:
                    continue
                if fields[0] == '[':
                    section = fields[1]
                    molecules += section == 'moleculetype'
                elif section == 'atoms' and molecules == 1 and fields[0].isdigit():
                    names.append(fields[4])
        cls.names = names
        by_name = {atom['atom']: np.array(atom['coord'], dtype=float) / g_var.sf
                   for group in cls.frag.values() for bead in group.values() for atom in bead.values()}
        cls.xyz = np.array([by_name[n] for n in names])
        cls.groups = at_mod.fragment_groups(cls.frag)
        number = {atom['atom']: n for group in cls.frag.values() for bead in group.values() for n, atom in bead.items()}
        cls.group_names = {g: [at_mod.frag_atom_name(cls.frag, n) for n in atoms] for g, atoms in cls.groups.items()}
        coords, bonded = at_mod.fragment_template(RES, cls.frag)
        cls.rings = [[at_mod.frag_atom_name(cls.frag, n) for n in ring] for ring in at_mod.six_membered_rings(bonded)]
        cls.number = number

    @classmethod
    def tearDownClass(cls):
        for name in [m for m in sys.modules if m == 'cg2at_lite' or m.startswith('cg2at_lite.')]:
            del sys.modules[name]
        os.chdir(ROOT)

    def write(self, file_name, xyz, box=None):
        with open(file_name, 'w') as pdb:
            if box:
                pdb.write('CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n'.format(box, box, box, 90, 90, 90))
            for k, (name, x) in enumerate(zip(self.names * (len(xyz) // len(self.names)), xyz)):
                pdb.write('ATOM  {:5d} {:<4s} {:<4s} {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n'.format(
                    k % 99999 + 1, name, RES, k // len(self.names) + 1, *x))
        return file_name

    def read(self, file_name):
        with open(file_name) as pdb:
            return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in pdb if l.startswith('ATOM')])

    def test_only_trapped_molecule_goes_back_to_its_placement(self):
        # two copies, placed (the bead fit) and "minimised": copy 1 has a sugar ring carbon squashed onto its
        # neighbour (trapped); copy 2 has an acyl-chain group pulled out, its bonds stretched (not a trap)
        index = {n: k for k, n in enumerate(self.names)}
        n = len(self.names)
        placed = np.vstack([self.xyz, self.xyz + np.array([30.0, 0, 0])])
        minimised = placed + 0.05                                   # every atom moved a little by the minimisation
        ring = self.rings[0]
        a, b = index[ring[0]], index[ring[1]]
        minimised[a] = minimised[b] + 0.6 * (minimised[a] - minimised[b])       # squashed to 60%
        chain = next(g for g, members in self.group_names.items()
                     if not set(members) & {x for r in self.rings for x in r} and len(members) >= 6)
        rows = [n + index[x] for x in self.group_names[chain]]
        outward = minimised[rows].mean(axis=0) - minimised[n:].mean(axis=0)
        minimised[rows] += 0.4 * outward / np.linalg.norm(outward)          # pulled out: its bonds stretched
        put_back = self.at_mod.replace_trapped_molecules(RES, self.write('min.pdb', minimised, box=200.0),
                                                         self.write('placed.pdb', placed, box=200.0), 'min_out.pdb')
        self.assertEqual(put_back, 1)
        out = self.read('min_out.pdb')
        np.testing.assert_allclose(out[:n], placed[:n], atol=1e-3)           # back to its own bead fit
        np.testing.assert_allclose(out[n:], minimised[n:], atol=1e-3)        # untouched

    def test_threaded_ring_is_unthreaded(self):
        index = {n: k for k, n in enumerate(self.names)}
        ring = self.rings[0]
        r = [index[n] for n in ring]
        # a second copy placed so one of its chain bonds passes through the first copy's ring centre
        chain = next(g for g, members in self.group_names.items()
                     if not set(members) & {a for rr in self.rings for a in rr} and len(members) >= 6)
        heavy = [n for n in self.group_names[chain] if not n.startswith('H')]
        coords, bonded = self.at_mod.fragment_template(RES, self.frag)
        a = self.number[heavy[0]]
        b = next(x for x in bonded[a] if self.at_mod.frag_atom_name(self.frag, x) in heavy)
        A, B = self.xyz[index[heavy[0]]], self.xyz[index[self.at_mod.frag_atom_name(self.frag, b)]]
        P = self.xyz[r]
        centre = P.mean(axis=0)
        normal = np.linalg.svd(P - centre)[2][2]
        # turn copy 2 so its A->B bond lies along the ring normal, then put the bond's midpoint on the centre
        u = (B - A) / np.linalg.norm(B - A)
        v = np.cross(u, normal)
        if np.linalg.norm(v) < 1e-6:
            rot = np.eye(3)
        else:
            angle = np.arccos(np.clip(np.dot(u, normal), -1, 1))
            rot = self.at_mod.rotation_about_axis(v / np.linalg.norm(v), angle)
        second = (self.xyz - (A + B) / 2).dot(rot)
        if np.dot((second[index[self.at_mod.frag_atom_name(self.frag, b)]] - second[index[heavy[0]]]), normal) < 0:
            second = -second
        second += centre
        both = np.vstack([self.xyz, second])
        found = self.at_mod.unthread_rings(RES, self.write('threaded.pdb', both, box=200.0), 'threaded_out.pdb')
        self.assertGreaterEqual(found, 1)
        self.assertEqual(self.at_mod.unthread_rings(RES, 'threaded_out.pdb', 'threaded_again.pdb'), 0)


if __name__ == '__main__':
    unittest.main()
