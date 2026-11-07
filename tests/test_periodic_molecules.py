"""
Per-molecule checks on molecules split across the periodic boundary: a POPE built from its fragment and
wrapped into a small box is not taken for distorted (its bonds across the boundary used to measure a box
length, so every lipid crossing it was replaced by a rigid fragment), it is written back unchanged in its
own periodic images, and a cis double bond turned trans is counted. No GROMACS needed.

Run from the repository root with:  python -m unittest discover tests
"""

import math
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = tempfile.mkdtemp(prefix='cg2at_pbc_')
RES = 'POPE'
BOX = 50.0


def tearDownModule():
    shutil.rmtree(TMP, ignore_errors=True)


class TestPeriodicMolecules(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # g_var reads the command line when it is first imported
        os.chdir(TMP)
        sys.argv = ['cg2at', '-c', 'x.gro', '-fg', 'martini_3-0_charmm36']
        sys.path.insert(0, os.path.join(ROOT, 'src'))
        from cg2at_lite.bin import gen, at_mod, g_var
        cls.gen, cls.at_mod, cls.g_var = gen, at_mod, g_var
        g_var.np_directories = [[os.path.join(ROOT, 'src', 'cg2at_lite', 'fragments', 'martini_3-0_charmm36',
                                              'non_protein') + '/']]
        location = gen.fragment_location(RES)
        gen.get_fragment_topology(RES, location)
        g_var.res_top[RES].setdefault('atom_masses', {})
        cls.frag, _ = at_mod.get_atomistic(location, RES)
        # the molecule in itp atom order, at the fragment's own geometry
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
        by_name = {atom['atom']: np.array(atom['coord'], dtype=float) / g_var.sf
                   for group in cls.frag.values() for bead in group.values() for atom in bead.values()}
        cls.names = names
        cls.xyz = np.array([by_name[n] for n in names])
        # centred on a box corner, so the molecule is split across all three boundaries
        cls.xyz -= cls.xyz.mean(axis=0)

    @classmethod
    def tearDownClass(cls):
        # g_var keeps the command line it read on import: drop the package so another test module
        # importing it reads its own
        for name in [m for m in sys.modules if m == 'cg2at_lite' or m.startswith('cg2at_lite.')]:
            del sys.modules[name]
        os.chdir(ROOT)

    def write(self, file_name, xyz, box=True):
        with open(file_name, 'w') as pdb:
            if box:
                pdb.write('CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f} P 1           1\n'.format(BOX, BOX, BOX, 90, 90, 90))
            for k, (name, x) in enumerate(zip(self.names, xyz)):
                pdb.write('ATOM  {:5d} {:<4s} {:<4s} {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n'.format(
                    k + 1, name, RES, 1, *x))
        return file_name

    def wrapped(self):
        return np.mod(self.xyz, BOX)

    def test_split_molecule_is_not_distorted(self):
        split = self.write('split.pdb', self.wrapped())
        self.assertGreater(np.ptp(self.wrapped(), axis=0).max(), BOX / 2)  # it is split
        self.assertEqual(self.at_mod.replace_distorted_molecules(RES, split, 'split_out.pdb'), 0)
        with open(split) as a, open('split_out.pdb') as b:
            self.assertEqual(a.read(), b.read())

    def test_split_molecule_without_box_is_flagged(self):
        # without a CRYST1 record nothing is made whole: the split molecule is still caught as broken
        split = self.write('split_nobox.pdb', self.wrapped(), box=False)
        self.assertEqual(self.at_mod.replace_distorted_molecules(RES, split, 'split_nobox_out.pdb'), 1)

    def test_trans_double_bond_is_counted(self):
        at_mod, frag = self.at_mod, self.frag
        self.assertEqual(at_mod.trans_double_bonds(RES, self.write('cis.pdb', self.wrapped())), (0, 1, 0))
        # turn C29=C210 trans: the C210 side rotated 180 degrees about the bond
        coords, bonded = at_mod.fragment_template(RES, frag)
        number = {at_mod.frag_atom_name(frag, n): n for n in coords}
        a, b = number['C29'], number['C210']
        side, stack = {b}, [b]
        while stack:
            for n in bonded[stack.pop()]:
                if n != a and n not in side:
                    side.add(n)
                    stack.append(n)
        side |= {h for h, parent in at_mod.hydrogen_parents(RES, frag).items() if parent in side}
        side_names = {at_mod.frag_atom_name(frag, n) for n in side}
        xyz = self.xyz.copy()
        A, B = xyz[self.names.index('C29')], xyz[self.names.index('C210')]
        k = (B - A) / np.linalg.norm(B - A)
        for i, name in enumerate(self.names):
            if name in side_names:
                v = xyz[i] - A
                xyz[i] = A + 2 * np.dot(v, k) * k - v
        self.assertEqual(at_mod.trans_double_bonds(RES, self.write('trans.pdb', np.mod(xyz, BOX))), (1, 1, 1))

    def test_box_vectors(self):
        box = self.at_mod.pdb_box(['CRYST1   50.000   60.000   70.000  90.00  90.00 120.00 P 1           1\n'])
        self.assertAlmostEqual(box[1][0], 60 * math.cos(math.radians(120)), 3)
        d = self.at_mod.minimum_image(np.array([49.0, 0.5, 69.0]), box)
        self.assertLess(np.linalg.norm(d), 2.0)
        self.assertIsNone(self.at_mod.pdb_box(['ATOM      1 ...\n']))


if __name__ == '__main__':
    unittest.main()
