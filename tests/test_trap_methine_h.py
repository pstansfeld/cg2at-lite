"""
A stereocentre whose methine hydrogen has been pushed through to the face of its three heavy neighbours is
trapped, although the handedness of the heavy atoms alone still reads right (KLA CB2 in a packed LPS leaflet).
No GROMACS needed. Run from the repository root with:  python -m unittest discover tests
"""

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = tempfile.mkdtemp(prefix='cg2at_methine_')
RES = 'KLA'


def tearDownModule():
    shutil.rmtree(TMP, ignore_errors=True)


class TestMethineHydrogen(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.chdir(TMP)
        sys.argv = ['cg2at', '-c', 'x.gro', '-fg', 'martini_3-0_charmm36']
        sys.path.insert(0, os.path.join(ROOT, 'src'))
        from cg2at_lite.bin import gen, at_mod, g_var
        cls.at_mod = at_mod
        g_var.np_directories = [[os.path.join(ROOT, 'src', 'cg2at_lite', 'fragments', 'martini_3-0_charmm36',
                                              'non_protein') + '/']]
        location = gen.fragment_location(RES)
        gen.get_fragment_topology(RES, location, ring_groups=True)
        g_var.res_top[RES].setdefault('atom_masses', {})
        cls.frag, _ = at_mod.get_atomistic(location, RES)
        cls.molecule = {n: {'coord': np.array(atom['coord'], dtype=float) / g_var.sf}
                        for group in cls.frag.values() for bead in group.values() for n, atom in bead.items()}

    def centre(self, name):
        for entry in self.at_mod.stereo_template(RES, self.frag):
            if self.at_mod.frag_atom_name(self.frag, entry[0]) == name:
                return entry
        self.fail(name + ' is not a stereocentre')

    def test_fragment_is_not_trapped(self):
        self.assertFalse(self.at_mod.trapped(self.molecule, RES, self.frag))
        self.assertEqual(self.at_mod.trap_sites(self.molecule, RES, self.frag), set())

    def test_hydrogen_through_is_trapped(self):
        c, three, fourth, carried, sign = self.centre('CB2')
        molecule = {n: {'coord': atom['coord'].copy()} for n, atom in self.molecule.items()}
        a, b, d = [molecule[x]['coord'] for x in three]
        normal = np.cross(b - a, d - a)
        normal /= np.linalg.norm(normal)
        h = molecule[fourth]['coord']
        molecule[fourth]['coord'] = h - 2 * np.dot(h - a, normal) * normal   # reflect through the heavy-atom plane
        v = [molecule[x]['coord'] - molecule[c]['coord'] for x in three]
        self.assertEqual(np.sign(np.dot(v[0], np.cross(v[1], v[2]))), sign)   # heavy atoms still read right
        self.assertTrue(self.at_mod.trapped(molecule, RES, self.frag))
        self.assertIn(c, self.at_mod.trap_sites(molecule, RES, self.frag))


if __name__ == '__main__':
    unittest.main()
