"""
Stereocentres between rigid fragment groups keep the fragment's handedness. CG beads carry no chirality,
so the beads around a junction can be in the mirror-image arrangement (Martini KLA: about half the
molecules at the Kdo II anomeric carbon, CR2); each group then fits its beads but the centre between
them comes out inverted. Here KLA at its fragment geometry has its Kdo II group turned about CR2 so CR2 is
inverted, and correct_group_chirality must put the group back where the fragment has it relative to
Kdo I, with the groups beyond it following. No GROMACS needed.

Run from the repository root with:  python -m unittest discover tests
"""

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = tempfile.mkdtemp(prefix='cg2at_chiral_')
RES = 'KLA'


def tearDownModule():
    shutil.rmtree(TMP, ignore_errors=True)


def handedness(xyz, centre, near):
    v = [xyz[x] - xyz[centre] for x in near]
    return np.sign(np.dot(v[0], np.cross(v[1], v[2])))


class TestGroupChirality(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # g_var reads the command line when it is first imported
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
        cls.groups = {group: [atom for bead in cls.frag[group].values() for atom in bead] for group in cls.frag}
        cls.xyz = {number: np.array(atom['coord'], dtype=float) for group in cls.frag.values()
                   for bead in group.values() for number, atom in bead.items()}
        cls.number = {atom['atom']: number for group in cls.frag.values() for bead in group.values()
                      for number, atom in bead.items()}
        cls.of_group = {atom: g for g, atoms in cls.groups.items() for atom in atoms}

    @classmethod
    def tearDownClass(cls):
        # drop the package so another test module importing it reads its own command line
        for name in [m for m in sys.modules if m == 'cg2at_lite' or m.startswith('cg2at_lite.')]:
            del sys.modules[name]
        os.chdir(ROOT)

    def test_inverted_junction_is_mended(self):
        at_mod, n = self.at_mod, self.number
        cr2, oq4, or6, cr3 = n['CR2'], n['OQ4'], n['OR6'], n['CR3']
        kdo2 = self.of_group[cr2]
        self.assertEqual(self.of_group[or6], kdo2)
        self.assertNotEqual(self.of_group[oq4], kdo2)
        family, order, checks, _ = at_mod.group_junctions(RES, self.frag, self.groups)
        self.assertIn(cr2, [c for c, near, sign in checks[kdo2]])
        beyond = family[kdo2][1]
        self.assertTrue(beyond)                       # the carboxylate and side chain of Kdo II
        # turn the Kdo II group (and the groups beyond it) 180 deg about the axis through CR2 that bisects
        # its two ring bonds: each group stays rigid, OQ4 (Kdo I) ends up on the other face of CR2
        a, b = self.xyz[or6] - self.xyz[cr2], self.xyz[cr3] - self.xyz[cr2]
        axis = a / np.linalg.norm(a) + b / np.linalg.norm(b)
        rot = at_mod.rotation_about_axis(axis, np.pi)
        residue = {x: {'coord': p.copy()} for x, p in self.xyz.items()}
        for g in [kdo2] + beyond:
            for x in self.groups[g]:
                residue[x]['coord'] = (self.xyz[x] - self.xyz[cr2]).dot(rot) + self.xyz[cr2]
        near = (oq4, or6, cr3)
        placed = {x: v['coord'] for x, v in residue.items()}
        self.assertNotEqual(handedness(placed, cr2, near), handedness(self.xyz, cr2, near))
        residue, mended = at_mod.correct_group_chirality(residue, RES, self.frag, self.groups)
        self.assertGreaterEqual(mended, 1)
        placed = {x: v['coord'] for x, v in residue.items()}
        self.assertEqual(handedness(placed, cr2, near), handedness(self.xyz, cr2, near))
        # Kdo II is back where the fragment has it relative to Kdo I, and the groups beyond came with it
        for g in [kdo2] + beyond:
            moved = np.array([placed[x] for x in self.groups[g]])
            frag = np.array([self.xyz[x] for x in self.groups[g]])
            self.assertLess(np.sqrt(np.mean(np.sum((moved - frag) ** 2, axis=1))), 0.05)

    def test_correct_placement_is_left_alone(self):
        residue = {x: {'coord': p.copy()} for x, p in self.xyz.items()}
        residue, mended = self.at_mod.correct_group_chirality(residue, RES, self.frag, self.groups)
        self.assertEqual(mended, 0)
        for x, p in self.xyz.items():
            np.testing.assert_allclose(residue[x]['coord'], p)


if __name__ == '__main__':
    unittest.main()
