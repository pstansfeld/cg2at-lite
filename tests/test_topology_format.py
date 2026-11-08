"""
Topology files written by cg2at-lite stay readable above atom 9999: the restraint file's atom indices used to be
written in 5-wide columns with nothing between them, so a large PGL molecule (39,723 atoms in a 20 x 20 nm,
three-layer PG layer) gave "100041000310001" and grompp stopped. No GROMACS needed.

Run from the repository root with:  python -m unittest discover tests
"""

import os
import shutil
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = tempfile.mkdtemp(prefix='cg2at_top_')


def tearDownModule():
    shutil.rmtree(TMP, ignore_errors=True)


class TestTopologyFormat(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # g_var reads the command line when it is first imported
        os.chdir(TMP)
        sys.argv = ['cg2at', '-c', 'x.gro', '-fg', 'martini_3-0_charmm36']
        sys.path.insert(0, os.path.join(ROOT, 'src'))
        from cg2at_lite.bin import gro
        cls.gro = gro

    @classmethod
    def tearDownClass(cls):
        # drop the package so another test module importing it reads its own command line
        for name in [m for m in sys.modules if m == 'cg2at_lite' or m.startswith('cg2at_lite.')]:
            del sys.modules[name]
        os.chdir(ROOT)

    def test_restraint_indices_above_9999(self):
        restraints = [[10004, 10003, 10001, 9991, 3.69, 15.0, 1000.0],
                      [99993, 100013, 1000009, 9991, 50.88, 15.0, 'CHIRAL_K']]
        file_name = os.path.join(TMP, 'big_cisres.itp')
        self.gro.write_cisres_file(file_name, restraints)
        rows = []
        with open(file_name) as itp:
            for line in itp:
                fields = line.split(';')[0].split()
                if fields and fields[0].isdigit():
                    rows.append(fields)
        self.assertEqual(len(rows), 2)
        for row, restraint in zip(rows, restraints):
            self.assertEqual([int(x) for x in row[:4]], restraint[:4])
            self.assertEqual(row[4], '1')
            self.assertAlmostEqual(float(row[5]), restraint[4], 2)


if __name__ == '__main__':
    unittest.main()
