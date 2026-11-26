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
from unittest.mock import patch

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

    def test_gromacs_command_uses_arguments_and_stdin(self):
        args, stdin = self.gro._split_gromacs_command(
            '/opt/GROMACS/gmx pdb2gmx -f chain.pdb -ter << EOF\n1\n0\nEOF'
        )
        self.assertEqual(args, ['/opt/GROMACS/gmx', 'pdb2gmx', '-f', 'chain.pdb', '-ter'])
        self.assertEqual(stdin, '1\n0')

        from cg2at_lite.bin import g_var
        g_var.args.v = 0
        g_var.gmx_version = True

        class Process:
            returncode = 0

            def communicate(self, input=None):
                self.input = input
                return b'', b''

        process = Process()
        with patch.object(self.gro.subprocess, 'Popen', return_value=process) as popen:
            _complete, error = self.gro.gromacs(
                ['/opt/GROMACS/gmx pdb2gmx -f chain.pdb -ter << EOF\n1\n0\nEOF',
                 os.path.join(TMP, 'missing.pdb')]
            )

        self.assertFalse(error)
        call = popen.call_args
        self.assertEqual(call.args[0], ['/opt/GROMACS/gmx', 'pdb2gmx', '-f', 'chain.pdb', '-ter'])
        self.assertNotIn('shell', call.kwargs)
        self.assertEqual(process.input, b'1\n0')

    def test_gromacs_mdrun_threads_disable_pin_and_set_environment(self):
        from cg2at_lite.bin import g_var
        g_var.args.v = 0
        g_var.gmx_version = True
        g_var.mdrun_ntomp = 3

        class Process:
            returncode = 0

            def communicate(self, input=None):
                return b'', b''

        with patch.object(self.gro.subprocess, 'Popen', return_value=Process()) as popen:
            _complete, error = self.gro.gromacs([
                '/opt/GROMACS/gmx mdrun -ntmpi 1 -pin on -deffnm test',
                os.path.join(TMP, 'missing.pdb'),
            ])

        self.assertFalse(error)
        call = popen.call_args
        self.assertEqual(
            call.args[0],
            ['/opt/GROMACS/gmx', 'mdrun', '-ntmpi', '1', '-pin', 'off', '-deffnm', 'test',
             '-ntomp', '3'],
        )
        self.assertEqual(call.kwargs['env']['OMP_NUM_THREADS'], '3')
        g_var.mdrun_ntomp = None

    def test_gromacs_nonzero_exit_is_reported_as_error(self):
        from cg2at_lite.bin import g_var
        g_var.args.v = 0
        g_var.gmx_version = True
        g_var.mdrun_ntomp = None

        class Process:
            returncode = 2

            def communicate(self, input=None):
                return b'', b'unrelated stderr'

        with patch.object(self.gro.subprocess, 'Popen', return_value=Process()):
            _complete, error = self.gro.gromacs([
                '/opt/GROMACS/gmx grompp -f test.mdp',
                os.path.join(TMP, 'missing.tpr'),
            ])

        self.assertTrue(error)

    def test_gromacs_diagnostic_and_version_streams_are_preserved(self):
        from cg2at_lite.bin import g_var
        g_var.args.v = 0
        g_var.mdrun_ntomp = None
        if hasattr(g_var, 'gmx_version'):
            del g_var.gmx_version

        class VersionProcess:
            returncode = 0

            def communicate(self, input=None):
                return b'GROMACS version: 2024.1\n', b''

        with patch.object(self.gro.subprocess, 'Popen', return_value=VersionProcess()):
            _complete, error = self.gro.gromacs([
                '/opt/GROMACS/gmx -version',
                os.path.join(TMP, 'version.txt'),
            ])
        self.assertFalse(error)
        self.assertTrue(g_var.gmx_version)

        class DiagnosticProcess:
            returncode = 0

            def communicate(self, input=None):
                return b'', b'Fatal error: invalid input'

        with patch.object(self.gro.subprocess, 'Popen', return_value=DiagnosticProcess()):
            _complete, error = self.gro.gromacs([
                '/opt/GROMACS/gmx grompp',
                os.path.join(TMP, 'missing.tpr'),
            ])
        self.assertTrue(error)


if __name__ == '__main__':
    unittest.main()
