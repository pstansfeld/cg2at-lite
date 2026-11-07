"""
Fast tests of the peptidoglycan (PGL) back-mapping on a small relaxed MemPrO layer (tests/data/pg_layer):
links read from the CG topology against links inferred from bead distances, the charge of every residue
against the CG model, the stereocentres of the fitted fragments, a grompp check of the written PGL
topology, and the refusal of an unrelaxed layer. Nothing is minimised (a few seconds in all).

Needs GROMACS: gmx on the PATH, or the gmx executable in the CG2AT_GMX environment variable.
Run from the repository root with:  python -m unittest discover tests
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'tests', 'data', 'pg_layer')
GMX = shutil.which(os.environ.get('CG2AT_GMX', 'gmx'))
TMP = tempfile.mkdtemp(prefix='cg2at_pg_')


def tearDownModule():
    shutil.rmtree(TMP, ignore_errors=True)


def read_pgl_itp(file_name):
    """(atoms [(resnr, resname, name, charge)], bonds [(i, j)] 0-based) of a PGL itp."""
    atoms, bonds, section = [], [], None
    with open(file_name) as itp:
        for line in itp:
            fields = line.split(';')[0].split()
            if not fields:
                continue
            if fields[0] == '[':
                section = fields[1]
            elif section == 'atoms':
                atoms.append((int(fields[2]), fields[3], fields[4], float(fields[6])))
            elif section == 'bonds':
                bonds.append((int(fields[0]) - 1, int(fields[1]) - 1))
    return atoms, bonds


@unittest.skipIf(GMX is None, 'GROMACS (gmx) not found; set CG2AT_GMX')
class TestPeptidoglycanLayer(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # g_var reads the command line when it is first imported
        os.chdir(TMP)
        sys.argv = ['cg2at', '-c', os.path.join(DATA, 'pg_layer.gro'), '-pgl', os.path.join(DATA, 'PGL.itp'),
                    '-ff', 'charmm36-cg2at', '-fg', 'martini_3-0_charmm36', '-w', 'tip3p', '-loc', 'run', '-gmx', GMX]
        sys.path.insert(0, os.path.join(ROOT, 'src'))
        from cg2at_lite.bin import gen, gro, at_mod, at_mod_pg, read_in, g_var
        cls.g_var, cls.at_mod_pg = g_var, at_mod_pg
        gen.fetch_forcefield_water_info()
        gen.check_input_flag()
        gen.find_gromacs()
        gen.read_database_directories()
        gen.forcefield_selection()
        gen.fragment_selection()
        gen.check_water_molecules()
        gen.fetch_fragment_multi()
        gen.fetch_fragment_single()
        gro.collect_input()
        box = read_in.read_initial_cg_pdb(test=True)
        g_var.box_vec = box
        read_in.real_box_vectors(box)
        read_in.fix_pbc(box, box, np.array([0, 0, 0]))
        at_mod.sanity_check()
        at_mod_pg.collect_pg_residues()
        cls.links, cls.cg_charge = at_mod_pg.links_from_topology(g_var.input_directory + 'PGL.itp')
        cls.inferred, _ = at_mod_pg.links_from_distances()
        cls.states = at_mod_pg.residue_states(cls.links)
        os.makedirs(g_var.working_dir + 'PGL', exist_ok=True)
        os.chdir(g_var.working_dir + 'PGL')
        g_var.pg_atoms, cls.coords = at_mod_pg.fit_pg_layer(cls.states)
        at_mod_pg.write_pg_itp(cls.states, cls.links, cls.cg_charge)
        at_mod_pg.write_pg_cisres(cls.states, cls.links)
        cls.pgl_dir = g_var.working_dir + 'PGL/'
        os.chdir(TMP)

    def test_links_from_topology_and_inferred_agree(self):
        kinds = {kind: sum(link[4] == kind for link in self.links) for kind in ['glycosidic', 'stem', 'cross']}
        residues = [resname for resname, beads in self.g_var.pg_residues]
        self.assertEqual(kinds['stem'], sum(r in self.at_mod_pg.STEMS for r in residues))
        self.assertEqual(kinds['stem'], residues.count('NAM'))
        self.assertGreater(kinds['cross'], 0)
        # at most one glycosidic bond per sugar (one per sugar in a strand closed across the box)
        self.assertLessEqual(kinds['glycosidic'], residues.count('NAG') + residues.count('NAM'))
        self.assertGreater(kinds['glycosidic'], residues.count('NAG'))
        # the layer has a strand closed across the box (MemPrO 3b45ad5: last sugar C1 -> first sugar O4)
        glycosidic = [link for link in self.links if link[4] == 'glycosidic']
        ends = [(link[0], link[1]) for link in glycosidic] + [(link[2], link[3]) for link in glycosidic]
        self.assertEqual(len(ends), len(set(ends)))
        self.assertTrue(any(abs(link[0] - link[2]) > 1 for link in glycosidic))
        self.assertEqual(sorted(link[:4] for link in self.links), sorted(link[:4] for link in self.inferred))
        # both kinds of cross-link, and a tetrapeptide that only accepts a 3-3 link
        donors = {residues[link[0]] for link in self.links if link[4] == 'cross'}
        self.assertEqual(donors, {'UPEP', 'UUPEP'})
        self.assertTrue(any(residues[k] == 'UPEP' and state['donor'] is None and state['acceptor'] is not None
                            for k, state in enumerate(self.states)))

    def test_charges_match_the_cg_model(self):
        # MemPrO charges a cross-link acceptor's bead 7 and a non-donor tetrapeptide's bead 8, so every residue
        # carries its all-atom charge (sugars at the end of an open strand differ by the 14bb +-0.18 only)
        atoms, bonds = read_pgl_itp(self.pgl_dir + 'PGL.itp')
        aa = {}
        for resnr, resname, name, charge in atoms:
            aa[resnr - 1] = aa.get(resnr - 1, 0.0) + charge
        for k, charge in aa.items():
            self.assertLess(abs(charge - self.cg_charge[k]), 0.2, msg=self.g_var.pg_residues[k][0] + str(k + 1))
        self.assertEqual(self.g_var.pg_charge_offset, 0)
        self.assertAlmostEqual(sum(aa.values()), sum(self.cg_charge.values()), places=3)

    def test_charge_offset_reported_for_an_older_cg_model(self):
        # a CG topology from MemPrO before the non-donor tetrapeptide fix (bead 8 neutral) is 1 e more positive
        # per such tetrapeptide; the offset is reported
        free = [k for k, state in enumerate(self.states)
                if self.g_var.pg_residues[k][0] == 'UPEP' and state['donor'] is None]
        old = dict(self.cg_charge)
        for k in free:
            old[k] += 1
        self.at_mod_pg.report_charge(old)
        try:
            self.assertEqual(self.g_var.pg_charge_offset, -len(free))
        finally:
            self.at_mod_pg.report_charge(self.cg_charge)

    def test_fitted_fragments_keep_their_stereocentres(self):
        inverted, open_rings, rmsd = self.at_mod_pg.report_geometry(None, 'fitted PGL', np.array(self.coords))
        self.assertEqual(inverted, [])
        self.assertLess(rmsd, 2.0)

    def test_pgl_topology_passes_grompp(self):
        work = os.path.join(TMP, 'grompp')
        os.makedirs(work, exist_ok=True)
        for name in ['PGL.itp', 'PGL_cisres.itp']:
            shutil.copy(self.pgl_dir + name, work)
        ff = os.path.join(ROOT, 'src', 'cg2at_lite', 'forcefields', 'charmm36-cg2at.ff', 'forcefield.itp')
        with open(os.path.join(work, 'topol.top'), 'w') as top:
            top.write('#include "' + ff + '"\n#include "PGL.itp"\n#ifdef CISRES\n#include "PGL_cisres.itp"\n#endif\n'
                      '[ system ]\nPG\n\n[ molecules ]\nPGL 1\n')
        with open(os.path.join(work, 'em.mdp'), 'w') as mdp:
            mdp.write('define = -DCISRES\nintegrator = steep\nnsteps = 0\ncutoff-scheme = Verlet\nperiodic-molecules = yes\n'
                      'coulombtype = PME\nrcoulomb = 1.2\nrvdw = 1.2\n')
        box = self.g_var.box_vec.split()[1:4]
        with open(os.path.join(work, 'pgl.pdb'), 'w') as pdb:
            pdb.write('CRYST1{:9.3f}{:9.3f}{:9.3f}  90.00  90.00  90.00 P 1           1\n'.format(*[float(x) for x in box]))
            for number, (atom, xyz) in enumerate(zip(self.g_var.pg_atoms, self.coords), 1):
                pdb.write('ATOM  {:5d} {:4s} {:4s} {:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n'.format(
                    number % 100000, atom['name'], atom['resname'][:4], atom['resid'] % 10000, *xyz))
        out = subprocess.run([GMX, 'grompp', '-f', 'em.mdp', '-c', 'pgl.pdb', '-p', 'topol.top', '-o', 'em.tpr', '-maxwarn', '1'],
                             cwd=work, capture_output=True, text=True)
        log = out.stdout + out.stderr
        self.assertEqual(out.returncode, 0, log[-3000:])
        self.assertNotIn('No default', log)

    def test_unrelaxed_layer_refused(self):
        from cg2at_lite.bin.exceptions import InputError
        saved = self.g_var.pg_residues
        squeezed = []
        for resname, beads in saved:
            beads = {bead: dict(value) for bead, value in beads.items()}
            if resname in ['NAG', 'NAM']:
                b1 = beads['B1']['coord']
                beads['B3']['coord'] = b1 + 0.5 * self.at_mod_pg.min_image(beads['B3']['coord'] - b1)
            squeezed.append((resname, beads))
        self.g_var.pg_residues = squeezed
        try:
            with self.assertRaises(InputError):
                self.at_mod_pg.check_relaxed()
        finally:
            self.g_var.pg_residues = saved


if __name__ == '__main__':
    unittest.main()
