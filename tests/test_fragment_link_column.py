"""Every non-protein fragment's link column (PDB B-factor) matches the bead-crossing bonds of its .itp.

Regenerate after remapping a fragment with: python tools/fragment_link_column.py --write
"""
import importlib.util
import os

TOOL = os.path.join(os.path.dirname(__file__), '..', 'tools', 'fragment_link_column.py')
spec = importlib.util.spec_from_file_location('fragment_link_column', TOOL)
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)


def test_link_columns_match_itp_bonds():
    stale = [(folder, status) for folder, status in tool.check_all(tool.default_roots()) if status.startswith('fix')]
    assert not stale, 'link column out of date (run tools/fragment_link_column.py --write): %s' % stale[:10]
