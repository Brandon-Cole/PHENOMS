"""
Unit tests for phenoms.hbond and phenoms.io.
"""

import os
import pytest
import mdtraj as md

from phenoms.hbond import label_hbond
from phenoms.io import load_trajectory, load_and_select_residues


# Path to minimal fixture (2 residues, N and O each)
FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
MINIMAL_PDB = os.path.join(FIXTURES_DIR, "minimal.pdb")


class TestLabelHbond:
    """Test label_hbond with minimal topology: ALA1 N, ALA1 O, GLY2 N, GLY2 O (indices 0,1,2,3)."""

    @pytest.fixture
    def traj(self):
        return md.load(MINIMAL_PDB)

    def test_n_o_same_residue(self, traj):
        # Donor N (0), Acceptor O (1) -> ALA1 -- ALA1
        out = label_hbond((0, 0, 1), traj)
        assert out == "ALA1 -- ALA1"

    def test_n_o_across_residues(self, traj):
        # Donor N (0), Acceptor O (3) -> ALA1 -- GLY2
        out = label_hbond((0, 0, 3), traj)
        assert out == "ALA1 -- GLY2"

    def test_o_n_reverse(self, traj):
        # Donor O (1), Acceptor N (2) -> GLY2 -- ALA1 (N residue first in label)
        out = label_hbond((1, 0, 2), traj)
        assert out == "GLY2 -- ALA1"

    def test_non_n_o_returns_none(self, traj):
        # N-N would not be N-O; we don't have two N's as donor/acceptor in minimal, but O-O (1, 3) should return None
        out = label_hbond((1, 0, 3), traj)
        assert out is None


class TestIo:
    def test_load_trajectory(self):
        traj = load_trajectory(MINIMAL_PDB)
        assert traj.n_frames == 1
        assert traj.n_atoms == 4

    def test_load_and_select_residues_whole(self):
        traj = load_and_select_residues(MINIMAL_PDB, resid_range=None)
        assert traj.n_atoms == 4

    def test_load_and_select_residues_region(self):
        traj = load_and_select_residues(MINIMAL_PDB, resid_range=(1, 1))
        assert traj.n_atoms == 2  # ALA1 has N and O
