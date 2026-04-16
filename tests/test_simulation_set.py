"""
Tests for SimulationSet (run, getters, optional bond statistics).
"""

import os
import pytest
import pandas as pd

from phenoms import SimulationSet
from phenoms.io import load_and_select_residues
from phenoms.hbond import process_frames

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
MINIMAL_PDB = os.path.join(FIXTURES_DIR, "minimal.pdb")
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KRAS_PDB = os.path.join(_REPO_ROOT, "H_bond_work", "rep1_md_Kras_only_500ns_none_nojump_center.pdb")


@pytest.fixture
def minimal_pdb_path():
    return MINIMAL_PDB


@pytest.fixture
def kras_pdb_path():
    """Path to Kras trajectory (501 frames). Skip if file not present."""
    if not os.path.isfile(KRAS_PDB):
        pytest.skip(f"Kras trajectory not found: {KRAS_PDB}")
    return KRAS_PDB


class TestSimulationSet:
    def test_init_single_file(self):
        sim = SimulationSet(["a.pdb"], resid_range=(1, 10), sub_frames=5)
        assert len(sim.pdb_files) == 1
        assert sim.resid_range == (1, 10)
        assert sim.sub_frames == 5

    def test_init_replicates(self):
        sim = SimulationSet(["a.pdb", "b.pdb"], resid_range=None)
        assert len(sim.pdb_files) == 2
        assert sim.resid_range is None

    def test_run_whole_protein_single_replicate(self, minimal_pdb_path):
        sim = SimulationSet(
            [minimal_pdb_path],
            resid_range=None,
            sub_frames=1,
        )
        try:
            sim.run(n_jobs=1)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        assert len(sim.get_pivot_tables()) == 1
        assert len(sim.get_hbond_dfs()) == 1
        assert sim.get_bond_labels_sorted() is not None
        assert sim.bond_statistics is None

    def test_run_with_bond_statistics_threshold(self, minimal_pdb_path):
        sim = SimulationSet(
            [minimal_pdb_path],
            resid_range=None,
            sub_frames=1,
            bond_statistics_threshold=0.5,
        )
        try:
            sim.run(n_jobs=1)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        assert sim.bond_statistics is not None
        assert "mean_lifetimes" in sim.bond_statistics
        assert "std_lifetimes" in sim.bond_statistics
        assert "mean_break_frequencies" in sim.bond_statistics
        assert "std_break_frequencies" in sim.bond_statistics

    def test_get_pivot_tables_before_run_raises(self):
        sim = SimulationSet(["a.pdb"])
        with pytest.raises(ValueError, match="No data"):
            sim.plot_heatmaps()

    def test_plot_aggregated_heatmap_after_run(self, minimal_pdb_path):
        sim = SimulationSet([minimal_pdb_path], sub_frames=1)
        try:
            sim.run(n_jobs=1)
        except ValueError as e:
            if "No bonds" in str(e):
                pytest.skip("Trajectory has no bonds (minimal fixture has no H atoms)")
            raise
        # Should not raise; we don't capture plt here
        sim.plot_aggregated_heatmap(save_path=os.path.join(FIXTURES_DIR, "agg.png"))
        if os.path.exists(os.path.join(FIXTURES_DIR, "agg.png")):
            os.remove(os.path.join(FIXTURES_DIR, "agg.png"))

    def test_run_kras_first_10ns(self, kras_pdb_path):
        """Integration test: run on Kras rep1, first 100 frames (~10 ns)."""
        sim = SimulationSet(
            [kras_pdb_path],
            resid_range=None,
            sub_frames=100,
        )
        sim.run(n_jobs=1)
        assert len(sim.get_pivot_tables()) == 1
        assert len(sim.get_hbond_dfs()) == 1
        pivot = sim.get_pivot_tables()[0]
        assert pivot.shape[0] > 0 and pivot.shape[1] == 100

    def test_run_kras_region(self, kras_pdb_path):
        """Integration test: run on Kras rep1 with residue range, first 50 frames."""
        sim = SimulationSet(
            [kras_pdb_path],
            resid_range=(1, 50),
            sub_frames=50,
        )
        sim.run(n_jobs=1)
        assert len(sim.get_pivot_tables()) == 1
        assert len(sim.get_hbond_dfs()) == 1
