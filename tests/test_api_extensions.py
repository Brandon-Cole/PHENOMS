"""
Tests for SimulationSet API extensions (backbone_only, native trajectories).
"""

import os

import pytest

from phenoms import SimulationSet, ComparisonSet
from phenoms.io import normalize_topology_list

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
MINIMAL_PDB = os.path.join(FIXTURES_DIR, "minimal.pdb")


def test_backbone_only_default():
    sim = SimulationSet([MINIMAL_PDB], sub_frames=1)
    assert sim.backbone_only is True
    assert sim.input_kind == "pdb"


def test_all_bonds_mode_flag():
    sim = SimulationSet([MINIMAL_PDB], sub_frames=1, backbone_only=False)
    assert sim.backbone_only is False


def test_from_trajectories_requires_topology():
    with pytest.raises(ValueError, match="topology"):
        SimulationSet(trajectories=[MINIMAL_PDB])


def test_from_trajectories_shared_topology():
    sim = SimulationSet.from_trajectories(
        [MINIMAL_PDB],
        topology=MINIMAL_PDB,
        sub_frames=1,
        backbone_only=True,
    )
    assert sim.input_kind == "trajectory"
    assert sim.topologies == [MINIMAL_PDB]
    assert sim.pdb_files == [MINIMAL_PDB]


def test_mutually_exclusive_pdb_and_traj():
    with pytest.raises(ValueError, match="either pdb_files"):
        SimulationSet(pdb_files=[MINIMAL_PDB], trajectories=[MINIMAL_PDB], topology=MINIMAL_PDB)


def test_normalize_topology_list_broadcast():
    tops = normalize_topology_list(3, topology="/tmp/top.pdb")
    assert tops == ["/tmp/top.pdb"] * 3


def test_comparison_set_allows_construct_before_run():
    a = SimulationSet([MINIMAL_PDB], sub_frames=1)
    b = SimulationSet([MINIMAL_PDB], sub_frames=1)
    comp = ComparisonSet(a, b, label_a="a", label_b="b")
    with pytest.raises(ValueError, match="must be run first"):
        comp.compare()


def test_run_all_bonds_mode():
    sim = SimulationSet([MINIMAL_PDB], sub_frames=1, backbone_only=False)
    try:
        sim.run(n_jobs=1, use_rust=False)
    except ValueError as e:
        if "No bonds" in str(e) or "empty" in str(e).lower():
            pytest.skip("Minimal fixture cannot form H-bonds")
        raise
    assert sim.get_occupancy_tables() is not None
