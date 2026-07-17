from pathlib import Path

import pytest

from phenoms.preprocess import (
    detect_replicate_input,
    discover_replicate_dirs,
)


def test_detect_gromacs_requires_topology(tmp_path):
    rep = tmp_path / "rep1"
    rep.mkdir()
    (rep / "traj.xtc").write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="GROMACS replicate"):
        detect_replicate_input(rep)


def test_detect_openmm_pair(tmp_path):
    rep = tmp_path / "rep_openmm"
    rep.mkdir()
    (rep / "traj.dcd").write_text("", encoding="utf-8")
    (rep / "top.pdb").write_text("MODEL        1\nENDMDL\nEND\n", encoding="utf-8")
    out = detect_replicate_input(rep)
    assert out.engine == "openmm"
    assert out.trajectory_path.name == "traj.dcd"
    assert out.topology_path.name == "top.pdb"


def test_detect_amber_pair(tmp_path):
    rep = tmp_path / "rep_amber"
    rep.mkdir()
    (rep / "traj.nc").write_text("", encoding="utf-8")
    (rep / "sys.prmtop").write_text("", encoding="utf-8")
    out = detect_replicate_input(rep)
    assert out.engine == "amber"
    assert out.trajectory_path.name == "traj.nc"
    assert out.topology_path.name == "sys.prmtop"


def test_discover_single_replicate_dir(tmp_path):
    rep = tmp_path / "rep_single"
    rep.mkdir()
    (rep / "traj.nc").write_text("", encoding="utf-8")
    (rep / "sys.parm7").write_text("", encoding="utf-8")
    reps = discover_replicate_dirs(rep)
    assert reps == [rep.resolve()]


def test_discover_set_with_subdirs(tmp_path):
    root = tmp_path / "wt_set"
    rep_a = root / "rep_a"
    rep_b = root / "rep_b"
    rep_a.mkdir(parents=True)
    rep_b.mkdir(parents=True)
    (rep_a / "traj.dcd").write_text("", encoding="utf-8")
    (rep_a / "top.pdb").write_text("MODEL        1\nENDMDL\nEND\n", encoding="utf-8")
    (rep_b / "traj.nc").write_text("", encoding="utf-8")
    (rep_b / "sys.prmtop").write_text("", encoding="utf-8")
    reps = discover_replicate_dirs(root)
    assert reps == [rep_a.resolve(), rep_b.resolve()]
