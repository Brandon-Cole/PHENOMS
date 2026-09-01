from pathlib import Path

import pytest

from phenoms.io import load_trajectory
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


def test_detect_gromacs_tpr_preferred(tmp_path):
    rep = tmp_path / "rep_gmx"
    rep.mkdir()
    (rep / "traj.xtc").write_text("", encoding="utf-8")
    (rep / "topol.tpr").write_text("", encoding="utf-8")
    out = detect_replicate_input(rep)
    assert out.engine == "gromacs"
    assert out.topology_path.name == "topol.tpr"


def test_load_trajectory_with_tpr_topology():
    """
    .tpr topology is read via MDAnalysis (no GROMACS binary required) and handed
    to MDTraj as a converted PDB topology.
    """
    pytest.importorskip("MDAnalysis")
    datafiles = pytest.importorskip("MDAnalysisTests.datafiles")

    traj = load_trajectory(datafiles.XTC, top=datafiles.TPR)
    assert traj.n_frames > 0
    assert traj.n_atoms > 0
    assert traj.topology.n_bonds > 0


def test_load_trajectory_tpr_without_mdanalysis_raises(monkeypatch, tmp_path):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "MDAnalysis":
            raise ImportError("no MDAnalysis")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    tpr = tmp_path / "topol.tpr"
    tpr.write_text("", encoding="utf-8")
    xtc = tmp_path / "traj.xtc"
    xtc.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="phenoms\\[gromacs\\]"):
        load_trajectory(str(xtc), top=str(tpr))


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
