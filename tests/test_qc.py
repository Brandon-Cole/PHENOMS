"""Tests for phenoms.qc utilities."""

from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from phenoms.qc import (
    assess_series_convergence,
    check_mdp_key_consistency,
    parse_mdp,
    rmsd_convergence_report,
)


def test_parse_mdp_basic(tmp_path: Path):
    p = tmp_path / "a.mdp"
    p.write_text(
        "integrator = md ; comment\n"
        "dt = 0.002\n"
        "nsteps=250000\n",
        encoding="utf-8",
    )
    out = parse_mdp(p)
    assert out["integrator"] == "md"
    assert out["dt"] == "0.002"
    assert out["nsteps"] == "250000"


def test_check_mdp_consistency_mismatch(tmp_path: Path):
    a = tmp_path / "a.mdp"
    b = tmp_path / "b.mdp"
    a.write_text("integrator = md\ndt = 0.002\n", encoding="utf-8")
    b.write_text("integrator = md\ndt = 0.004\n", encoding="utf-8")
    rep = check_mdp_key_consistency([a, b], keys=("integrator", "dt"))
    assert rep["ok"] is False
    assert "dt" in rep["mismatches"]


def test_assess_series_convergence_true():
    x = np.concatenate([np.linspace(0.0, 0.2, 100), np.full(100, 0.2)])
    rep = assess_series_convergence(x, mean_tolerance=0.2, slope_tolerance=1e-3)
    assert rep["converged"] is True


def test_assess_series_convergence_false():
    x = np.linspace(0.0, 1.0, 200)
    rep = assess_series_convergence(x, mean_tolerance=0.05, slope_tolerance=1e-4)
    assert rep["converged"] is False


def test_rmsd_convergence_report_has_fields():
    # Minimal synthetic trajectory with 20 frames / 3 atoms.
    xyz = np.random.default_rng(0).normal(size=(20, 3, 3)).astype(np.float32) * 0.01
    top = md.Topology()
    chain = top.add_chain()
    res = top.add_residue("GLY", chain)
    top.add_atom("N", md.element.nitrogen, res)
    top.add_atom("CA", md.element.carbon, res)
    top.add_atom("C", md.element.carbon, res)
    traj = md.Trajectory(xyz=xyz, topology=top)
    rep = rmsd_convergence_report(traj)
    assert "converged" in rep
    assert rep["metric"] == "rmsd_to_first_frame_nm"
    assert rep["n_points"] == 20
