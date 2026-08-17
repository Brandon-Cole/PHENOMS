"""
Regression tests for Rust H-bond kernel speed optimizations.

Ensures spatial-pruned kernels return the same bond hits as the brute-force
reference and that downstream pivot tables are unchanged.
"""

from __future__ import annotations

import os
import time
from pathlib import Path

import mdtraj as md
import numpy as np
import pytest

from phenoms.analysis import create_pivot_table
from phenoms.hbond import (
    _process_frames_rust,
    _run_baker_hubbard_rs,
    _run_baker_hubbard_rs_brute,
    _topology_donor_acceptor_lists,
    process_frames,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
TINY_PROTEIN = FIXTURES_DIR / "tiny_protein.pdb"
_REPO_ROOT = FIXTURES_DIR.parent.parent
KRAS_PDB = _REPO_ROOT / "H_bond_work" / "rep1_md_Kras_only_500ns_none_nojump_center.pdb"


def _require_rust():
    if _run_baker_hubbard_rs is None or _run_baker_hubbard_rs_brute is None:
        pytest.skip("phenoms_hbond_rs not built (maturin develop --release)")


def _bond_rows(raw) -> np.ndarray:
    if raw is None or raw.size == 0:
        return np.zeros((0, 4), dtype=np.uint32)
    rows = np.asarray(raw, dtype=np.uint32)
    return rows[np.lexsort((rows[:, 3], rows[:, 2], rows[:, 1], rows[:, 0]))]


@pytest.fixture(scope="module")
def multi_frame_traj():
    if not TINY_PROTEIN.is_file():
        pytest.skip(f"fixture missing: {TINY_PROTEIN}")
    base = md.load(str(TINY_PROTEIN))
    n_frames = 40
    xyz = np.tile(base.xyz, (n_frames, 1, 1)).astype(np.float64)
    rng = np.random.default_rng(42)
    xyz += rng.normal(0.0, 0.002, size=xyz.shape)
    return md.Trajectory(xyz=xyz.astype(np.float32), topology=base.topology)


def test_default_kernel_matches_brute_on_small_search_space(multi_frame_traj):
    _require_rust()
    traj = multi_frame_traj
    n_frames = traj.n_frames
    n_atoms = traj.n_atoms
    xyz_flat = np.ascontiguousarray(traj.xyz.reshape(-1), dtype=np.float64)
    pairs, acceptors = _topology_donor_acceptor_lists(traj, backbone_only=True)

    brute = _run_baker_hubbard_rs_brute(
        xyz_flat, n_frames, n_atoms, pairs, acceptors
    )
    adaptive = _run_baker_hubbard_rs(
        xyz_flat, n_frames, n_atoms, pairs, acceptors
    )
    assert np.array_equal(_bond_rows(brute), _bond_rows(adaptive))


def test_spatial_kernel_matches_brute(multi_frame_traj):
    _require_rust()
    traj = multi_frame_traj
    n_frames = traj.n_frames
    n_atoms = traj.n_atoms
    xyz_flat = np.ascontiguousarray(traj.xyz.reshape(-1), dtype=np.float64)
    pairs, acceptors = _topology_donor_acceptor_lists(traj, backbone_only=True)
    assert pairs and acceptors

    brute = _run_baker_hubbard_rs_brute(
        xyz_flat, n_frames, n_atoms, pairs, acceptors
    )
    spatial = _run_baker_hubbard_rs(
        xyz_flat, n_frames, n_atoms, pairs, acceptors
    )
    assert np.array_equal(_bond_rows(brute), _bond_rows(spatial))


def test_spatial_kernel_matches_brute_all_bonds(multi_frame_traj):
    _require_rust()
    traj = multi_frame_traj
    n_frames = traj.n_frames
    n_atoms = traj.n_atoms
    xyz_flat = np.ascontiguousarray(traj.xyz.reshape(-1), dtype=np.float64)
    pairs, acceptors = _topology_donor_acceptor_lists(traj, backbone_only=False)
    assert pairs and acceptors

    brute = _run_baker_hubbard_rs_brute(
        xyz_flat, n_frames, n_atoms, pairs, acceptors
    )
    spatial = _run_baker_hubbard_rs(
        xyz_flat, n_frames, n_atoms, pairs, acceptors
    )
    assert np.array_equal(_bond_rows(brute), _bond_rows(spatial))


def test_process_frames_rust_matches_mdtraj_fallback(multi_frame_traj):
    _require_rust()
    traj = multi_frame_traj
    rust_df = _process_frames_rust(
        traj, sub_frames=traj.n_frames, backbone_only=True, n_jobs=4
    )
    md_df = process_frames(
        traj, sub_frames=traj.n_frames, n_jobs=1, use_rust=False
    )
    assert rust_df is not None
    rust_keys = rust_df[
        ["Frame", "Donor Index", "Hydrogen Index", "Acceptor Index"]
    ].astype(int)
    md_keys = md_df[
        ["Frame", "Donor Index", "Hydrogen Index", "Acceptor Index"]
    ].astype(int)
    assert rust_keys.equals(md_keys.reset_index(drop=True))


def test_create_pivot_table_matches_legacy_implementation(multi_frame_traj):
    _require_rust()
    traj = multi_frame_traj
    hbonds_df = _process_frames_rust(
        traj, sub_frames=traj.n_frames, backbone_only=True, n_jobs=1
    )
    assert hbonds_df is not None and not hbonds_df.empty
    labels = sorted(hbonds_df["Bond Label"].unique())
    new_pivot = create_pivot_table(hbonds_df, labels)
    legacy = hbonds_df.pivot_table(
        index="Bond Label", columns="Frame", aggfunc="size", fill_value=0
    ).reindex(labels, fill_value=0)
    assert new_pivot.equals(legacy)


@pytest.mark.skip(reason="optional local benchmark; run manually when profiling")
def test_spatial_kernel_faster_than_brute_on_larger_traj():
    _require_rust()
    if not KRAS_PDB.is_file():
        pytest.skip(f"Kras trajectory not found: {KRAS_PDB}")
    traj = md.load(str(KRAS_PDB))
    traj = traj[:120]
    n_frames = traj.n_frames
    n_atoms = traj.n_atoms
    xyz_flat = np.ascontiguousarray(traj.xyz.reshape(-1), dtype=np.float64)
    pairs, acceptors = _topology_donor_acceptor_lists(traj, backbone_only=True)
    assert pairs and acceptors

    t0 = time.perf_counter()
    _run_baker_hubbard_rs_brute(xyz_flat, n_frames, n_atoms, pairs, acceptors)
    brute_sec = time.perf_counter() - t0

    t0 = time.perf_counter()
    _run_baker_hubbard_rs(xyz_flat, n_frames, n_atoms, pairs, acceptors)
    spatial_sec = time.perf_counter() - t0

    # Spatial should not be slower; allow small timing noise on tiny runs.
    assert spatial_sec <= brute_sec * 1.05, (
        f"spatial ({spatial_sec:.4f}s) slower than brute ({brute_sec:.4f}s)"
    )
