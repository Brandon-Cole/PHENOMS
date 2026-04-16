#!/usr/bin/env python3
"""
Official kernel benchmark: Rust vs MDTraj vs MDAnalysis.

Default CSV output: under :func:`phenoms.outputs.default_output_root` in
``benchmarks/kernel/`` (override with ``--out-dir``). See also ``docker/README.md``
for a reproducible **1 CPU** Docker run.

Three-way timing on the **same** first-``N`` frames for each replicate PDB:

1. **Rust** — ``run_baker_hubbard`` (single-thread): raw (frame, donor, H, acceptor) hits,
   same geometry as MDTraj Baker–Hubbard (0.25 nm H···A, 120° at H).

2. **MDTraj** — **docs-style** single call: ``md.baker_hubbard(traj_n, periodic=False)``
   on ``traj_n = traj[:N]`` (defaults: ``freq=0.1``, ``exclude_water=True``, …).
   Returned rows are **persistent** (d,h,a) triplets, not the same statistic as Rust.

3. **MDAnalysis** — **docs-style** ``HydrogenBondAnalysis``: PDBs without bond/charge
   data use explicit ``donors_sel`` / ``hydrogens_sel`` / ``acceptors_sel`` and
   ``d_h_cutoff`` per MDAnalysis guidance. Criteria differ from Baker–Hubbard
   (e.g. default HBA angle / D–A cutoff); timings are **API representative**, not
   matched-geometry kernels.

Threading: before importing NumPy / MDAnalysis / MDTraj, this script sets standard
OpenMP/BLAS thread caps to **4** (via ``os.environ.setdefault``). Override by
exporting variables **before** launching the script.

Rust remains explicitly single-threaded (``run_baker_hubbard``, not the Rayon
variant).
"""

from __future__ import annotations

import os
import sys


def _configure_single_thread_defaults() -> None:
    """Set common BLAS/OpenMP thread env vars unless already set."""
    for key, val in (
        ("OMP_NUM_THREADS", "4"),
        ("OPENBLAS_NUM_THREADS", "4"),
        ("MKL_NUM_THREADS", "4"),
        ("NUMEXPR_NUM_THREADS", "4"),
        ("VECLIB_MAXIMUM_THREADS", "4"),
    ):
        os.environ.setdefault(key, val)


_configure_single_thread_defaults()

import argparse
import json
import time
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
from tqdm import tqdm

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from phenoms.hbond import (
    _run_baker_hubbard_rs,
    _run_baker_hubbard_rs_threads,
    _topology_donor_acceptor_lists,
)
from phenoms.io import load_and_select_residues
from phenoms.outputs import default_output_root


def _default_replicates(data_dir: Path) -> list[Path]:
    return [
        data_dir / "rep1_md_0_500_nolig_nojump_center.pdb",
        data_dir / "rep2_md_0_500_no_lig_nojump_center.pdb",
        data_dir / "rep3_md_0_500_no_lig_nojump_center.pdb",
    ]


def _write_benchmark_readme(out_dir: Path, *, n_frames: int) -> None:
    text = f"""Benchmark: kernel / API timing (replicates × first {n_frames} frames)

Rust
  Single-thread ``run_baker_hubbard``. Count = rows in output (one per geometric hit).

MDTraj
  One call per replicate: ``md.baker_hubbard(traj[:{n_frames}], periodic=False)``
  (other arguments left at **library defaults**, including ``freq=0.1`` and
  ``exclude_water=True``). Count = number of returned (donor, H, acceptor) triplets.

MDAnalysis
  ``HydrogenBondAnalysis`` with explicit selections for bondless/chargeless PDB
  inputs (protein N/H and O/N acceptors, ``d_h_cutoff=1.2`` Å). Count = rows in
  ``results.hbonds``. Default HBA distance/angle cutoffs apply unless changed in code.

Thread env (setdefault before imports): OMP_NUM_THREADS, OPENBLAS_NUM_THREADS,
MKL_NUM_THREADS, NUMEXPR_NUM_THREADS, VECLIB_MAXIMUM_THREADS → 1

Docker (1 CPU): see docker/README.md in the repository.
"""
    (out_dir / "BENCHMARK_README.txt").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Benchmark Rust vs MDTraj (baker_hubbard) vs MDAnalysis (HBA) on the same frame subset."
    )
    parser.add_argument(
        "--sub-frames",
        type=int,
        default=100,
        help="Use the first N frames of each replicate (default: 100).",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="Output directory (default: $PHENOMS_OUTPUT_DIR/benchmarks/kernel or ./phenom_outputs/benchmarks/kernel).",
    )
    parser.add_argument(
        "--rust-threads",
        type=int,
        default=4,
        help="Rust kernel threads (default: 4). Use with --cpus for parallel Rust timing.",
    )
    parser.add_argument(
        "--mda-workers",
        type=int,
        default=4,
        help="MDAnalysis HBA workers (default: 4).",
    )
    args = parser.parse_args()

    if _run_baker_hubbard_rs is None:
        print("Rust extension unavailable. Build with `maturin develop --release`.", file=sys.stderr)
        return 1
    if args.rust_threads > 1 and _run_baker_hubbard_rs_threads is None:
        print("Threaded Rust extension unavailable. Rebuild extension with current sources.", file=sys.stderr)
        return 1

    data_dir = _REPO_ROOT / "phenoms" / "test_data" / "2_24_25_hbond_heatmaps"
    pdb_files = _default_replicates(data_dir)
    missing = [p for p in pdb_files if not p.exists()]
    if missing:
        print("Missing input files:", file=sys.stderr)
        for p in missing:
            print(f"  {p}", file=sys.stderr)
        return 1

    out_dir = Path(args.out_dir) if args.out_dir else default_output_root() / "benchmarks" / "kernel"
    out_dir.mkdir(parents=True, exist_ok=True)
    n_cap = int(args.sub_frames)

    meta = {
        "thread_env": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "mdtraj_baker_hubbard_kwargs": {"periodic": False, "note": "all other args MDTraj defaults (freq=0.1, exclude_water=True, …)"},
        "mdanalysis_hba": {
            "donors_sel": "protein and element N",
            "hydrogens_sel": "protein and element H",
            "acceptors_sel": "protein and (element O or element N)",
            "d_h_cutoff_A": 1.2,
            "note": "PDB without bonds/charges; HBA default d_a_cutoff / angle apply.",
        },
    }
    (out_dir / "benchmark_config.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    rows = []
    for pdb in tqdm(pdb_files, desc="replicates", unit="replicate"):
        traj = load_and_select_residues(str(pdb), resid_range=None)
        n = min(n_cap, traj.n_frames)
        traj_n = traj[:n]

        donor_h_pairs, acceptors = _topology_donor_acceptor_lists(traj_n, backbone_only=False)
        if donor_h_pairs is None or acceptors is None:
            print(f"Skipping {pdb.name}: no donor/acceptor topology lists.")
            continue

        xyz_nm = traj_n.xyz.astype(np.float64)
        xyz_flat = np.ascontiguousarray(xyz_nm.reshape(-1))

        t0 = time.perf_counter()
        if args.rust_threads > 1:
            raw = _run_baker_hubbard_rs_threads(
                xyz_flat,
                xyz_nm.shape[0],
                xyz_nm.shape[1],
                donor_h_pairs,
                acceptors,
                int(args.rust_threads),
            )
        else:
            raw = _run_baker_hubbard_rs(
                xyz_flat,
                xyz_nm.shape[0],
                xyz_nm.shape[1],
                donor_h_pairs,
                acceptors,
            )
        rust_sec = time.perf_counter() - t0
        rust_hits = int(raw.shape[0]) if raw is not None else 0

        t0 = time.perf_counter()
        hb_md = md.baker_hubbard(traj_n, periodic=False)
        mdt_sec = time.perf_counter() - t0
        mdtraj_hits = int(len(hb_md))

        u = mda.Universe(str(pdb))
        hb_mda = HBA(
            universe=u,
            donors_sel="protein and element N",
            hydrogens_sel="protein and element H",
            acceptors_sel="protein and (element O or element N)",
            d_h_cutoff=1.2,
        )
        t0 = time.perf_counter()
        hb_mda.run(start=0, stop=n, step=1, verbose=False, n_workers=int(args.mda_workers), backend="multiprocessing")
        mda_sec = time.perf_counter() - t0
        mda_hbonds = hb_mda.results.hbonds
        mdanalysis_hits = int(len(mda_hbonds)) if mda_hbonds is not None else 0

        rows.append(
            {
                "replicate": pdb.stem,
                "frames": int(n),
                "rust_kernel_seconds": float(rust_sec),
                "mdtraj_baker_hubbard_seconds": float(mdt_sec),
                "mdanalysis_hba_seconds": float(mda_sec),
                "speedup_mdtraj_over_rust_x": float(mdt_sec / rust_sec) if rust_sec > 0 else float("inf"),
                "speedup_mda_over_rust_x": float(mda_sec / rust_sec) if rust_sec > 0 else float("inf"),
                "speedup_mda_over_mdtraj_x": float(mda_sec / mdt_sec) if mdt_sec > 0 else float("inf"),
                "rust_raw_hbond_rows": rust_hits,
                "mdtraj_persistent_triplets": mdtraj_hits,
                "mdanalysis_hba_rows": mdanalysis_hits,
                "rust_threads": int(args.rust_threads),
                "mda_workers": int(args.mda_workers),
            }
        )

    per_rep = pd.DataFrame(rows)
    per_path = out_dir / "kernel_per_replicate.csv"
    per_rep.to_csv(per_path, index=False)

    if not per_rep.empty:
        rust_total = float(per_rep["rust_kernel_seconds"].sum())
        mdt_total = float(per_rep["mdtraj_baker_hubbard_seconds"].sum())
        mda_total = float(per_rep["mdanalysis_hba_seconds"].sum())
        summary = pd.DataFrame(
            [
                {
                    "replicates": int(len(per_rep)),
                    "frames_per_replicate": int(per_rep["frames"].max()),
                    "rust_kernel_total_seconds": rust_total,
                    "mdtraj_baker_hubbard_total_seconds": mdt_total,
                    "mdanalysis_hba_total_seconds": mda_total,
                    "speedup_mdtraj_over_rust_x": mdt_total / rust_total if rust_total > 0 else float("inf"),
                    "speedup_mda_over_rust_x": mda_total / rust_total if rust_total > 0 else float("inf"),
                    "speedup_mda_over_mdtraj_x": mda_total / mdt_total if mdt_total > 0 else float("inf"),
                }
            ]
        )
    else:
        summary = pd.DataFrame()

    summary_path = out_dir / "kernel_summary.csv"
    summary.to_csv(summary_path, index=False)
    _write_benchmark_readme(out_dir, n_frames=n_cap if not per_rep.empty else n_cap)

    print(f"Kernel benchmark CSV: {per_path}")
    print(f"Kernel summary CSV:   {summary_path}")
    print(f"Notes:                {out_dir / 'BENCHMARK_README.txt'}")
    if not summary.empty:
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
