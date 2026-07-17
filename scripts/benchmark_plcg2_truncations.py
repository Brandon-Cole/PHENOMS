#!/usr/bin/env python3
"""
Benchmark Rust vs MDTraj vs MDAnalysis on PLCg2 test data after preprocessing.

Pipeline:
1) preprocess input trajectory (protein-only, image molecules, center, fit, frame cap)
2) export preprocessed trajectory artifacts
3) benchmark full protein + truncation windows from the same preprocessed data

Outputs:
- preprocessed trajectory files
- per-window benchmark CSV
- long-format CSV for curve plotting
- summary CSV
- PDF table
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import tempfile
import time
from pathlib import Path


def _configure_thread_defaults() -> None:
    for key, val in (
        ("OMP_NUM_THREADS", "4"),
        ("OPENBLAS_NUM_THREADS", "4"),
        ("MKL_NUM_THREADS", "4"),
        ("NUMEXPR_NUM_THREADS", "4"),
        ("VECLIB_MAXIMUM_THREADS", "4"),
    ):
        os.environ.setdefault(key, val)


_configure_thread_defaults()

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
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
from phenoms.outputs import default_output_root


def _preprocess_trajectory(
    traj_path: Path,
    top_path: Path,
    *,
    sub_frames: int,
    out_dir: Path,
    apply_imaging: bool,
    center: bool,
    fit: bool,
) -> tuple[md.Trajectory, Path, Path]:
    traj_raw = md.load(str(traj_path), top=str(top_path))
    if sub_frames <= 0:
        raise ValueError("--sub-frames must be > 0.")
    n = min(int(sub_frames), traj_raw.n_frames)
    traj = traj_raw[:n]
    protein_idx = traj.topology.select("protein")
    if len(protein_idx) == 0:
        raise ValueError("No protein atoms found after loading trajectory/topology.")
    traj = traj.atom_slice(protein_idx)

    if apply_imaging:
        try:
            traj.image_molecules(inplace=True)
        except Exception:
            # Best effort across formats; missing bond info can break imaging.
            pass
    if center:
        traj.center_coordinates()
    if fit and traj.n_frames > 1:
        traj.superpose(traj, frame=0, atom_indices=list(range(traj.n_atoms)))

    out_pdb = out_dir / f"plcg2_preprocessed_{n}f.pdb"
    out_xtc = out_dir / f"plcg2_preprocessed_{n}f.xtc"
    traj.save_pdb(str(out_pdb))
    traj.save_xtc(str(out_xtc))
    return traj, out_pdb, out_xtc


def _truncation_windows(
    topology: md.Topology,
    *,
    n_windows: int,
    truncation_step: int,
) -> list[dict]:
    if n_windows < 1:
        raise ValueError("--n-windows must be >= 1.")
    if truncation_step <= 0:
        raise ValueError("--truncation-step must be > 0.")

    residues = [r for r in topology.residues]
    if not residues:
        raise ValueError("No residues found in preprocessed trajectory.")
    total = len(residues)

    windows = [
        {
            "window": "Full",
            "start_pos_1b": 1,
            "end_pos_1b": total,
            "residue_indices": [r.index for r in residues],
            "residue_count": total,
        }
    ]

    for i in range(1, n_windows):
        keep = total - (i * truncation_step)
        if keep <= 0:
            break
        sub_res = residues[:keep]
        start_pos = 1
        end_pos = len(sub_res)
        windows.append(
            {
                "window": f"Trunc{i}_{start_pos}_{end_pos}",
                "start_pos_1b": start_pos,
                "end_pos_1b": end_pos,
                "residue_indices": [r.index for r in sub_res],
                "residue_count": len(sub_res),
            }
        )
    return windows


def _atom_indices_for_residues(traj: md.Trajectory, residue_indices: list[int]) -> np.ndarray:
    residue_set = set(residue_indices)
    atom_idx = [a.index for a in traj.topology.atoms if a.residue.index in residue_set]
    return np.array(atom_idx, dtype=int)


def _atom_indices_for_residues_in_top(topology: md.Topology, residue_indices: list[int]) -> np.ndarray:
    residue_set = set(residue_indices)
    atom_idx = [a.index for a in topology.atoms if a.residue.index in residue_set]
    return np.array(atom_idx, dtype=int)


def _run_methods(
    traj_n: md.Trajectory,
    mda_universe: mda.Universe,
    *,
    rust_threads: int,
    mda_workers: int,
    run_mdtraj: bool = True,
    mdtraj_frame_cap: int | None = None,
) -> dict[str, float | int]:
    donor_h_pairs, acceptors = _topology_donor_acceptor_lists(traj_n, backbone_only=False)
    if donor_h_pairs is None or acceptors is None:
        return {
            "rust_seconds": float("nan"),
            "mdtraj_seconds": float("nan"),
            "mdanalysis_seconds": float("nan"),
            "rust_hits": 0,
            "mdtraj_hits": 0,
            "mdanalysis_hits": 0,
        }

    xyz_nm = traj_n.xyz.astype(np.float64)
    xyz_flat = np.ascontiguousarray(xyz_nm.reshape(-1))

    t0 = time.perf_counter()
    if rust_threads > 1:
        raw = _run_baker_hubbard_rs_threads(
            xyz_flat,
            xyz_nm.shape[0],
            xyz_nm.shape[1],
            donor_h_pairs,
            acceptors,
            int(rust_threads),
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
    # Rust output can be large at 1000 frames; release it before next method.
    del raw, xyz_nm, xyz_flat
    gc.collect()

    hb_md = None
    mdt_sec = float("nan")
    mdtraj_hits = 0
    if run_mdtraj:
        mdtraj_input = traj_n if mdtraj_frame_cap is None else traj_n[: min(int(mdtraj_frame_cap), traj_n.n_frames)]
        t0 = time.perf_counter()
        hb_md = md.baker_hubbard(mdtraj_input, periodic=False)
        mdt_sec = time.perf_counter() - t0
        mdtraj_hits = int(len(hb_md))

    hb_mda = HBA(
        universe=mda_universe,
        # Match Rust donor/acceptor element classes (backbone_only=False): N/O/S.
        donors_sel="protein and (element N or element O or element S)",
        hydrogens_sel="protein and element H",
        acceptors_sel="protein and (element N or element O or element S)",
        d_h_cutoff=1.2,
        # Use BH-style angular criterion; MDA also requires a D-A cutoff.
        d_a_cutoff=3.0,
        d_h_a_angle_cutoff=120.0,
        update_selections=False,
    )
    t0 = time.perf_counter()
    hb_mda.run(
        start=0,
        stop=traj_n.n_frames,
        step=1,
        verbose=False,
        n_workers=int(mda_workers),
        backend="multiprocessing",
    )
    mda_sec = time.perf_counter() - t0
    mda_hbonds = hb_mda.results.hbonds
    mdanalysis_hits = int(len(mda_hbonds)) if mda_hbonds is not None else 0
    # Free potentially large MDA/MDTraj intermediate arrays aggressively.
    del hb_mda, mda_hbonds, hb_md

    return {
        "rust_seconds": float(rust_sec),
        "mdtraj_seconds": float(mdt_sec),
        "mdanalysis_seconds": float(mda_sec),
        "rust_hits": rust_hits,
        "mdtraj_hits": mdtraj_hits,
        "mdanalysis_hits": mdanalysis_hits,
    }


def _render_pdf_table(df: pd.DataFrame, out_pdf: Path) -> None:
    display = df[
        [
            "window",
            "residue_count",
            "frames",
            "rust_seconds",
            "mdtraj_seconds",
            "mdanalysis_seconds",
            "speedup_mdtraj_over_rust_x",
            "speedup_mda_over_rust_x",
        ]
    ].copy()
    for col in (
        "rust_seconds",
        "mdtraj_seconds",
        "mdanalysis_seconds",
        "speedup_mdtraj_over_rust_x",
        "speedup_mda_over_rust_x",
    ):
        display[col] = display[col].map(lambda v: f"{v:.3f}" if pd.notna(v) else "-")

    rows = display.values.tolist()
    fig_h = max(3.0, 1.2 + 0.38 * len(rows))
    fig, ax = plt.subplots(figsize=(13, fig_h))
    ax.axis("off")
    tbl = ax.table(
        cellText=rows,
        colLabels=[
            "Window",
            "Residues",
            "Frames",
            "Rust (s)",
            "MDTraj (s)",
            "MDAnalysis (s)",
            "MDTraj/Rust",
            "MDA/Rust",
        ],
        cellLoc="left",
        loc="center",
        colWidths=[0.15, 0.10, 0.08, 0.11, 0.11, 0.12, 0.11, 0.11],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1.0, 1.35)
    for c in range(8):
        cell = tbl[(0, c)]
        cell.set_text_props(weight="bold")
        cell.set_facecolor("#EDEDED")
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser(
        description="PLCg2 benchmark with full protein + 9 equal truncations."
    )
    p.add_argument(
        "--traj",
        type=str,
        default=str(_REPO_ROOT / "plcg2_data" / "wtPLCg2-1.3500ns.nc"),
        help="Trajectory path (default: plcg2_data/wtPLCg2-1.3500ns.nc).",
    )
    p.add_argument(
        "--top",
        type=str,
        default=str(_REPO_ROOT / "plcg2_data" / "wtPLCg2.prmtop"),
        help="Topology path (default: plcg2_data/wtPLCg2.prmtop).",
    )
    p.add_argument(
        "--preprocessed-pdb",
        type=str,
        default=None,
        help="If set, skip preprocess and benchmark directly from a preprocessed multi-frame protein-only PDB.",
    )
    p.add_argument(
        "--preprocessed-traj",
        type=str,
        default=None,
        help="If set with --preprocessed-top, skip preprocess and benchmark directly from preprocessed trajectory (e.g., .xtc).",
    )
    p.add_argument(
        "--preprocessed-top",
        type=str,
        default=None,
        help="Topology file for --preprocessed-traj (e.g., single-frame preprocessed .pdb).",
    )
    p.add_argument(
        "--sub-frames",
        type=int,
        default=500,
        help="Use first N frames after preprocessing (default: 500; ~0.5 microseconds if 1 ns/frame).",
    )
    p.add_argument(
        "--n-windows",
        type=int,
        default=10,
        help="Total windows including Full (default: 10).",
    )
    p.add_argument(
        "--truncation-step",
        type=int,
        default=100,
        help="Residues removed from the C-terminus per truncation (default: 100).",
    )
    p.add_argument(
        "--rust-threads",
        type=int,
        default=4,
        help="Rust kernel threads (default: 4).",
    )
    p.add_argument(
        "--mda-workers",
        type=int,
        default=4,
        help="MDAnalysis workers (default: 4).",
    )
    p.add_argument(
        "--skip-mdtraj",
        action="store_true",
        help="Skip MDTraj baker_hubbard stage (MDTraj columns become NA).",
    )
    p.add_argument(
        "--mdtraj-sub-frames",
        type=int,
        default=None,
        help="Optional frame cap only for MDTraj stage (e.g., 250) while Rust/MDA use --sub-frames.",
    )
    p.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="Output directory (default: $PHENOMS_OUTPUT_DIR/benchmarks/plcg2_truncations).",
    )
    p.add_argument(
        "--no-image",
        action="store_true",
        help="Disable imaging/PBC cleanup during preprocessing.",
    )
    p.add_argument(
        "--no-center",
        action="store_true",
        help="Disable centering during preprocessing.",
    )
    p.add_argument(
        "--no-fit",
        action="store_true",
        help="Disable frame superposition during preprocessing.",
    )
    args = p.parse_args()

    if _run_baker_hubbard_rs is None:
        print("Rust extension unavailable. Build with `maturin develop --release`.", file=sys.stderr)
        return 1
    if args.rust_threads > 1 and _run_baker_hubbard_rs_threads is None:
        print("Threaded Rust extension unavailable. Rebuild extension with current sources.", file=sys.stderr)
        return 1

    out_dir = Path(args.out_dir).expanduser().resolve() if args.out_dir else default_output_root() / "benchmarks" / "plcg2_truncations"
    out_dir.mkdir(parents=True, exist_ok=True)
    prep_pdb: Path | None = None
    prep_xtc: Path | None = None

    if args.preprocessed_pdb and (args.preprocessed_traj or args.preprocessed_top):
        print("Use either --preprocessed-pdb OR (--preprocessed-traj and --preprocessed-top), not both.", file=sys.stderr)
        return 1

    if args.preprocessed_traj and not args.preprocessed_top:
        print("--preprocessed-top is required with --preprocessed-traj.", file=sys.stderr)
        return 1

    if args.preprocessed_pdb:
        prep_pdb = Path(args.preprocessed_pdb).expanduser().resolve()
        if not prep_pdb.exists():
            print(f"Preprocessed PDB not found: {prep_pdb}", file=sys.stderr)
            return 1
        top_for_windows = md.load_topology(str(prep_pdb))
        source_traj_path = prep_pdb
        source_top_path = prep_pdb
    elif args.preprocessed_traj:
        prep_traj = Path(args.preprocessed_traj).expanduser().resolve()
        prep_top = Path(args.preprocessed_top).expanduser().resolve()
        if not prep_traj.exists() or not prep_top.exists():
            print("Preprocessed traj/top not found.", file=sys.stderr)
            print(f"  traj: {prep_traj} (exists={prep_traj.exists()})", file=sys.stderr)
            print(f"  top:  {prep_top} (exists={prep_top.exists()})", file=sys.stderr)
            return 1
        top_for_windows = md.load_topology(str(prep_top))
        source_traj_path = prep_traj
        source_top_path = prep_top
        prep_pdb = prep_top
    else:
        traj_path = Path(args.traj).expanduser().resolve()
        top_path = Path(args.top).expanduser().resolve()
        if not traj_path.exists() or not top_path.exists():
            print("Missing PLCg2 input files.", file=sys.stderr)
            print(f"  traj: {traj_path} (exists={traj_path.exists()})", file=sys.stderr)
            print(f"  top:  {top_path} (exists={top_path.exists()})", file=sys.stderr)
            return 1
        traj_n, prep_pdb, prep_xtc = _preprocess_trajectory(
            traj_path,
            top_path,
            sub_frames=int(args.sub_frames),
            out_dir=out_dir,
            apply_imaging=not args.no_image,
            center=not args.no_center,
            fit=not args.no_fit,
        )
        top_for_windows = traj_n.topology
        source_traj_path = prep_xtc
        source_top_path = prep_pdb

    n = int(args.sub_frames)
    windows = _truncation_windows(
        top_for_windows,
        n_windows=int(args.n_windows),
        truncation_step=int(args.truncation_step),
    )
    rows: list[dict] = []
    for w in windows:
        window_name = w["window"]
        residue_indices = w["residue_indices"]
        atom_idx = _atom_indices_for_residues_in_top(top_for_windows, residue_indices)
        if len(atom_idx) == 0:
            continue
        sub_all = md.load(str(source_traj_path), top=str(source_top_path), atom_indices=atom_idx)
        sub = sub_all[: min(int(args.sub_frames), sub_all.n_frames)]
        # Avoid storing one trajectory file per truncation; MDAnalysis consumes
        # a short-lived temp PDB for this window only.
        with tempfile.TemporaryDirectory(prefix="plcg2_window_") as td:
            window_pdb = Path(td) / f"{window_name}.pdb"
            sub.save_pdb(str(window_pdb))
            mda_u = mda.Universe(str(window_pdb))
            m = _run_methods(
                sub,
                mda_u,
                rust_threads=int(args.rust_threads),
                mda_workers=int(args.mda_workers),
                run_mdtraj=not bool(args.skip_mdtraj),
                mdtraj_frame_cap=args.mdtraj_sub_frames,
            )
        sub_n_frames = int(sub.n_frames)
        sub_n_atoms = int(sub.n_atoms)
        del sub_all, sub, mda_u
        gc.collect()
        rust_sec = m["rust_seconds"]
        mdt_sec = m["mdtraj_seconds"]
        mda_sec = m["mdanalysis_seconds"]
        rows.append(
            {
                "window": window_name,
                "residue_start_pos_1b": int(w["start_pos_1b"]),
                "residue_end_pos_1b": int(w["end_pos_1b"]),
                "residue_count": int(w["residue_count"]),
                "frames": sub_n_frames,
                "atoms": sub_n_atoms,
                "rust_seconds": float(rust_sec),
                "mdtraj_seconds": float(mdt_sec),
                "mdanalysis_seconds": float(mda_sec),
                "frames_rust": sub_n_frames,
                "frames_mda": sub_n_frames,
                "frames_mdtraj": int(min(sub_n_frames, int(args.mdtraj_sub_frames))) if (not args.skip_mdtraj and args.mdtraj_sub_frames is not None) else (sub_n_frames if not args.skip_mdtraj else 0),
                "speedup_mdtraj_over_rust_x": float(mdt_sec / rust_sec) if rust_sec and rust_sec > 0 else float("nan"),
                "speedup_mda_over_rust_x": float(mda_sec / rust_sec) if rust_sec and rust_sec > 0 else float("nan"),
                "rust_hits": int(m["rust_hits"]),
                "mdtraj_hits": (float("nan") if args.skip_mdtraj else int(m["mdtraj_hits"])),
                "mdanalysis_hits": int(m["mdanalysis_hits"]),
                "Rust_H_Bond_Total": int(m["rust_hits"]),
                "MDTraj_H_Bond_Total": (float("nan") if args.skip_mdtraj else int(m["mdtraj_hits"])),
                "MDAnalysis_H_Bond_Total": int(m["mdanalysis_hits"]),
                "rust_threads": int(args.rust_threads),
                "mda_workers": int(args.mda_workers),
            }
        )

    per_window = pd.DataFrame(rows)
    per_window_path = out_dir / "plcg2_truncations_per_window.csv"
    per_window.to_csv(per_window_path, index=False)

    long_rows = []
    for _, r in per_window.iterrows():
        long_rows.extend(
            [
                {
                    "window": r["window"],
                    "residue_count": int(r["residue_count"]),
                    "method": "rust",
                    "seconds": float(r["rust_seconds"]),
                    "Total_HBonds_Method": int(r["rust_hits"]),
                },
                {
                    "window": r["window"],
                    "residue_count": int(r["residue_count"]),
                    "method": "mdtraj",
                    "seconds": float(r["mdtraj_seconds"]),
                    "Total_HBonds_Method": (None if pd.isna(r["mdtraj_hits"]) else int(r["mdtraj_hits"])),
                },
                {
                    "window": r["window"],
                    "residue_count": int(r["residue_count"]),
                    "method": "mdanalysis",
                    "seconds": float(r["mdanalysis_seconds"]),
                    "Total_HBonds_Method": int(r["mdanalysis_hits"]),
                },
            ]
        )
    long_df = pd.DataFrame(long_rows)
    long_path = out_dir / "plcg2_truncations_long.csv"
    long_df.to_csv(long_path, index=False)

    summary = pd.DataFrame(
        [
            {
                "n_windows": int(len(per_window)),
                "frames_per_window": int(n),
                "rust_total_seconds": float(per_window["rust_seconds"].sum()),
                "mdtraj_total_seconds": float(per_window["mdtraj_seconds"].sum()),
                "mdanalysis_total_seconds": float(per_window["mdanalysis_seconds"].sum()),
            }
        ]
    )
    summary_path = out_dir / "plcg2_truncations_summary.csv"
    summary.to_csv(summary_path, index=False)

    pdf_path = out_dir / "plcg2_truncations_table.pdf"
    _render_pdf_table(per_window, pdf_path)

    config = {
        "traj": None if args.preprocessed_pdb else str(Path(args.traj).expanduser().resolve()),
        "top": None if args.preprocessed_pdb else str(Path(args.top).expanduser().resolve()),
        "preprocessed_input_pdb": str(prep_pdb) if prep_pdb is not None else None,
        "preprocessed_input_traj": str(Path(args.preprocessed_traj).expanduser().resolve()) if args.preprocessed_traj else None,
        "preprocessed_input_top": str(Path(args.preprocessed_top).expanduser().resolve()) if args.preprocessed_top else None,
        "sub_frames": int(args.sub_frames),
        "n_windows": int(args.n_windows),
        "truncation_step": int(args.truncation_step),
        "preprocessed_pdb": str(prep_pdb) if prep_pdb is not None else None,
        "preprocessed_xtc": str(prep_xtc) if prep_xtc is not None else None,
        "preprocess": {
            "ran_in_script": bool(not args.preprocessed_pdb),
            "protein_only": True if not args.preprocessed_pdb else None,
            "apply_imaging": bool(not args.no_image) if not args.preprocessed_pdb else None,
            "center": bool(not args.no_center) if not args.preprocessed_pdb else None,
            "fit": bool(not args.no_fit) if not args.preprocessed_pdb else None,
        },
        "rust_threads": int(args.rust_threads),
        "mda_workers": int(args.mda_workers),
        "skip_mdtraj": bool(args.skip_mdtraj),
        "mdtraj_sub_frames": int(args.mdtraj_sub_frames) if args.mdtraj_sub_frames is not None else None,
        "thread_env": {
            k: os.environ.get(k)
            for k in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        },
    }
    (out_dir / "benchmark_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    if prep_pdb is not None:
        print(f"Preprocessed PDB: {prep_pdb}")
    if prep_xtc is not None:
        print(f"Preprocessed XTC: {prep_xtc}")
    print(f"Per-window CSV:   {per_window_path}")
    print(f"Long CSV:       {long_path}")
    print(f"Summary CSV:    {summary_path}")
    print(f"PDF table:      {pdf_path}")
    if not per_window.empty:
        print(per_window.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
