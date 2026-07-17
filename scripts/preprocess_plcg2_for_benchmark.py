#!/usr/bin/env python3
"""
Preprocess PLCg2 trajectory on host before Docker benchmarking.

Outputs a protein-only, centered, fitted trajectory capped to N frames.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md


def main() -> int:
    p = argparse.ArgumentParser(description="Preprocess PLCg2 for benchmark.")
    repo = Path(__file__).resolve().parent.parent
    p.add_argument("--traj", default=str(repo / "plcg2_data" / "wtPLCg2-1.3500ns.nc"))
    p.add_argument("--top", default=str(repo / "plcg2_data" / "wtPLCg2.prmtop"))
    p.add_argument("--sub-frames", type=int, default=500)
    p.add_argument("--out-dir", default=str(repo / "phenom_outputs" / "benchmarks" / "plcg2_preprocessed"))
    p.add_argument("--no-image", action="store_true")
    p.add_argument("--no-center", action="store_true")
    p.add_argument("--no-fit", action="store_true")
    args = p.parse_args()

    traj_path = Path(args.traj).expanduser().resolve()
    top_path = Path(args.top).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    traj = md.load(str(traj_path), top=str(top_path))
    n = min(int(args.sub_frames), traj.n_frames)
    if n <= 0:
        raise ValueError("--sub-frames must be > 0")
    traj = traj[:n]
    protein_idx = traj.topology.select("protein")
    if len(protein_idx) == 0:
        raise ValueError("No protein atoms found.")
    traj = traj.atom_slice(protein_idx)

    if not args.no_image:
        try:
            traj.image_molecules(inplace=True)
        except Exception:
            pass
    if not args.no_center:
        traj.center_coordinates()
    if not args.no_fit and traj.n_frames > 1:
        traj.superpose(traj, frame=0, atom_indices=list(range(traj.n_atoms)))

    out_pdb = out_dir / f"plcg2_preprocessed_{n}f.pdb"
    out_top = out_dir / "plcg2_preprocessed_top.pdb"
    out_xtc = out_dir / f"plcg2_preprocessed_{n}f.xtc"
    traj.save_pdb(str(out_pdb))
    traj[0].save_pdb(str(out_top))
    traj.save_xtc(str(out_xtc))
    print(f"Wrote preprocessed PDB: {out_pdb}")
    print(f"Wrote preprocessed topology PDB: {out_top}")
    print(f"Wrote preprocessed XTC: {out_xtc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
