#!/usr/bin/env python3
"""
Preprocess GROMACS outputs into a centered, fitted, protein-only multi-frame PDB.

This mirrors a common 2-step trjconv workflow:
1) center + remove PBC jumps on the full system
2) fit on protein and write frame-thinned protein-only PDB
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def _pick_one(run_dir: Path, pattern: str) -> Path | None:
    hits = sorted(run_dir.glob(pattern))
    return hits[0] if hits else None


def _run_gmx(cmd: list[str], stdin_text: str) -> None:
    subprocess.run(
        cmd,
        input=stdin_text.encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )


def main() -> int:
    p = argparse.ArgumentParser(description="Export centered/fitted protein-only PDB frames from GROMACS outputs.")
    p.add_argument("--run-dir", required=True, help="Directory containing *_500ns.xtc/tpr (or *_500ns.gro).")
    p.add_argument("--out-dir", required=True, help="Output directory.")
    p.add_argument("--frame-dt-ps", type=int, default=1000, help="Frame interval in ps for output PDB (default 1000).")
    p.add_argument("--start-ps", type=int, default=1000, help="Start time in ps (default 1000).")
    p.add_argument("--end-ps", type=int, default=500000, help="End time in ps (default 500000).")
    p.add_argument("--xtc-pattern", default="*_500ns.xtc")
    p.add_argument("--tpr-pattern", default="*_500ns.tpr")
    p.add_argument("--gro-pattern", default="*_500ns.gro")
    p.add_argument("--gmx-bin", default="gmx", help="GROMACS executable (default: gmx).")
    args = p.parse_args()

    run_dir = Path(args.run_dir).expanduser().resolve()
    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    run_name = run_dir.name

    traj = _pick_one(run_dir, args.xtc_pattern)
    tpr = _pick_one(run_dir, args.tpr_pattern)
    gro = _pick_one(run_dir, args.gro_pattern)

    out_pdb = out_dir / f"{run_name}_protein_centered_fit_nopbc.pdb"
    tmp_xtc = run_dir / f"{run_name}_centered_tmp.xtc"

    if traj is not None and tpr is not None:
        _run_gmx(
            [
                args.gmx_bin,
                "trjconv",
                "-f",
                str(traj),
                "-s",
                str(tpr),
                "-o",
                str(tmp_xtc),
                "-pbc",
                "mol",
                "-ur",
                "compact",
                "-center",
            ],
            "Protein\nSystem\n",
        )
        _run_gmx(
            [
                args.gmx_bin,
                "trjconv",
                "-f",
                str(tmp_xtc),
                "-s",
                str(tpr),
                "-o",
                str(out_pdb),
                "-fit",
                "rot+trans",
                "-b",
                str(args.start_ps),
                "-e",
                str(args.end_ps),
                "-dt",
                str(args.frame_dt_ps),
            ],
            "Protein\nProtein\n",
        )
        if tmp_xtc.exists():
            tmp_xtc.unlink()
    elif gro is not None:
        _run_gmx(
            [
                args.gmx_bin,
                "trjconv",
                "-f",
                str(gro),
                "-s",
                str(gro),
                "-o",
                str(out_pdb),
            ],
            "Protein\n",
        )
    else:
        raise FileNotFoundError(
            f"No matching trajectory/topology files in {run_dir} "
            f"(patterns: {args.xtc_pattern}, {args.tpr_pattern}, {args.gro_pattern})."
        )

    print(f"Wrote: {out_pdb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
