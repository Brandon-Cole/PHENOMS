#!/usr/bin/env python3
"""
Rust vs MDTraj parity check for the difference-plot underlying logic.

Computes ligand − no-ligand protection difference:
- Occupancy-based imputation: impute_threshold = 0.4 (in compare_two_sets)
- Delta threshold (visual logic): diff_threshold = 0.2 (20% delta)

Then compares, for each residue number:
- Difference_rust vs Difference_mdtraj
- Whether the residue would be labeled (|Difference| >= 0.2) in each backend

This is meant to explain why lifetime/protection-derived plots might drift
between backends: bond label sets often differ slightly, and that changes
per-residue occupancy/differences.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from phenoms.comparison import compare_two_sets
from phenoms.simulation import SimulationSet, default_n_jobs


def _select_test_pdbs(data_dir: Path, use_all_reps: bool) -> tuple[list[Path], list[Path]]:
    """
    Default subset: 1 replicate per set for speed.
    If use_all_reps, use 3 reps per set like the normal pipeline.
    """
    if use_all_reps:
        set_a = [
            data_dir / "rep1_md_0_500_nolig_nojump_center.pdb",
            data_dir / "rep2_md_0_500_no_lig_nojump_center.pdb",
            data_dir / "rep3_md_0_500_no_lig_nojump_center.pdb",
        ]
        set_b = sorted(data_dir.glob("rep*_md_0_300_55085_nojump_center.pdb"))
        return set_a, set_b

    # Quick parity subset (1 rep per set)
    set_a = [data_dir / "rep1_md_0_500_nolig_nojump_center.pdb"]
    set_b = [data_dir / "rep2_md_0_300_55085_nojump_center.pdb"]
    return set_a, set_b


def _run_backend(
    set_a_pdbs: list[Path],
    set_b_pdbs: list[Path],
    *,
    sub_frames: int | None,
    use_rust: bool,
    n_jobs: int | None,
):
    sim_a = SimulationSet([str(p) for p in set_a_pdbs], sub_frames=sub_frames).run(
        n_jobs=n_jobs,
        use_rust=use_rust,
    )
    sim_b = SimulationSet([str(p) for p in set_b_pdbs], sub_frames=sub_frames).run(
        n_jobs=n_jobs,
        use_rust=use_rust,
    )

    # Protection difference per donor residue.
    # With flip_difference=True: Difference = lig − no_lig
    comp = compare_two_sets(
        sim_a.get_pivot_tables(),
        sim_b.get_pivot_tables(),
        label_a="no_lig",
        label_b="lig",
        impute_threshold=0.4,
        flip_difference=True,
        donor_aggregation="sum",
    )
    return sim_a, sim_b, comp


def main() -> int:
    parser = argparse.ArgumentParser(description="Rust vs MDTraj difference-plot parity check.")
    parser.add_argument(
        "--sub-frames",
        type=int,
        default=100,
        help="Frames per replicate to process. Use 0 or omit for full trajectory.",
    )
    parser.add_argument(
        "--use-all-reps",
        action="store_true",
        help="Use 3 replicates per set (slower).",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="MDTraj workers. Omit for default (CPUs minus 2).",
    )
    parser.add_argument(
        "--diff-threshold",
        type=float,
        default=0.2,
        help="Delta threshold used for labeling / delta-imputation (20% default).",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parent.parent
    data_dir = repo / "phenoms" / "test_data" / "2_24_25_hbond_heatmaps"
    if not data_dir.exists():
        print(f"Missing test data directory: {data_dir}", file=sys.stderr)
        return 1

    if args.sub_frames in (0,):
        sub_frames = None
    else:
        sub_frames = args.sub_frames

    n_jobs = args.n_jobs if args.n_jobs is not None else default_n_jobs()

    set_a_pdbs, set_b_pdbs = _select_test_pdbs(data_dir, use_all_reps=args.use_all_reps)
    missing = [p for p in set_a_pdbs + set_b_pdbs if not p.exists()]
    if missing:
        print("Missing test PDB(s):", file=sys.stderr)
        for p in missing:
            print(f"  {p}", file=sys.stderr)
        return 1

    print("Running parity check on:")
    print("  set_a (no_lig):")
    for p in set_a_pdbs:
        print(f"    - {p.name}")
    print("  set_b (lig):")
    for p in set_b_pdbs:
        print(f"    - {p.name}")
    print(f"  sub_frames={sub_frames if sub_frames is not None else 'all'}, n_jobs={n_jobs}")

    # Ensure Rust extension import is possible (otherwise, hard-fail the parity check).
    try:
        import phenoms_hbond_rs  # noqa: F401
    except ImportError as e:
        print("Rust extension not found. Run `maturin develop`.", file=sys.stderr)
        return 1

    sim_a_rust, sim_b_rust, comp_rust = _run_backend(
        set_a_pdbs,
        set_b_pdbs,
        sub_frames=sub_frames,
        use_rust=True,
        n_jobs=n_jobs,
    )
    sim_a_md, sim_b_md, comp_md = _run_backend(
        set_a_pdbs,
        set_b_pdbs,
        sub_frames=sub_frames,
        use_rust=False,
        n_jobs=n_jobs,
    )

    # Bond label set comparison (quick explanation for why values drift).
    bonds_a_rust = set(sim_a_rust.get_bond_labels_sorted())
    bonds_b_rust = set(sim_b_rust.get_bond_labels_sorted())
    bonds_a_md = set(sim_a_md.get_bond_labels_sorted())
    bonds_b_md = set(sim_b_md.get_bond_labels_sorted())

    print("\nBond label set overlap (reason for drift):")
    def _overlap(name: str, set_r: set[str], set_m: set[str]) -> None:
        inter = set_r & set_m
        only_r = set_r - set_m
        only_m = set_m - set_r
        print(
            f"  {name}: rust={len(set_r)} mdtraj={len(set_m)} overlap={len(inter)} "
            f"only_rust={len(only_r)} only_mdtraj={len(only_m)}"
        )
        if len(inter) == 0:
            print("    (no overlap; differences expected)")

    _overlap("set_a (no_lig)", bonds_a_rust, bonds_a_md)
    _overlap("set_b (lig)", bonds_b_rust, bonds_b_md)

    # Compare difference values.
    # Columns from compare_two_sets: no_lig, lig, Difference, Residue Number.
    merged = comp_rust.merge(
        comp_md,
        on="Residue Number",
        suffixes=("_rust", "_mdtraj"),
        how="outer",
    ).sort_values("Residue Number")

    merged["Difference_rust"] = merged["Difference_rust"].astype(float)
    merged["Difference_mdtraj"] = merged["Difference_mdtraj"].astype(float)
    merged["abs_diff"] = (merged["Difference_rust"] - merged["Difference_mdtraj"]).abs()

    # Apply the delta threshold logic (visual labeling / delta-impute).
    merged["label_rust"] = merged["Difference_rust"].abs() >= args.diff_threshold
    merged["label_mdtraj"] = merged["Difference_mdtraj"].abs() >= args.diff_threshold
    merged["label_mismatch"] = merged["label_rust"] ^ merged["label_mdtraj"]

    n_mismatch = int(merged["label_mismatch"].sum())
    max_abs = float(merged["abs_diff"].max()) if len(merged) else 0.0
    mean_abs = float(merged["abs_diff"].mean()) if len(merged) else 0.0

    print("\nDifference curve parity (lig − no_lig):")
    print(f"  residues compared: {len(merged)}")
    print(f"  label mismatches at |delta| >= {args.diff_threshold}: {n_mismatch}")
    print(f"  abs(Diff_rust − Diff_mdtraj): mean={mean_abs:.4f} max={max_abs:.4f}")

    if len(merged) >= 2:
        corr = float(np.corrcoef(merged["Difference_rust"], merged["Difference_mdtraj"])[0, 1])
        print(f"  correlation of raw Difference values: {corr:.4f}")

    # Show top residues by mismatch.
    top = merged.sort_values("abs_diff", ascending=False).head(10)[
        ["Residue Number", "Difference_rust", "Difference_mdtraj", "abs_diff", "label_rust", "label_mdtraj"]
    ]
    print("\nTop 10 residues by abs difference:")
    for _, row in top.iterrows():
        print(
            f"  resid {int(row['Residue Number'])}: rust={row['Difference_rust']:.3f} "
            f"mdtraj={row['Difference_mdtraj']:.3f} abs={row['abs_diff']:.3f} "
            f"label_rust={bool(row['label_rust'])} label_mdtraj={bool(row['label_mdtraj'])}"
        )

    out_dir = repo / "phenoms" / "test_outputs" / "rust_mdtraj_parity"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"parity_difference_subframes_{'all' if sub_frames is None else sub_frames}.csv"
    merged.to_csv(out_csv, index=False)
    print(f"\nWrote comparison table: {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

