#!/usr/bin/env python3
"""
Renumber multiple mobile PDBs to a common reference numbering scheme.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def main() -> int:
    repo = Path(__file__).resolve().parent.parent
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    os.environ.setdefault("MPLBACKEND", "Agg")

    from phenoms.cleanup import renumber_many_to_reference

    p = argparse.ArgumentParser(description="Renumber many mobile PDBs to a reference numbering scheme.")
    p.add_argument("--reference", required=True, help="Reference PDB (target numbering).")
    p.add_argument(
        "--mobile-dir",
        required=True,
        help="Directory containing mobile PDBs to renumber.",
    )
    p.add_argument(
        "--mobile-glob",
        default="*.pdb",
        help="Glob pattern for mobile files inside --mobile-dir (default: *.pdb).",
    )
    p.add_argument("--output-dir", required=True, help="Directory for renumbered outputs.")
    p.add_argument("--chain", default=None, help="Chain ID if needed (e.g. A).")
    p.add_argument("--suffix", default="_renumbered_to_ref.pdb")
    p.add_argument("--match", type=int, default=2)
    p.add_argument("--mismatch", type=int, default=-1)
    p.add_argument("--gap", type=int, default=-2)
    p.add_argument(
        "--fill-unmapped-mobile",
        action="store_true",
        help="Assign consecutive ref-like numbering to mobile-only gaps/overhangs.",
    )
    args = p.parse_args()

    mobile_dir = Path(args.mobile_dir).expanduser().resolve()
    mobiles = sorted(mobile_dir.glob(args.mobile_glob))
    if not mobiles:
        raise FileNotFoundError(f"No files matched {args.mobile_glob} in {mobile_dir}")

    reports = renumber_many_to_reference(
        args.reference,
        mobiles,
        args.output_dir,
        chain_id=args.chain,
        suffix=args.suffix,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
        fill_unmapped_mobile=args.fill_unmapped_mobile,
    )
    print(f"Renumbered {len(reports)} files to {args.output_dir}")
    for stem, rep in reports.items():
        print(
            f"- {stem}: mapped={rep.n_mapped}, mismatches={rep.n_mismatch_at_aligned_positions}, "
            f"filled_unmapped={rep.n_filled_unmapped_mobile}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
