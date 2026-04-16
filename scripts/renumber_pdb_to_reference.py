#!/usr/bin/env python3
"""
Renumber a mobile PDB (e.g. mutant trajectory) to match reference residue numbers
via sequence alignment. Use before PHENOMS when WT/mut use different PDB numbering.

Example (EGFR WT vs Δ747–749, same chain implicit):

    conda activate phenoms
    cd /path/to/PHENOMS
    python scripts/renumber_pdb_to_reference.py \\
        --reference protein_pdb_exports/3_19_26_wt_rep_1_protein_centered_fit_nopbc.pdb \\
        --mobile protein_pdb_exports/3_19_26_del747-749_rep1_protein_centered_fit_nopbc.pdb \\
        --output phenoms/test_outputs/renumbered/del_rep1_renumbered.pdb

If multiple chains exist, pass --chain A (or the appropriate chain id).
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

    from phenoms.cleanup import align_and_renumber_pdb

    p = argparse.ArgumentParser(
        description="Align mobile PDB sequence to reference and rewrite residue numbers."
    )
    p.add_argument("--reference", required=True, help="Reference PDB (target numbering).")
    p.add_argument("--mobile", required=True, help="Mobile PDB to renumber.")
    p.add_argument("--output", required=True, help="Output PDB path.")
    p.add_argument(
        "--chain",
        default=None,
        help="Chain ID to align (required if PDB has multiple chains).",
    )
    p.add_argument("--match", type=int, default=2, help="NW match score (default 2).")
    p.add_argument("--mismatch", type=int, default=-1, help="NW mismatch score (default -1).")
    p.add_argument("--gap", type=int, default=-2, help="NW gap penalty (default -2).")
    p.add_argument(
        "--fill-unmapped-mobile",
        action="store_true",
        help=(
            "Assign consecutive resSeq to mobile residues aligned to gaps in the reference "
            "(e.g. WT-only residues vs mut deletion, or N/C overhang)."
        ),
    )
    args = p.parse_args()

    rep = align_and_renumber_pdb(
        args.reference,
        args.mobile,
        args.output,
        chain_id=args.chain,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
        fill_unmapped_mobile=args.fill_unmapped_mobile,
    )
    print("Renumbering report:")
    print(f"  reference residues: {rep.n_reference_residues}")
    print(f"  mobile residues:    {rep.n_mobile_residues}")
    print(f"  mapped residues:    {rep.n_mapped}")
    print(f"  mismatches (aligned columns, different AA): {rep.n_mismatch_at_aligned_positions}")
    print(f"  mobile vs ref-gap columns (before fill): {rep.n_mobile_only_gaps_in_reference}")
    print(f"  filled (new numbers for those + overhangs): {rep.n_filled_unmapped_mobile}")
    print(f"  wrote: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
