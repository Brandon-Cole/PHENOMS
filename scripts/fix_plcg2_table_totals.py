#!/usr/bin/env python3
"""
Regenerate PLCg2 benchmark PDF table including per-method H-bond totals.

Usage:
  python scripts/fix_plcg2_table_totals.py \
    --benchmark-dir phenom_outputs/benchmarks/plcg2_docker_4cpu_100f_all
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _fmt_float(v: float) -> str:
    return f"{v:.3f}" if pd.notna(v) else "-"


def _fmt_int_or_dash(v: float) -> str:
    return str(int(v)) if pd.notna(v) else "-"


def render_table(benchmark_dir: Path) -> Path:
    per_window = benchmark_dir / "plcg2_truncations_per_window.csv"
    if not per_window.exists():
        raise FileNotFoundError(f"Missing file: {per_window}")

    df = pd.read_csv(per_window)

    # Ensure compatibility with older outputs that do not have renamed columns yet.
    if "Rust_H_Bond_Total" not in df.columns:
        df["Rust_H_Bond_Total"] = df.get("rust_hits")
    if "MDTraj_H_Bond_Total" not in df.columns:
        df["MDTraj_H_Bond_Total"] = df.get("mdtraj_hits")
    if "MDAnalysis_H_Bond_Total" not in df.columns:
        df["MDAnalysis_H_Bond_Total"] = df.get("mdanalysis_hits")

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
            "Rust_H_Bond_Total",
            "MDTraj_H_Bond_Total",
            "MDAnalysis_H_Bond_Total",
        ]
    ].copy()

    for col in (
        "rust_seconds",
        "mdtraj_seconds",
        "mdanalysis_seconds",
        "speedup_mdtraj_over_rust_x",
        "speedup_mda_over_rust_x",
    ):
        display[col] = display[col].map(_fmt_float)

    for col in ("Rust_H_Bond_Total", "MDTraj_H_Bond_Total", "MDAnalysis_H_Bond_Total"):
        display[col] = display[col].map(_fmt_int_or_dash)

    rows = display.values.tolist()
    fig_h = max(3.2, 1.2 + 0.38 * len(rows))
    fig, ax = plt.subplots(figsize=(19, fig_h))
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
            "Rust_Total",
            "MDTraj_Total",
            "MDAnalysis_Total",
        ],
        cellLoc="left",
        loc="center",
        colWidths=[0.17, 0.07, 0.06, 0.08, 0.09, 0.10, 0.08, 0.08, 0.10, 0.10, 0.12],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10.5)
    tbl.scale(1.0, 1.45)
    for c in range(11):
        cell = tbl[(0, c)]
        cell.set_text_props(weight="bold")
        cell.set_facecolor("#EDEDED")

    out_pdf = benchmark_dir / "plcg2_truncations_table.pdf"
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return out_pdf


def main() -> int:
    p = argparse.ArgumentParser(description="Regenerate PLCg2 table PDF with H-bond total columns.")
    p.add_argument("--benchmark-dir", type=str, required=True, help="Benchmark output directory path.")
    args = p.parse_args()

    out = render_table(Path(args.benchmark_dir).expanduser().resolve())
    print(f"Wrote: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

