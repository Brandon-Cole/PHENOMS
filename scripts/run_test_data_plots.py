#!/usr/bin/env python3
"""
Run PHENOMS on local test PDBs (phenoms/test_data/2_24_25_hbond_heatmaps/) and write plots.

Activate the **phenoms** conda env, ``cd`` to the repo root, then:

    python scripts/run_test_data_plots.py --resid-range 50 70

No separate package install needed for this script (repo root is on ``sys.path``).

Omit ``--sub-frames`` to use the full trajectory per PDB.

Output: ``$PHENOMS_OUTPUT_DIR/test_plots/<tag>/`` or ``./phenom_outputs/test_plots/<tag>/`` by default
(override with ``--output-dir``).
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
import json
from typing import Tuple


def main() -> int:
    os.environ.setdefault("MPLBACKEND", "Agg")

    parser = argparse.ArgumentParser(description="Generate plots from phenoms test_data PDBs.")
    parser.add_argument(
        "--sub-frames",
        type=int,
        default=None,
        metavar="N",
        help="Max frames per replicate. Omit for entire trajectory.",
    )
    parser.add_argument(
        "--resid-range",
        nargs=2,
        type=int,
        default=None,
        metavar=("START", "END"),
        help="Residue range (1-based, inclusive) to analyze. Default: whole protein.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=None,
        help="MDTraj parallel workers. Omit for all CPUs minus 2.",
    )
    parser.add_argument(
        "--no-rust",
        action="store_true",
        help="Use MDTraj Baker–Hubbard instead of the Rust extension (slower).",
    )
    parser.add_argument(
        "--export-pdb-only",
        action="store_true",
        help="Only export the B-factor PDB from the saved comparison CSV artifact (fast).",
    )
    parser.add_argument(
        "--comparison-csv",
        type=str,
        default=None,
        help="Path to comparison CSV artifact used for --export-pdb-only. Default: generated under out_dir.",
    )
    parser.add_argument(
        "--export-connectivity-graph-html",
        action="store_true",
        help="Export an interactive HTML residue connectivity graph (bond-level) next to other outputs.",
    )
    parser.add_argument(
        "--export-community-graphs-html",
        action="store_true",
        help="Export community-aware interactive connectivity graphs for set A, set B, and set difference.",
    )
    parser.add_argument(
        "--community-resolution",
        type=float,
        default=1.0,
        help="Resolution parameter for greedy modularity community detection (default: 1.0).",
    )
    parser.add_argument(
        "--export-bfactor-viewer-html",
        action="store_true",
        help="Export a standalone HTML viewer that colors the B-factor PDB by B-factor (no PyMOL required).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="Base directory for this run (default: $PHENOMS_OUTPUT_DIR/test_plots or ./phenom_outputs/test_plots).",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parent.parent
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    from phenoms.outputs import default_output_root

    data_dir = repo / "phenoms" / "test_data" / "2_24_25_hbond_heatmaps"
    resid_range = tuple(args.resid_range) if args.resid_range is not None else None
    out_tag = (
        f"resid_{resid_range[0]}_{resid_range[1]}" if resid_range is not None else "whole_protein"
    )
    base = Path(args.output_dir).expanduser().resolve() if args.output_dir else default_output_root() / "test_plots"
    out_dir = base / out_tag
    out_dir.mkdir(parents=True, exist_ok=True)
    # Remove legacy top-level artifacts from older script layouts to keep outputs organized.
    legacy_files = [
        "difference_lig_minus_no_lig_labeled.png",
        "difference_lig_minus_no_lig_unlabeled.png",
        "difference_autocorr_threshold_labeled.png",
        "difference_autocorr_threshold_unlabeled.png",
        "pca_two_groups.png",
        "community_graph_diff.html",
        "community_graph_set_a.html",
        "community_graph_set_b.html",
        "connectivity_graph_diff.html",
        "community_nodes_diff.csv",
        "community_nodes_set_a.csv",
        "community_nodes_set_b.csv",
        "community_summary_diff.csv",
        "community_summary_set_a.csv",
        "community_summary_set_b.csv",
        "comparison_lig_minus_no_lig.csv",
    ]
    for rel in legacy_files:
        p = out_dir / rel
        if p.exists() and p.is_file():
            p.unlink()
    comparison_csv_default = out_dir / "difference_analysis" / "raw_data" / "comparison_lig_minus_no_lig.csv"
    comparison_csv_path = Path(args.comparison_csv) if args.comparison_csv else comparison_csv_default

    def method_dirs(method_name: str) -> Tuple[Path, Path]:
        method_dir = out_dir / method_name
        raw_dir = method_dir / "raw_data"
        method_dir.mkdir(parents=True, exist_ok=True)
        raw_dir.mkdir(parents=True, exist_ok=True)
        return method_dir, raw_dir

    set_a_pdbs = [
        data_dir / "rep1_md_0_500_nolig_nojump_center.pdb",
        data_dir / "rep2_md_0_500_no_lig_nojump_center.pdb",
        data_dir / "rep3_md_0_500_no_lig_nojump_center.pdb",
    ]
    set_b_pdbs = sorted(data_dir.glob("rep*_md_0_300_55085_nojump_center.pdb"))
    missing = [p for p in set_a_pdbs + set_b_pdbs if not p.exists()]
    if missing:
        print("Missing test PDBs (copy from H_bond_work first):", file=sys.stderr)
        for p in missing:
            print(f"  {p}", file=sys.stderr)
        return 1

    if args.export_pdb_only:
        from phenoms.structure import write_pdb_bfactors
        from phenoms.viewer_bfactor_html import export_bfactor_colored_pdb_viewer_html
        bfactor_dir, _bfactor_raw = method_dirs("bfactor_analysis")

        if not comparison_csv_path.exists():
            print(
                f"Comparison CSV not found: {comparison_csv_path}. Run the script once without --export-pdb-only first.",
                file=sys.stderr,
            )
            return 1

        import pandas as pd

        comparison_df = pd.read_csv(comparison_csv_path)
        # Color B-factors on a static ligand structure (model/frame 0).
        reference_pdb = str(set_b_pdbs[0])
        pdb_out = bfactor_dir / "bfactors_lig_minus_no_lig_difference.pdb"
        write_pdb_bfactors(
            reference_pdb,
            comparison_df,
            str(pdb_out),
            value_column=None,
            model_index=0,
        )
        print(f"Wrote B-factor PDB: {pdb_out}")

        if args.export_bfactor_viewer_html:
            viewer_html = bfactor_dir / "bfactor_viewer_lig_minus_no_lig_difference.html"
            export_bfactor_colored_pdb_viewer_html(
                str(pdb_out),
                str(viewer_html),
                title="Ligand − no-ligand Δ (colored by B-factor)",
            )
            print(f"Wrote B-factor HTML viewer: {viewer_html}")

        return 0

    sub = None if args.sub_frames is None else args.sub_frames
    if sub == 0:
        sub = None

    use_rust = not args.no_rust
    if use_rust:
        try:
            import phenoms_hbond_rs  # noqa: F401
        except ImportError:
            print(
                "Rust extension not found. Run: maturin develop  (or use --no-rust)",
                file=sys.stderr,
            )
            return 1
        print("H-bond backend: Rust (phenoms_hbond_rs)")
    else:
        print("H-bond backend: MDTraj (--no-rust)")

    from phenoms.simulation import SimulationSet, default_n_jobs
    from phenoms.comparison_set import ComparisonSet
    from phenoms.plotting import suggest_difference_threshold_autocorr

    n_disp = args.n_jobs if args.n_jobs is not None else default_n_jobs()
    print("Set A (0_500 no lig):", len(set_a_pdbs), "replicates")
    print("Set B (0_300 55085):", len(set_b_pdbs), "replicates")
    print(f"sub_frames={'all' if sub is None else sub}, n_jobs={n_disp} (MDTraj fallback only)")

    try:
        sim_a = SimulationSet(
            [str(p) for p in set_a_pdbs],
            sub_frames=sub,
            resid_range=resid_range,
        ).run(n_jobs=args.n_jobs, use_rust=use_rust)
        sim_b = SimulationSet(
            [str(p) for p in set_b_pdbs],
            sub_frames=sub,
            resid_range=resid_range,
        ).run(n_jobs=args.n_jobs, use_rust=use_rust)
    except Exception as e:
        print(f"Run failed: {e}", file=sys.stderr)
        return 1

    if not sim_a.get_bond_labels_sorted() or not sim_b.get_bond_labels_sorted():
        print("No backbone N–O bonds found in one or both sets.", file=sys.stderr)
        return 1

    set_a_dir, set_a_raw = method_dirs("set_a_heatmaps")
    set_b_dir, set_b_raw = method_dirs("set_b_heatmaps")
    sim_a.plot_heatmaps(save_dir=str(set_a_dir))
    sim_b.plot_heatmaps(save_dir=str(set_b_dir))

    # Raw replicate-level artifacts per set.
    for i, hb in enumerate(sim_a.get_hbond_dfs(), start=1):
        hb.to_csv(set_a_raw / f"hbonds_rep{i}.csv", index=False)
    for i, pt in enumerate(sim_a.get_pivot_tables(), start=1):
        pt.to_csv(set_a_raw / f"pivot_rep{i}.csv")
    for i, hb in enumerate(sim_b.get_hbond_dfs(), start=1):
        hb.to_csv(set_b_raw / f"hbonds_rep{i}.csv", index=False)
    for i, pt in enumerate(sim_b.get_pivot_tables(), start=1):
        pt.to_csv(set_b_raw / f"pivot_rep{i}.csv")

    comp = ComparisonSet(
        sim_a,
        sim_b,
        label_a="no_lig_0_500",
        label_b="55085_set",
    )
    # Use ligand − no-ligand sign convention (positive values = higher protection in liganded set).
    comp.compare(flip_difference=True)
    comparison_heatmap_dir, comparison_heatmap_raw = method_dirs("comparison_aligned_heatmaps")
    comp.plot_heatmaps_both(save_dir=str(comparison_heatmap_dir))

    # Autocorr-derived effective delta threshold (used for label selection; also returned by plot meta).
    difference_dir, difference_raw = method_dirs("difference_analysis")
    comp_df = comp.get_comparison_df()
    comp_df.to_csv(comparison_csv_path, index=False)
    print(f"Wrote comparison CSV artifact: {comparison_csv_path}")
    autocorr_eff_t = suggest_difference_threshold_autocorr(comp_df, kappa=1.96)
    print(f"Autocorr effective diff threshold (for labeling): {autocorr_eff_t}")
    (difference_raw / "thresholds.json").write_text(
        json.dumps(
            {
                "manual_threshold": 0.2,
                "autocorr_effective_threshold": float(autocorr_eff_t),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    # Labeled + unlabeled variants (visual only). Underlying comparison/imputation remains unchanged.
    comp.plot_difference(
        diff_threshold_mode="manual",
        diff_threshold=0.2,
        impute_small_differences=False,
        show_labels=True,
        interpolate=True,
        save_path=str(difference_dir / "difference_lig_minus_no_lig_labeled.png"),
        title="Donor protection difference (ligand − no ligand) (labeled, manual threshold=0.2) (Entire protein)",
    )
    comp.plot_difference(
        diff_threshold_mode="manual",
        diff_threshold=0.2,
        impute_small_differences=False,
        show_labels=False,
        interpolate=True,
        save_path=str(difference_dir / "difference_lig_minus_no_lig_unlabeled.png"),
        title="Donor protection difference (ligand − no ligand) (unlabeled, manual threshold=0.2) (Entire protein)",
    )
    comp.plot_difference(
        diff_threshold_mode="autocorr",
        diff_threshold=autocorr_eff_t,
        impute_small_differences=False,
        show_labels=True,
        interpolate=True,
        save_path=str(difference_dir / "difference_autocorr_threshold_labeled.png"),
        title=f"Donor protection difference (ligand − no ligand) (labeled, autocorr threshold={autocorr_eff_t:.3f}) (Entire protein)",
    )
    comp.plot_difference(
        diff_threshold_mode="autocorr",
        diff_threshold=autocorr_eff_t,
        impute_small_differences=False,
        show_labels=False,
        interpolate=True,
        save_path=str(difference_dir / "difference_autocorr_threshold_unlabeled.png"),
        title=f"Donor protection difference (ligand − no ligand) (unlabeled, autocorr threshold={autocorr_eff_t:.3f}) (Entire protein)",
    )

    # Export visualization PDB: write the ligand − no-ligand Difference values into B-factors.
    # (Residues not present in the comparison mapping get B-factor=0.0.)
    # Color B-factors on a static ligand structure (model/frame 0).
    reference_pdb = str(set_b_pdbs[0])
    bfactor_dir, bfactor_raw = method_dirs("bfactor_analysis")
    pdb_out = bfactor_dir / "bfactors_lig_minus_no_lig_difference.pdb"
    from phenoms.structure import write_pdb_bfactors

    try:
        # Color B-factors on a static ligand structure (model/frame 0).
        write_pdb_bfactors(
            reference_pdb,
            comp_df,
            str(pdb_out),
            value_column=None,
            model_index=0,
        )
        print(f"Wrote B-factor PDB: {pdb_out}")
        comp_df.to_csv(bfactor_raw / "comparison_values_for_bfactors.csv", index=False)
    except Exception as e:
        print(f"Warning: failed to write B-factor PDB: {e}", file=sys.stderr)

    if args.export_connectivity_graph_html:
        connectivity_dir, connectivity_raw = method_dirs("connectivity_diff")
        graph_html = connectivity_dir / "connectivity_graph_diff.html"
        comp.export_connectivity_graph_html(
            str(graph_html),
            graph_mode="diff",
            top_k_edges=30,
            min_abs_edge_delta=None,
            directed=False,
        )
        comp_df.to_csv(connectivity_raw / "comparison_diff_input.csv", index=False)
        print(f"Wrote connectivity graph HTML: {graph_html}")

    if args.export_community_graphs_html:
        community_specs = [
            ("set_a", "community_set_a"),
            ("set_b", "community_set_b"),
            ("diff", "community_diff"),
        ]
        for mode, dirname in community_specs:
            community_dir, community_raw = method_dirs(dirname)
            html_path = community_dir / f"community_graph_{mode}.html"
            nodes_csv = community_raw / f"community_nodes_{mode}.csv"
            summary_csv = community_raw / f"community_summary_{mode}.csv"
            comp.export_connectivity_community_graph_html(
                str(html_path),
                graph_mode=mode,
                top_k_edges=40,
                min_abs_edge_delta=None,
                directed=False,
                resolution=float(args.community_resolution),
                community_nodes_csv_path=str(nodes_csv),
                community_summary_csv_path=str(summary_csv),
            )
            print(f"Wrote community graph HTML ({mode}): {html_path}")
            print(f"Wrote community node CSV ({mode}): {nodes_csv}")
            print(f"Wrote community summary CSV ({mode}): {summary_csv}")

    if args.export_bfactor_viewer_html:
        from phenoms.viewer_bfactor_html import export_bfactor_colored_pdb_viewer_html

        if pdb_out.exists():
            viewer_html = bfactor_dir / "bfactor_viewer_lig_minus_no_lig_difference.html"
            export_bfactor_colored_pdb_viewer_html(
                str(pdb_out),
                str(viewer_html),
                title="Ligand − no-ligand Δ (colored by B-factor)",
            )
            print(f"Wrote B-factor HTML viewer: {viewer_html}")

    import matplotlib.pyplot as plt

    plt.close("all")
    comp.run_pca(plot=True, title="PCA (no_lig_0_500 vs 55085_0_300)")
    pca_dir, pca_raw = method_dirs("pca_analysis")
    plt.savefig(pca_dir / "pca_two_groups.png", dpi=150, bbox_inches="tight")
    plt.close()
    comp_df.to_csv(pca_raw / "comparison_input_for_pca_context.csv", index=False)

    # Raw comparison-aligned heatmap matrix references.
    for i, pt in enumerate(comp.set_a.get_plot_pivot_tables(), start=1):
        pt.to_csv(comparison_heatmap_raw / f"set_a_pivot_rep{i}.csv")
    for i, pt in enumerate(comp.set_b.get_plot_pivot_tables(), start=1):
        pt.to_csv(comparison_heatmap_raw / f"set_b_pivot_rep{i}.csv")

    print("Wrote plots under:", out_dir)
    for p in sorted(out_dir.rglob("*.png")):
        print(" ", p.relative_to(out_dir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
