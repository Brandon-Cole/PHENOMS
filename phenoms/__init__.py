"""
PHENOMS: Backbone H-bond network analysis for MD simulations (HDX-MS style).

Public API:
- SimulationSet: main class for single or replicate trajectories (region or whole protein).
- ComparisonSet: two-run differential analysis and exports.
- default_output_root: where CSVs/benchmarks go by default (``PHENOMS_OUTPUT_DIR`` or ./phenom_outputs).
- detect_hbonds / detect_hbonds_with_occupancy: direct Rust-backed H-bond detection (all bonds or backbone-only).
- load_trajectory, load_and_select_residues: I/O helpers.
- run_backbone_hbond_analysis: one-shot function for region-based analysis.
- cleanup (optional): align_and_renumber_pdb / renumber_many_to_reference for mismatched PDB numbering.
"""

from phenoms.io import load_trajectory, load_and_select_residues
from phenoms.outputs import default_output_root
from phenoms.hbond import (
    label_hbond,
    label_hbond_all,
    process_frame,
    process_frame_all,
    process_frames,
    process_frames_all,
    hbond_occupancy_table,
    export_hbond_occupancy_csv,
)
from phenoms.analysis import (
    extract_residue_numbers,
    create_pivot_table,
    calculate_bond_statistics,
)
from phenoms.plotting import (
    plot_heatmap,
    plot_heatmap_with_legend,
    plot_top_hbonds,
    plot_aggregated_heatmap,
    plot_bond_lifetimes_with_error_bars,
    plot_break_frequencies_with_error_bars,
)
from phenoms.simulation import SimulationSet, default_n_jobs
from phenoms.comparison_set import ComparisonSet
from phenoms.comparison import (
    compare_two_sets,
    aggregate_pivot_by_donor,
    differentially_protected,
    prefer_difference_clipped_column,
)
from phenoms.plotting import plot_difference, suggest_difference_threshold_autocorr
from phenoms.dimensionality_reduction import run_manifold_suite
try:
    from phenoms.structure import write_pdb_bfactors
except ImportError:
    write_pdb_bfactors = None

from phenoms.cleanup import (
    RenumberReport,
    ResidueKey,
    align_and_renumber_pdb,
    build_mobile_to_reference_map,
    ordered_residue_run,
    renumber_many_to_reference,
    renumber_pdb_file,
)

__all__ = [
    "SimulationSet",
    "default_n_jobs",
    "default_output_root",
    "ComparisonSet",
    "compare_two_sets",
    "aggregate_pivot_by_donor",
    "differentially_protected",
    "prefer_difference_clipped_column",
    "plot_difference",
    "suggest_difference_threshold_autocorr",
    "run_manifold_suite",
    "load_trajectory",
    "load_and_select_residues",
    "label_hbond",
    "label_hbond_all",
    "process_frame",
    "process_frame_all",
    "process_frames",
    "process_frames_all",
    "hbond_occupancy_table",
    "export_hbond_occupancy_csv",
    "extract_residue_numbers",
    "create_pivot_table",
    "calculate_bond_statistics",
    "plot_heatmap",
    "plot_heatmap_with_legend",
    "plot_top_hbonds",
    "plot_aggregated_heatmap",
    "plot_bond_lifetimes_with_error_bars",
    "plot_break_frequencies_with_error_bars",
    "run_backbone_hbond_analysis",
    "write_pdb_bfactors",
    "RenumberReport",
    "ResidueKey",
    "align_and_renumber_pdb",
    "build_mobile_to_reference_map",
    "ordered_residue_run",
    "renumber_many_to_reference",
    "renumber_pdb_file",
    "detect_hbonds",
    "detect_hbonds_with_occupancy",
]

__version__ = "0.1.0"


def run_backbone_hbond_analysis(
    pdb_files,
    resid_range=None,
    sub_frames=None,
    n_jobs=None,
    plot_heatmaps=True,
    bond_statistics_threshold=None,
    output_dir=None,
):
    """
    One-shot: load replicates, run N–O H-bond analysis, build pivot tables,
    optionally plot heatmaps and/or compute bond statistics (whole-protein style).

    Parameters
    ----------
    pdb_files : list of str
        Paths to PDB files (single or replicates).
    resid_range : tuple (int, int) or None
        Residue range; None = whole protein.
    sub_frames : int or None
        Frames per replicate; omit or ``None`` = full trajectory.
    n_jobs : int or None
        MDTraj workers; omit or ``None`` = all CPUs minus two.
    plot_heatmaps : bool
        If True, show heatmap per replicate.
    bond_statistics_threshold : float or None
        If set (e.g. 0.5), compute mean/std lifetime and break frequency across
        replicates and show bar plots (whole-protein workflow).
    output_dir : str, pathlib.Path, or None
        If set, CSV artifacts are written under ``output_dir/raw_data/`` after ``run()``.

    Returns
    -------
    sim : SimulationSet
        After run(); use sim.get_pivot_tables(), sim.get_hbond_dfs(), etc.
    """
    sim = SimulationSet(
        pdb_files=pdb_files,
        resid_range=resid_range,
        sub_frames=sub_frames,
        bond_statistics_threshold=bond_statistics_threshold,
        output_dir=output_dir,
    )
    sim.run(n_jobs=n_jobs)
    if plot_heatmaps:
        sim.plot_heatmaps()
    if bond_statistics_threshold is not None and sim.bond_statistics is not None:
        sim.plot_bond_statistics()
        sim.plot_aggregated_heatmap()
    return sim


def detect_hbonds(
    pdb_file,
    *,
    sub_frames=None,
    n_jobs=None,
    use_rust=True,
    backbone_only=False,
):
    """
    Detect hydrogen bonds from a trajectory PDB path with optional Rust acceleration.

    Parameters
    ----------
    pdb_file : str
        Input trajectory/topology path.
    sub_frames : int or None
        Number of frames to process; None uses all.
    n_jobs : int or None
        Parallel jobs for MDTraj fallback; defaults to :func:`default_n_jobs`.
    use_rust : bool
        If True, use Rust backend when available.
    backbone_only : bool
        If True, keep only backbone-relevant N-O labels. If False (default), return
        all detected donor/acceptor H-bonds.
    """
    if n_jobs is None:
        n_jobs = default_n_jobs()
    traj = load_trajectory(pdb_file)
    if backbone_only:
        return process_frames(
            traj,
            sub_frames=sub_frames,
            n_jobs=n_jobs,
            use_rust=use_rust,
        )
    return process_frames_all(
        traj,
        sub_frames=sub_frames,
        n_jobs=n_jobs,
        use_rust=use_rust,
    )


def detect_hbonds_with_occupancy(
    pdb_file,
    *,
    sub_frames=None,
    n_jobs=None,
    use_rust=True,
    backbone_only=False,
    output_csv_path=None,
):
    """
    Detect H-bonds and return both frame-level and occupancy summaries.
    Optionally writes occupancy summary CSV.
    """
    hbonds_df = detect_hbonds(
        pdb_file,
        sub_frames=sub_frames,
        n_jobs=n_jobs,
        use_rust=use_rust,
        backbone_only=backbone_only,
    )
    occupancy_df = hbond_occupancy_table(hbonds_df)
    if output_csv_path is not None:
        occupancy_df.to_csv(output_csv_path, index=False)
    return hbonds_df, occupancy_df
