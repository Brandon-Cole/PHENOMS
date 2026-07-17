"""
PHENOMS: Backbone H-bond network analysis for MD simulations (HDX-MS style).

High-level API:
- SimulationSet / ComparisonSet: replicate analysis and differential comparison
- simulation_set_from_dir / comparison_sets_from_dirs: prep + analyze from engine folders
- detect_hbonds / detect_hbonds_with_occupancy: direct detection helpers
- run_backbone_hbond_analysis: one-shot PDB workflow

Also see the ``phenoms`` CLI (``phenoms --help``) for prep / run / compare.
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
from phenoms.qc import (
    parse_mdp,
    check_mdp_key_consistency,
    assess_series_convergence,
    rmsd_convergence_report,
)
from phenoms.preprocess import (
    ReplicateInput,
    detect_replicate_input,
    discover_replicate_dirs,
    normalize_replicate_to_pdb,
    prepare_set_from_dir,
)

__all__ = [
    # High-level workflow
    "SimulationSet",
    "ComparisonSet",
    "default_n_jobs",
    "default_output_root",
    "run_backbone_hbond_analysis",
    "simulation_set_from_dir",
    "comparison_sets_from_dirs",
    "detect_hbonds",
    "detect_hbonds_with_occupancy",
    # Comparison helpers
    "compare_two_sets",
    "aggregate_pivot_by_donor",
    "differentially_protected",
    "prefer_difference_clipped_column",
    "plot_difference",
    "suggest_difference_threshold_autocorr",
    "run_manifold_suite",
    # I/O + detection
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
    # Plotting
    "plot_heatmap",
    "plot_heatmap_with_legend",
    "plot_top_hbonds",
    "plot_aggregated_heatmap",
    "plot_bond_lifetimes_with_error_bars",
    "plot_break_frequencies_with_error_bars",
    # Structure / cleanup / QC / preprocess
    "write_pdb_bfactors",
    "RenumberReport",
    "ResidueKey",
    "align_and_renumber_pdb",
    "build_mobile_to_reference_map",
    "ordered_residue_run",
    "renumber_many_to_reference",
    "renumber_pdb_file",
    "parse_mdp",
    "check_mdp_key_consistency",
    "assess_series_convergence",
    "rmsd_convergence_report",
    "ReplicateInput",
    "detect_replicate_input",
    "discover_replicate_dirs",
    "normalize_replicate_to_pdb",
    "prepare_set_from_dir",
]

__version__ = "0.2.0"


def run_backbone_hbond_analysis(
    pdb_files,
    resid_range=None,
    sub_frames=None,
    n_jobs=None,
    plot_heatmaps=True,
    bond_statistics_threshold=None,
    output_dir=None,
    *,
    backbone_only=True,
):
    """
    One-shot: load replicates, run H-bond analysis, build pivot tables,
    optionally plot heatmaps and/or compute bond statistics.

    Defaults to backbone N–O mode (``backbone_only=True``). Pass
    ``backbone_only=False`` for all-bond detection.
    """
    sim = SimulationSet(
        pdb_files=pdb_files,
        resid_range=resid_range,
        sub_frames=sub_frames,
        bond_statistics_threshold=bond_statistics_threshold,
        output_dir=output_dir,
        backbone_only=backbone_only,
    )
    sim.run(n_jobs=n_jobs)
    if plot_heatmaps:
        sim.plot_heatmaps()
    if bond_statistics_threshold is not None and sim.bond_statistics is not None:
        sim.plot_bond_statistics()
        sim.plot_aggregated_heatmap()
    return sim


def detect_hbonds(
    path,
    *,
    top=None,
    sub_frames=None,
    n_jobs=None,
    use_rust=True,
    backbone_only=True,
):
    """
    Detect hydrogen bonds from a PDB or native trajectory path.

    Parameters
    ----------
    path : str
        Trajectory / PDB path.
    top : str or None
        Topology for native MD formats (``.xtc``/``.nc``/…).
    sub_frames : int or None
        Number of frames to process; None uses all.
    n_jobs : int or None
        Parallel jobs for MDTraj fallback; defaults to :func:`default_n_jobs`.
    use_rust : bool
        If True, use Rust backend when available.
    backbone_only : bool
        If True (default), keep only backbone N–O labels. If False, return all
        detected donor/acceptor H-bonds.
    """
    if n_jobs is None:
        n_jobs = default_n_jobs()
    traj = load_trajectory(path, top=top)
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
    path,
    *,
    top=None,
    sub_frames=None,
    n_jobs=None,
    use_rust=True,
    backbone_only=True,
    output_csv_path=None,
):
    """
    Detect H-bonds and return both frame-level and occupancy summaries.
    Optionally writes occupancy summary CSV.
    """
    hbonds_df = detect_hbonds(
        path,
        top=top,
        sub_frames=sub_frames,
        n_jobs=n_jobs,
        use_rust=use_rust,
        backbone_only=backbone_only,
    )
    occupancy_df = hbond_occupancy_table(hbonds_df)
    if output_csv_path is not None:
        occupancy_df.to_csv(output_csv_path, index=False)
    return hbonds_df, occupancy_df


def simulation_set_from_dir(
    input_dir,
    prepared_dir,
    *,
    resid_range=None,
    sub_frames=None,
    bond_statistics_threshold=None,
    output_dir=None,
    frame_dt_ps=1000.0,
    start_ps=None,
    end_ps=None,
    apply_imaging=True,
    center=True,
    fit=True,
    backbone_only=True,
):
    """
    Build a SimulationSet from a replicate directory or set directory.

    Normalizes engine-native inputs (GROMACS/OpenMM/AMBER) to multi-frame PDBs,
    then returns a ready-to-``.run()`` SimulationSet. Backbone-only is the default.
    """
    pdb_files = prepare_set_from_dir(
        input_dir,
        prepared_dir,
        frame_dt_ps=frame_dt_ps,
        start_ps=start_ps,
        end_ps=end_ps,
        apply_imaging=apply_imaging,
        center=center,
        fit=fit,
    )
    return SimulationSet(
        pdb_files=pdb_files,
        resid_range=resid_range,
        sub_frames=sub_frames,
        bond_statistics_threshold=bond_statistics_threshold,
        output_dir=output_dir,
        backbone_only=backbone_only,
    )


def comparison_sets_from_dirs(
    dir_a,
    dir_b,
    prepared_dir_a,
    prepared_dir_b,
    *,
    resid_range=None,
    sub_frames=None,
    bond_statistics_threshold=None,
    output_dir_a=None,
    output_dir_b=None,
    frame_dt_ps=1000.0,
    start_ps=None,
    end_ps=None,
    apply_imaging=True,
    center=True,
    fit=True,
    label_a="set_a",
    label_b="set_b",
    backbone_only=True,
):
    """
    Build two SimulationSet objects from two class directories (e.g., WT and MUT).

    Returns (set_a, set_b, comparison_set).
    """
    set_a = simulation_set_from_dir(
        dir_a,
        prepared_dir_a,
        resid_range=resid_range,
        sub_frames=sub_frames,
        bond_statistics_threshold=bond_statistics_threshold,
        output_dir=output_dir_a,
        frame_dt_ps=frame_dt_ps,
        start_ps=start_ps,
        end_ps=end_ps,
        apply_imaging=apply_imaging,
        center=center,
        fit=fit,
        backbone_only=backbone_only,
    )
    set_b = simulation_set_from_dir(
        dir_b,
        prepared_dir_b,
        resid_range=resid_range,
        sub_frames=sub_frames,
        bond_statistics_threshold=bond_statistics_threshold,
        output_dir=output_dir_b,
        frame_dt_ps=frame_dt_ps,
        start_ps=start_ps,
        end_ps=end_ps,
        apply_imaging=apply_imaging,
        center=center,
        fit=fit,
        backbone_only=backbone_only,
    )
    return set_a, set_b, ComparisonSet(set_a, set_b, label_a=label_a, label_b=label_b)
