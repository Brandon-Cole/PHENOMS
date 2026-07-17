"""
Click CLI for PHENOMS.

Python API remains primary; this wraps common prep / run / compare workflows.
"""

from __future__ import annotations

from pathlib import Path

import click

from phenoms import (
    SimulationSet,
    __version__,
    comparison_sets_from_dirs,
    default_output_root,
    prepare_set_from_dir,
)


def _parse_resid_range(values):
    if values is None:
        return None
    if len(values) != 2:
        raise click.BadParameter("resid-range needs exactly two integers: START END")
    return (int(values[0]), int(values[1]))


@click.group()
@click.version_option(__version__, prog_name="phenoms")
def main():
    """PHENOMS: MD hydrogen-bond / HDX-style analysis."""


@main.command("prep")
@click.option("--input-dir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--prepared-dir", required=True, type=click.Path(file_okay=False))
@click.option("--frame-dt-ps", default=1000.0, show_default=True, type=float)
@click.option("--start-ps", default=None, type=float)
@click.option("--end-ps", default=None, type=float)
@click.option("--no-imaging", is_flag=True, help="Skip best-effort PBC imaging.")
@click.option("--no-center", is_flag=True, help="Skip centering.")
@click.option("--no-fit", is_flag=True, help="Skip fitting/superposition.")
def prep_cmd(input_dir, prepared_dir, frame_dt_ps, start_ps, end_ps, no_imaging, no_center, no_fit):
    """Normalize GROMACS/OpenMM/AMBER replicate folders to multi-frame PDBs."""
    paths = prepare_set_from_dir(
        input_dir,
        prepared_dir,
        frame_dt_ps=frame_dt_ps,
        start_ps=start_ps,
        end_ps=end_ps,
        apply_imaging=not no_imaging,
        center=not no_center,
        fit=not no_fit,
    )
    click.echo(f"Prepared {len(paths)} replicate PDB(s) under {prepared_dir}")
    for p in paths:
        click.echo(f"  {p}")


@main.command("run")
@click.option("--pdb", "pdb_files", multiple=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--traj", "trajectories", multiple=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--topology", default=None, type=click.Path(exists=True, dir_okay=False))
@click.option("--topology-each", "topologies", multiple=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--output-dir", default=None, type=click.Path(file_okay=False))
@click.option("--sub-frames", default=None, type=int)
@click.option("--resid-range", nargs=2, type=int, default=None)
@click.option("--all-bonds", is_flag=True, help="Detect all H-bonds (default is backbone N-O only).")
@click.option("--bond-statistics-threshold", default=None, type=float)
@click.option("--n-jobs", default=None, type=int)
@click.option("--no-rust", is_flag=True, help="Force MDTraj fallback.")
@click.option("--qc", is_flag=True, help="Enable RMSD QC (fail-fast on nonconvergence).")
def run_cmd(
    pdb_files,
    trajectories,
    topology,
    topologies,
    output_dir,
    sub_frames,
    resid_range,
    all_bonds,
    bond_statistics_threshold,
    n_jobs,
    no_rust,
    qc,
):
    """Run H-bond analysis on PDB replicates or native trajectories."""
    if bool(pdb_files) == bool(trajectories):
        raise click.UsageError("Provide either --pdb (one or more) or --traj (one or more), not both/neither.")
    if trajectories and not topology and not topologies:
        raise click.UsageError("Native --traj inputs require --topology or --topology-each.")

    out = Path(output_dir) if output_dir else (default_output_root() / "cli_run")
    resid = _parse_resid_range(resid_range)

    if pdb_files:
        sim = SimulationSet(
            pdb_files=list(pdb_files),
            resid_range=resid,
            sub_frames=sub_frames,
            bond_statistics_threshold=bond_statistics_threshold,
            output_dir=out,
            backbone_only=not all_bonds,
        )
    else:
        sim = SimulationSet.from_trajectories(
            list(trajectories),
            topology=topology,
            topologies=list(topologies) if topologies else None,
            resid_range=resid,
            sub_frames=sub_frames,
            bond_statistics_threshold=bond_statistics_threshold,
            output_dir=out,
            backbone_only=not all_bonds,
        )

    sim.run(n_jobs=n_jobs, use_rust=not no_rust, qc=qc)
    click.echo(f"Wrote artifacts under {out}")
    click.echo(f"Bonds (union): {len(sim.get_bond_labels_sorted())}")
    click.echo(f"Mode: {'all-bonds' if all_bonds else 'backbone-only'}")


@main.command("compare")
@click.option("--dir-a", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--dir-b", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("--prepared-dir-a", required=True, type=click.Path(file_okay=False))
@click.option("--prepared-dir-b", required=True, type=click.Path(file_okay=False))
@click.option("--output-dir", default=None, type=click.Path(file_okay=False))
@click.option("--label-a", default="set_a", show_default=True)
@click.option("--label-b", default="set_b", show_default=True)
@click.option("--sub-frames", default=None, type=int)
@click.option("--resid-range", nargs=2, type=int, default=None)
@click.option("--frame-dt-ps", default=1000.0, show_default=True, type=float)
@click.option("--all-bonds", is_flag=True, help="Detect all H-bonds (default backbone N-O).")
@click.option("--n-jobs", default=None, type=int)
def compare_cmd(
    dir_a,
    dir_b,
    prepared_dir_a,
    prepared_dir_b,
    output_dir,
    label_a,
    label_b,
    sub_frames,
    resid_range,
    frame_dt_ps,
    all_bonds,
    n_jobs,
):
    """Prep two class directories, run analysis, and export comparison artifacts."""
    out = Path(output_dir) if output_dir else (default_output_root() / "cli_compare")
    resid = _parse_resid_range(resid_range)
    set_a, set_b, comp = comparison_sets_from_dirs(
        dir_a,
        dir_b,
        prepared_dir_a,
        prepared_dir_b,
        resid_range=resid,
        sub_frames=sub_frames,
        output_dir_a=out / label_a,
        output_dir_b=out / label_b,
        frame_dt_ps=frame_dt_ps,
        label_a=label_a,
        label_b=label_b,
        backbone_only=not all_bonds,
    )
    set_a.run(n_jobs=n_jobs)
    set_b.run(n_jobs=n_jobs)
    comp.compare()
    comp.export_comparison_artifacts(out / "comparison")
    click.echo(f"Wrote comparison artifacts under {out / 'comparison'}")


if __name__ == "__main__":
    main()
