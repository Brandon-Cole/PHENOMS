"""
SimulationSet: single or replicate MD runs, region or whole-protein H-bond analysis.
"""

import json
import os
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def default_n_jobs():
    """Parallel workers for MDTraj frame loop: all CPUs minus two (minimum 1)."""
    return max(1, (os.cpu_count() or 4) - 2)

from phenoms.io import load_and_select_residues
from phenoms.hbond import hbond_occupancy_table, process_frames
from phenoms.analysis import extract_residue_numbers, create_pivot_table, calculate_bond_statistics, fluctuating_bonds
from phenoms.plotting import (
    plot_heatmap,
    plot_heatmap_with_legend,
    plot_top_hbonds,
    plot_aggregated_heatmap,
    plot_bond_lifetimes_with_error_bars,
    plot_break_frequencies_with_error_bars,
    plot_difference,
)


class SimulationSet:
    """
    Single or replicate MD simulations with backbone N–O H-bond analysis.

    - resid_range=None: whole protein (no plot filtering).
    - resid_range=(start, end): filter bonds for heatmaps and manifold (PCA/t-SNE/Isomap),
      while running the underlying H-bond analysis on the entire protein.
    - bond_statistics_threshold: if set, compute per-replicate bond lifetime and
      break frequency, then mean ± std across replicates (whole-protein style).
    """

    def __init__(
        self,
        pdb_files,
        resid_range=None,
        sub_frames=None,
        bond_statistics_threshold=None,
        output_dir=None,
    ):
        """
        Parameters
        ----------
        pdb_files : list of str or str
            One or more PDB paths. If single path, pass [path].
        resid_range : tuple (int, int) or None
            Residue range for analysis; None = whole protein.
        sub_frames : int or None
            Number of frames to analyze per replicate. Omit or ``None`` = entire trajectory.
        bond_statistics_threshold : float or None
            If set, compute bond statistics (lifetime, break frequency) and
            mean ± std across replicates; used for bar plots and aggregated heatmap.
        output_dir : str, pathlib.Path, or None
            If set, after a successful :meth:`run`, per-replicate CSVs and a manifest are
            written under ``output_dir/raw_data/`` (see :meth:`export_run_artifacts`).
        """
        if isinstance(pdb_files, str):
            pdb_files = [pdb_files]
        self.pdb_files = list(pdb_files)
        self.resid_range = resid_range
        self.sub_frames = sub_frames
        self.bond_statistics_threshold = bond_statistics_threshold
        self.output_dir = Path(output_dir).expanduser().resolve() if output_dir is not None else None

        self._hbond_dfs = []
        self._pivot_tables = []
        self._bond_labels_sorted = []
        self._bond_statistics = None

    def run(self, n_jobs=None, use_rust=True):
        """
        Load all trajectories, run Baker–Hubbard (N–O only), build pivot tables.
        If bond_statistics_threshold was set, compute per-replicate stats and mean ± std.

        Parameters
        ----------
        n_jobs : int or None
            MDTraj parallel workers per replicate (Rust path ignores this). ``None`` =
            :func:`default_n_jobs` (all CPUs minus two).
        use_rust : bool
            If True, use Rust extension when available. Set False if Rust returns no bonds
            for your topology (MDTraj fallback).
        """
        if n_jobs is None:
            n_jobs = default_n_jobs()
        all_hbonds_dfs = []
        all_bond_labels = set()

        for pdb_file in tqdm(self.pdb_files, desc="Replicates", unit="replicate"):
            trajectory = load_and_select_residues(pdb_file, resid_range=None)
            hbonds_df = process_frames(
                trajectory,
                sub_frames=self.sub_frames,
                n_jobs=n_jobs,
                use_rust=use_rust,
            )
            all_bond_labels.update(hbonds_df["Bond Label"].unique())
            all_hbonds_dfs.append(hbonds_df)

        bond_labels_sorted = sorted(
            all_bond_labels,
            key=lambda x: extract_residue_numbers(x),
        )
        self._bond_labels_sorted = bond_labels_sorted

        all_pivot_tables = []
        for hbonds_df in all_hbonds_dfs:
            pivot = create_pivot_table(hbonds_df, bond_labels_sorted)
            all_pivot_tables.append(pivot)

        self._hbond_dfs = all_hbonds_dfs
        self._pivot_tables = all_pivot_tables

        if self.bond_statistics_threshold is not None:
            self._compute_bond_statistics()

        if self.output_dir is not None:
            self.export_run_artifacts(self.output_dir)

        return self

    def export_run_artifacts(self, output_dir):
        """
        Write per-replicate H-bond tables, Polars-backed occupancy summaries, pivots,
        and a small ``manifest.json`` under ``output_dir/raw_data/``.

        Safe to call again after :meth:`run` with a different path.
        """
        if not self._hbond_dfs:
            raise ValueError("No data. Run .run() first.")
        output_dir = Path(output_dir).expanduser().resolve()
        raw = output_dir / "raw_data"
        raw.mkdir(parents=True, exist_ok=True)
        for i, df in enumerate(self._hbond_dfs):
            stem = Path(self.pdb_files[i]).stem
            df.to_csv(raw / f"{stem}_hbonds.csv", index=False)
            hbond_occupancy_table(df).to_csv(raw / f"{stem}_occupancy.csv", index=False)
        for i, pivot in enumerate(self._pivot_tables):
            stem = Path(self.pdb_files[i]).stem
            pivot.to_csv(raw / f"{stem}_pivot.csv")
        manifest = {
            "pdb_files": [str(p) for p in self.pdb_files],
            "resid_range": self.resid_range,
            "sub_frames": self.sub_frames,
            "bond_statistics_threshold": self.bond_statistics_threshold,
        }
        (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    def _compute_bond_statistics(self):
        """Fill self._bond_statistics with mean/std lifetime and break frequency."""
        threshold = self.bond_statistics_threshold
        all_lifetimes = []
        all_breaks = []
        for pivot in self._pivot_tables:
            life, breaks = calculate_bond_statistics(pivot, threshold)
            all_lifetimes.append(life)
            all_breaks.append(breaks)
        lifetimes_df = pd.DataFrame(all_lifetimes)
        breaks_df = pd.DataFrame(all_breaks)
        self._bond_statistics = {
            "mean_lifetimes": lifetimes_df.mean(axis=0),
            "std_lifetimes": lifetimes_df.std(axis=0),
            "mean_break_frequencies": breaks_df.mean(axis=0),
            "std_break_frequencies": breaks_df.std(axis=0),
        }

    def get_hbond_dfs(self):
        """Return list of per-replicate H-bond DataFrames (Bond Label, Frame, ...)."""
        return self._hbond_dfs

    def get_pivot_tables(self):
        """Return list of per-replicate pivot tables (bond label x frame)."""
        return self._pivot_tables

    def get_bond_labels_sorted(self):
        """Return sorted list of all bond labels (union across replicates)."""
        return self._bond_labels_sorted

    def _plot_region_str(self) -> str:
        if self.resid_range is None:
            return "Entire protein"
        start, end = self.resid_range
        return f"Residues {start}-{end}"

    def _filter_bond_labels_for_plot(self, bond_labels):
        """
        Filter bond labels to those fully contained within `self.resid_range`.

        A bond label is included iff BOTH donor and acceptor residue numbers are
        within the range. Donor residue number is the first integer in the bond label.
        """
        if self.resid_range is None:
            return list(bond_labels)
        start, end = self.resid_range
        out = []
        for b in bond_labels:
            d, a = extract_residue_numbers(b)
            if start <= d <= end and start <= a <= end:
                out.append(b)
        return out

    def get_plot_bond_labels_sorted(self):
        """Bond labels (sorted) to use for heatmaps and manifold plots."""
        if not self._bond_labels_sorted:
            return []
        return self._filter_bond_labels_for_plot(self._bond_labels_sorted)

    def get_plot_pivot_tables(self):
        """Reindex each replicate pivot table to the plot bond list."""
        plot_labels = self.get_plot_bond_labels_sorted()
        return [pt.reindex(plot_labels, fill_value=0) for pt in self._pivot_tables]

    @property
    def bond_statistics(self):
        """None or dict with mean_lifetimes, std_lifetimes, mean_break_frequencies, std_break_frequencies."""
        return self._bond_statistics

    def plot_heatmaps(self, save_dir=None, use_legend=False):
        """
        Plot one heatmap per replicate (whole protein or region, depending on resid_range).
        Bond rows = union of bonds found in this set's replicates. For comparison of two sets
        with aligned bond lists, use ComparisonSet.plot_heatmaps_both() instead.
        """
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        plot_pivots = self.get_plot_pivot_tables()
        region_str = self._plot_region_str()
        for i, pivot in enumerate(plot_pivots):
            pdb_path = self.pdb_files[i] if i < len(self.pdb_files) else None
            stem = Path(pdb_path).stem if pdb_path else f"Replicate_{i + 1}"
            title_name = f"{stem} ({region_str})"
            path = f"{save_dir}/{stem}.png" if save_dir else None
            if use_legend:
                plot_heatmap_with_legend(
                    pivot,
                    f"H-Bond Occupation over Time - {title_name}",
                    save_path=path,
                )
            else:
                plot_heatmap(pivot, title_name, save_path=path)

    def plot_top_hbonds(self, threshold=0.5, save_path=None):
        """Bar plot of top H-bonds by average lifetime (presence fraction > threshold)."""
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        plot_top_hbonds(
            self._pivot_tables,
            threshold=threshold,
            sub_frames=self.sub_frames,
            save_path=save_path,
        )

    def plot_bond_statistics(
        self,
        lifetime_title="Average Lifetime of Hydrogen Bonds with Standard Deviation (Threshold = 50%)",
        break_title="Frequency of Breaks of Hydrogen Bonds with Standard Deviation",
        save_lifetime_path=None,
        save_break_path=None,
    ):
        """Bar plots for mean lifetime and break frequency with error bars (requires bond_statistics_threshold)."""
        if self._bond_statistics is None:
            raise ValueError("Bond statistics not computed. Set bond_statistics_threshold and run .run().")
        plot_bond_lifetimes_with_error_bars(
            self._bond_statistics["mean_lifetimes"],
            self._bond_statistics["std_lifetimes"],
            title=lifetime_title,
            save_path=save_lifetime_path,
        )
        plot_break_frequencies_with_error_bars(
            self._bond_statistics["mean_break_frequencies"],
            self._bond_statistics["std_break_frequencies"],
            title=break_title,
            save_path=save_break_path,
        )

    def plot_aggregated_heatmap(self, save_path=None):
        """Average presence of each bond across all replicates."""
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        plot_aggregated_heatmap(
            self.get_plot_pivot_tables(),
            save_path=save_path,
            region_str=self._plot_region_str(),
        )

    def get_fluctuating_bonds(self, quantile=0.9):
        """Bonds with highest variance (fluctuation) across frames. For single-set highlighting."""
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        return fluctuating_bonds(self._pivot_tables, quantile=quantile)

    def _default_replicate_labels(self):
        return [
            self.pdb_files[i] if i < len(self.pdb_files) else f"rep{i + 1}"
            for i in range(len(self._pivot_tables))
        ]

    def run_pca(self, group_labels=None, plot=True, title="PCA"):
        """PCA on this set’s replicates (one point per replicate). Same bond space as two-set case."""
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_pca as _run_pca, plot_manifold
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        data = aggregate_replicate_data(self.get_plot_pivot_tables())
        scores, var_ratio = _run_pca(data, n_components=2)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self._plot_region_str()})")
        return scores, var_ratio

    def run_tsne(self, group_labels=None, perplexity=2, plot=True, title="t-SNE", random_state=42):
        """t-SNE on this set’s replicates."""
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_tsne, plot_manifold
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        data = aggregate_replicate_data(self.get_plot_pivot_tables())
        scores = run_tsne(data, perplexity=perplexity, random_state=random_state)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self._plot_region_str()})")
        return scores

    def run_isomap(self, group_labels=None, n_neighbors=5, plot=True, title="Isomap"):
        """Isomap on this set’s replicates."""
        from phenoms.dimensionality_reduction import aggregate_replicate_data, run_isomap, plot_manifold
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        data = aggregate_replicate_data(self.get_plot_pivot_tables())
        scores = run_isomap(data, n_neighbors=n_neighbors)
        if plot:
            plot_manifold(scores, group_labels, title=f"{title} ({self._plot_region_str()})")
        return scores

    def run_manifold_suite(self, group_labels=None, perplexity=2, n_neighbors=5, random_state=42, plot=True):
        """
        PCA, t-SNE, and Isomap in one call (single-set replicates as points).
        group_labels: e.g. one label per replicate; default = PDB filenames.
        """
        from phenoms.dimensionality_reduction import run_manifold_suite as _suite
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        group_labels = group_labels or self._default_replicate_labels()
        return _suite(
            self.get_plot_pivot_tables(),
            group_labels,
            perplexity=perplexity,
            n_neighbors=n_neighbors,
            random_state=random_state,
            plot=plot,
        )

    def write_structure_bfactors(self, pdb_path, output_path, metric="variance"):
        """
        Write PDB with B-factors = per-residue metric (for PyMOL/Chimera).
        metric 'variance': highlight fluctuating residues; 'mean': mean occupancy.
        Requires biopython. Single-set only (no comparison).
        """
        from phenoms.comparison import per_donor_metric_from_pivots
        from phenoms.structure import write_pdb_bfactors
        if not self._pivot_tables:
            raise ValueError("No data. Run .run() first.")
        donor_df = per_donor_metric_from_pivots(
            self._pivot_tables, metric=metric, donor_aggregation="sum"
        )
        write_pdb_bfactors(pdb_path, donor_df, output_path, value_column="Average")
