"""
Heatmaps, top-H-bond bar plot, aggregated heatmap, and bond statistics bar plots.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch


def plot_heatmap(pivot_table, file_name, save_path=None, figsize=(12, 8)):
    """
    H-bond occupation over time (bond labels vs frames). Optionally save to file.

    Parameters
    ----------
    pivot_table : pd.DataFrame
    file_name : str
        Used in title.
    save_path : str or None
        If set, save figure to this path and do not show.
    figsize : tuple
    """
    plt.figure(figsize=figsize)
    sns.heatmap(pivot_table, cmap="bwr", cbar=False)
    plt.title(f"H-Bond Occupation over Time - {file_name}", weight="bold")
    plt.xlabel("Frame (ns)", weight="bold", labelpad=10)
    plt.ylabel("Bond Label", weight="bold", labelpad=10)
    n_cols = pivot_table.shape[1]
    step = max(1, n_cols // 20)
    plt.xticks(
        ticks=np.arange(0, n_cols, step),
        labels=np.arange(0, n_cols, step),
        rotation=90,
        fontsize=8,
    )
    plt.yticks(rotation=0)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_heatmap_with_legend(pivot_table, title, save_path=None, figsize=(12, 8)):
    """
    Same as plot_heatmap but with a legend for Absent (0) / Present (1).
    """
    plt.figure(figsize=figsize)
    sns.heatmap(pivot_table, cmap="bwr", cbar=False)
    legend_handles = [
        Patch(color="blue", label="Absent (0)"),
        Patch(color="red", label="Present (1)"),
    ]
    plt.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=10,
        bbox_to_anchor=(1.15, 1),
        frameon=False,
    )
    plt.title(title, weight="bold")
    plt.xlabel("Frame (ns)", weight="bold", labelpad=10)
    plt.ylabel("Bond Label", weight="bold", labelpad=10)
    n_cols = pivot_table.shape[1]
    step = max(1, n_cols // 20)
    plt.xticks(
        ticks=np.arange(0, n_cols, step),
        labels=np.arange(0, n_cols, step),
        rotation=90,
        fontsize=8,
    )
    plt.yticks(rotation=0)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_top_hbonds(all_pivot_tables, threshold, sub_frames=None, save_path=None):
    """
    Bar chart of top H-bonds by average lifetime (presence fraction > threshold).
    """
    if not all_pivot_tables:
        return
    combined = pd.concat(all_pivot_tables, axis=1).fillna(0)
    total_frames = sub_frames or combined.shape[1]
    presence_frac = (combined > 0).sum(axis=1) / total_frames
    combined["Average Lifetime"] = combined.mean(axis=1)
    significant = combined.loc[presence_frac > threshold].sort_values(
        by="Average Lifetime", ascending=False
    )
    if significant.empty:
        return
    plt.figure(figsize=(10, 6))
    significant["Average Lifetime"].plot(kind="bar", color="skyblue")
    plt.title(f"Top Hydrogen Bonds (Threshold: {threshold})")
    plt.xlabel("Hydrogen Bond")
    plt.ylabel("Average Lifetime")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_bond_lifetimes_with_error_bars(
    mean_lifetimes,
    std_lifetimes,
    title="Average Lifetime of Hydrogen Bonds with Standard Deviation (Threshold = 50%)",
    save_path=None,
):
    """Bar plot of mean bond lifetimes with std error bars across replicates."""
    if mean_lifetimes.empty:
        return
    sorted_indices = mean_lifetimes.index
    plt.figure(figsize=(12, 6))
    mean_lifetimes[sorted_indices].plot(
        kind="bar",
        yerr=std_lifetimes.reindex(sorted_indices),
        capsize=4,
        color="skyblue",
    )
    plt.title(title, weight="bold")
    plt.xlabel("Bond Label", weight="bold")
    plt.ylabel("Average Lifetime in ns (frames)", weight="bold")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_break_frequencies_with_error_bars(
    mean_break_frequencies,
    std_break_frequencies,
    title="Frequency of Breaks of Hydrogen Bonds with Standard Deviation",
    save_path=None,
):
    """Bar plot of mean break frequency with std error bars across replicates."""
    if mean_break_frequencies.empty:
        return
    sorted_indices = mean_break_frequencies.index
    plt.figure(figsize=(12, 6))
    mean_break_frequencies[sorted_indices].plot(
        kind="bar",
        yerr=std_break_frequencies.reindex(sorted_indices),
        capsize=4,
        color="lightcoral",
    )
    plt.title(title)
    plt.xlabel("Bond Label")
    plt.ylabel("Frequency of Breaks")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def plot_aggregated_heatmap(
    all_pivot_tables,
    save_path=None,
    figsize=(12, 8),
    region_str="Entire protein",
):
    """
    Average presence of each bond across replicates (one row heatmap).
    """
    if not all_pivot_tables:
        return
    concatenated = pd.concat(all_pivot_tables, axis=1)
    average_presence = concatenated.mean(axis=1)
    sorted_labels = concatenated.index
    avg_df = average_presence.reindex(sorted_labels).reset_index()
    avg_df.columns = ["Bond Label", "Average Presence"]
    plt.figure(figsize=figsize)
    sns.heatmap(
        avg_df[["Average Presence"]].T,
        cmap="bwr",
        cbar=True,
        annot=True,
        fmt=".2f",
    )
    plt.title(f"Average Presence of Hydrogen Bonds Over All Replicates ({region_str})")
    plt.xlabel("Bond Label")
    plt.ylabel("Average Presence")
    plt.xticks(
        ticks=np.arange(avg_df.shape[0]) + 0.5,
        labels=avg_df["Bond Label"],
        rotation=45,
        ha="right",
    )
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def suggest_difference_threshold_autocorr(comparison_df, y_column="Difference_clipped", kappa=1.96):
    """
    Estimate a |difference| threshold from sequence structure (lag-1 autocorrelation).

    Treats the per-residue difference series (in table order) as roughly AR(1); the
    innovation standard deviation scales local noise. Returns
    ``clip(kappa * sigma_innovation, 0.05, 0.35)``.

    Use with ``plot_difference(..., diff_threshold_mode="autocorr")`` or call directly.

    Parameters
    ----------
    comparison_df : pd.DataFrame
        Must be ordered along sequence (e.g. from compare_two_sets).
    y_column : str
    kappa : float
        Multiplier on the estimated noise scale (default ~95% under Gaussian heuristic).

    Returns
    -------
    float
    """
    col = y_column if y_column in comparison_df.columns else "Difference"
    y = comparison_df[col].astype(float).values
    n = len(y)
    if n < 3:
        return 0.2
    y = y - np.mean(y)
    rho = np.corrcoef(y[:-1], y[1:])[0, 1]
    if not np.isfinite(rho):
        rho = 0.0
    rho = float(np.clip(rho, -0.99, 0.99))
    var_y = float(np.var(y, ddof=1)) if n > 1 else 0.01
    sigma_eps = np.sqrt(max(var_y * (1.0 - rho**2), 1e-8))
    return float(np.clip(kappa * sigma_eps, 0.05, 0.35))


def plot_difference(
    comparison_df,
    x_column="Residue Number",
    y_column="Difference_clipped",
    title="Change in Protection (Set B − Set A)",
    diff_threshold=0.2,
    diff_threshold_mode="manual",
    impute_small_differences=True,
    show_labels=True,
    interpolate=True,
    interpolation_points=900,
    autocorr_kappa=1.96,
    save_path=None,
    figsize=(12, 8),
    return_meta=False,
):
    """
    Line plot of per-donor difference vs residue.

    By default ``y_column="Difference_clipped"`` (values in [-1, 1]); use
    ``y_column="Difference"`` for raw donor deltas.

    Parameters
    ----------
    diff_threshold : float or None
        Used when ``diff_threshold_mode="manual"``: residues with |y| above this get
        text labels. When imputing, plotted y is set to 0 where |y| < threshold.
        Default is ``0.2`` (20% delta).
    diff_threshold_mode : str
        ``"manual"`` — use ``diff_threshold``.
        ``"autocorr"`` — threshold from :func:`suggest_difference_threshold_autocorr`
        (``diff_threshold`` ignored except if autocorr fails).
    impute_small_differences : bool
        If True (default), plot values with |y| below the effective threshold as 0
        (smoothed / noise-suppressed line). If False, plot raw differences; labels
        still use the effective threshold when set.
    show_labels : bool
        If True (default), add red residue-number text for residues with
        |difference| >= effective threshold.
    interpolate : bool
        If True, draw a dashed smoothed line (visual only) using interpolation.
        This does not change the underlying data or thresholding.
    interpolation_points : int
        Number of x points to draw the interpolated line.
    autocorr_kappa : float
        Passed to autocorr threshold only.
    return_meta : bool
        If True, return dict with ``effective_threshold``, ``diff_threshold_mode``,
        ``impute_small_differences``.
    """
    if diff_threshold_mode not in ("manual", "autocorr"):
        raise ValueError("diff_threshold_mode must be 'manual' or 'autocorr'")

    y_col = y_column if y_column in comparison_df.columns else "Difference"

    if diff_threshold_mode == "autocorr":
        eff_t = suggest_difference_threshold_autocorr(
            comparison_df, y_column=y_col, kappa=autocorr_kappa
        )
    else:
        eff_t = diff_threshold

    plt.figure(figsize=figsize)
    x = comparison_df[x_column].astype(float)
    y_raw = comparison_df[y_col].astype(float)

    if impute_small_differences and eff_t is not None:
        y_plot = y_raw.where(y_raw.abs() >= eff_t, 0.0)
    else:
        y_plot = y_raw

    # Points + (optional) interpolated smoothed line.
    # The underlying curve uses a dashed style so labeled threshold plots are visually distinct.
    plt.plot(x, y_plot, marker="o", color="skyblue", linestyle="--", linewidth=2, markersize=4)

    if interpolate and eff_t is not None and len(x) >= 3:
        # Interpolation is visualization-only.
        x_fine = np.linspace(float(x.min()), float(x.max()), int(interpolation_points))
        y_fine = np.interp(x_fine, x.to_numpy(), y_plot.to_numpy())
        # Interpolation is visualization-only; use a dotted line for clarity.
        plt.plot(x_fine, y_fine, color="skyblue", linestyle=":", linewidth=2)

    if show_labels and eff_t is not None and "Donor Residue" in comparison_df.columns:
        mask = (y_raw > eff_t) | (y_raw < -eff_t)
        high = comparison_df.loc[mask]
        for idx, row in high.iterrows():
            yr = float(row[y_col])
            yp = float(y_plot.loc[idx])
            plt.text(
                float(row[x_column]),
                yp + 0.02 if yr > 0 else yp - 0.02,
                row["Donor Residue"],
                ha="left",
                va="bottom" if yr > 0 else "top",
                fontsize=9,
                color="red",
            )
    plt.xlabel("Donor Residue Number", fontsize=16, fontweight="bold")
    plt.ylabel(
        "Δ protection (clipped ±1)" if y_col == "Difference_clipped" else "Difference",
        fontsize=16,
        fontweight="bold",
    )
    plt.title(title, fontweight="bold")
    plt.grid(False)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()

    meta = {
        "effective_threshold": eff_t,
        "diff_threshold_mode": diff_threshold_mode,
        "impute_small_differences": impute_small_differences,
    }
    if return_meta:
        return meta
    return None
