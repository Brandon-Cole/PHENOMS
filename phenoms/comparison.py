"""
Two-set comparison: differential protection (occupation difference) with optional threshold.
Donor-residue aggregation and comparison DataFrame for plotting / PDB B-factors.
"""

import pandas as pd
import numpy as np


def bond_label_to_donor_residue(bond_label):
    """'ALA1 -- GLY2' -> 'ALA1' (first residue = donor convention)."""
    return bond_label.split(" -- ")[0].strip()


def aggregate_pivot_by_donor(pivot_tables, total_frames=None, donor_aggregation="sum"):
    """
    Aggregate bond-level occupancy to per-donor-residue average protection.

    For each donor residue, average the mean occupancy of all bonds that have
    that donor (first residue in bond label). Uses mean across replicates then
    mean across bonds per donor.

    Parameters
    ----------
    pivot_tables : list of pd.DataFrame
        Each: index = Bond Label, columns = Frame, values = 0/1 or count.
    total_frames : int or None
        If None, use number of columns (frames) from first pivot.

    donor_aggregation : str
        How to combine multiple bond labels that share the same donor residue.
        - "mean": average across bonds for that donor (current behavior historically)
        - "sum": sum across bonds for that donor (matches your old notebook-style weighting)
        - "max": maximum across bonds for that donor (donor protected if any bond is present)

    Returns
    -------
    pd.DataFrame
        Columns: Donor Residue, Average (mean protection 0-1).
    """
    if not pivot_tables:
        return pd.DataFrame(columns=["Donor Residue", "Average"])

    combined = pd.concat(pivot_tables, axis=1)
    mean_occ = (combined > 0).astype(float).mean(axis=1)
    donor_series = mean_occ.index.map(bond_label_to_donor_residue)
    df = pd.DataFrame(
        {
            "Bond Label": mean_occ.index,
            "Mean Occupancy": mean_occ.values,
            "Donor Residue": donor_series,
        }
    )

    donor_aggregation = donor_aggregation.lower()
    if donor_aggregation == "mean":
        donor_metric = df.groupby("Donor Residue", as_index=False)["Mean Occupancy"].mean()
    elif donor_aggregation == "sum":
        donor_metric = df.groupby("Donor Residue", as_index=False)["Mean Occupancy"].sum()
    elif donor_aggregation == "max":
        donor_metric = df.groupby("Donor Residue", as_index=False)["Mean Occupancy"].max()
    else:
        raise ValueError("donor_aggregation must be one of: 'mean', 'sum', 'max'")

    donor_metric = donor_metric.rename(columns={"Mean Occupancy": "Average"})
    return donor_metric


def compare_two_sets(
    pivot_tables_a,
    pivot_tables_b,
    label_a="set_a",
    label_b="set_b",
    impute_threshold=0.4,
    flip_difference=False,
    donor_aggregation="sum",
):
    """
    Compare two simulation sets: per-donor average protection and difference.

    Parameters
    ----------
    pivot_tables_a, pivot_tables_b : list of pd.DataFrame
        Pivot tables (bond x frame) for each set.
    label_a, label_b : str
        Column names for the two averages.
    impute_threshold : float or None
        Where both average_a and average_b <= impute_threshold, set Difference to 0.
        Set to None to skip imputation.
    flip_difference : bool
        If True, Difference = average_b - average_a (e.g. Ligand - No Ligand).

    Returns
    -------
    pd.DataFrame
        Donor Residue, {label_a}, {label_b}, Difference (raw), Difference_clipped
        (``clip`` to [-1, 1] for visualization / PDB), Residue Number. Sorted by residue.
    """
    donor_a = aggregate_pivot_by_donor(
        pivot_tables_a, donor_aggregation=donor_aggregation
    )
    donor_b = aggregate_pivot_by_donor(
        pivot_tables_b, donor_aggregation=donor_aggregation
    )
    donor_a = donor_a.rename(columns={"Average": label_a})
    donor_b = donor_b.rename(columns={"Average": label_b})
    merged = donor_a.merge(donor_b, on="Donor Residue", how="outer").fillna(0)
    merged["Difference"] = merged[label_a] - merged[label_b]
    if flip_difference:
        merged["Difference"] = -merged["Difference"]
    if impute_threshold is not None:
        low_both = (merged[label_a] <= impute_threshold) & (merged[label_b] <= impute_threshold)
        merged.loc[low_both, "Difference"] = 0
    merged["Difference_clipped"] = merged["Difference"].astype(float).clip(-1.0, 1.0)
    merged["Residue Number"] = merged["Donor Residue"].str.extract(r"(\d+)", expand=False).astype(int)
    merged = merged.sort_values("Residue Number").reset_index(drop=True)
    out_cols = ["Donor Residue", label_a, label_b, "Difference", "Difference_clipped", "Residue Number"]
    return merged[out_cols]


def prefer_difference_clipped_column(comparison_df: pd.DataFrame) -> str:
    """Return ``Difference_clipped`` if present (new CSVs), else ``Difference``."""
    if "Difference_clipped" in comparison_df.columns:
        return "Difference_clipped"
    return "Difference"


def differentially_protected(
    comparison_df,
    threshold=0.2,
    *,
    difference_column: str = "Difference_clipped",
):
    """Filter to rows where abs(difference_column) > threshold (default: clipped difference)."""
    col = difference_column
    if col not in comparison_df.columns:
        col = "Difference"
    return comparison_df[(comparison_df[col] > threshold) | (comparison_df[col] < -threshold)]


def per_donor_metric_from_pivots(
    pivot_tables, metric="variance", donor_aggregation="sum"
):
    """
    Per-donor metric for single-set structure labeling (fluctuation or mean occupancy).
    metric 'variance': mean across bonds (with that donor) of variance across frames; high = fluctuating.
    metric 'mean': same as aggregate_pivot_by_donor (mean occupancy).
    """
    if metric == "mean":
        return aggregate_pivot_by_donor(
            pivot_tables, donor_aggregation=donor_aggregation
        )
    combined = pd.concat(pivot_tables, axis=1)
    binary = (combined > 0).astype(float)
    var_per_bond = binary.var(axis=1)
    donor_series = var_per_bond.index.map(bond_label_to_donor_residue)
    df = pd.DataFrame({"Bond Label": var_per_bond.index, "Variance": var_per_bond.values, "Donor Residue": donor_series})
    donor_avg = df.groupby("Donor Residue", as_index=False)["Variance"].mean()
    donor_avg = donor_avg.rename(columns={"Variance": "Average"})
    donor_avg["Residue Number"] = donor_avg["Donor Residue"].str.extract(r"(\d+)", expand=False).astype(int)
    return donor_avg
