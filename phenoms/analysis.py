"""
Pivot tables, bond label sorting, and bond statistics (lifetime, break frequency).
"""

import numpy as np
import pandas as pd


def extract_residue_numbers(bond_label):
    """
    Parse 'Residue1 -- Residue2' into (res1, res2) as integers.

    Parameters
    ----------
    bond_label : str
        e.g. 'Ala1 -- Gly2'

    Returns
    -------
    tuple (int, int)
    """
    residues = bond_label.split(" -- ")
    first = int("".join(filter(str.isdigit, residues[0])))
    second = int("".join(filter(str.isdigit, residues[1])))
    return first, second


def create_pivot_table(hbonds_df, bond_labels_sorted):
    """
    Pivot table: rows = bond labels (consistent order), columns = frames, values = count (0/1+).

    Parameters
    ----------
    hbonds_df : pd.DataFrame
        Must have 'Bond Label' and 'Frame'.
    bond_labels_sorted : list of str
        Full set of bond labels in desired order.

    Returns
    -------
    pd.DataFrame
    """
    hbonds_pivot = hbonds_df.pivot_table(
        index="Bond Label",
        columns="Frame",
        aggfunc="size",
        fill_value=0,
    )
    hbonds_pivot = hbonds_pivot.reindex(bond_labels_sorted, fill_value=0)
    return hbonds_pivot


def _average_bond_lifetime(row):
    """Mean length of consecutive runs of 1s in a binary row (frames)."""
    ones_indices = np.where(row == 1)[0]
    if len(ones_indices) == 0:
        return 0.0
    consecutive_lengths = []
    current_length = 1
    for i in range(1, len(ones_indices)):
        if ones_indices[i] == ones_indices[i - 1] + 1:
            current_length += 1
        else:
            consecutive_lengths.append(current_length)
            current_length = 1
    consecutive_lengths.append(current_length)
    return float(np.mean(consecutive_lengths))


def _frequency_of_breaks(row):
    """Number of frames where bond is absent (0)."""
    return int(np.sum(row == 0))


def calculate_bond_statistics(hbonds_pivot, threshold):
    """
    For bonds with mean occupancy > threshold, compute average lifetime and break frequency.

    Parameters
    ----------
    hbonds_pivot : pd.DataFrame
        Rows = bond labels, columns = frames, values = 0/1 (or counts).
    threshold : float
        Keep only bonds with mean(row) > threshold.

    Returns
    -------
    average_bond_lifetimes : pd.Series
    break_frequencies : pd.Series
    """
    valid_bonds = hbonds_pivot.index[hbonds_pivot.mean(axis=1) > threshold]
    filtered_pivot = hbonds_pivot.loc[valid_bonds]
    if filtered_pivot.empty:
        return pd.Series(dtype=float), pd.Series(dtype=float)
    # Treat as binary: any value > 0 -> 1 for lifetime/break logic
    binary = (filtered_pivot > 0).astype(int)
    average_bond_lifetimes = binary.apply(_average_bond_lifetime, axis=1)
    break_frequencies = binary.apply(_frequency_of_breaks, axis=1)
    return average_bond_lifetimes, break_frequencies


def fluctuating_bonds(pivot_tables, quantile=0.9, metric="variance"):
    """
    Bonds that fluctuate most (high variance across frames). For single-set highlighting.

    Parameters
    ----------
    pivot_tables : list of pd.DataFrame
    quantile : float in (0, 1]
        Keep bonds with variance >= this quantile (e.g. 0.9 = top 10%).
    metric : str
        'variance' = variance of occupancy across frames (mean across replicates).

    Returns
    -------
    pd.Series
        Bond labels (index) with variance (or metric) value, sorted descending.
    """
    combined = pd.concat(pivot_tables, axis=1)
    binary = (combined > 0).astype(float)
    var_per_bond = binary.var(axis=1)
    if var_per_bond.empty:
        return pd.Series(dtype=float)
    thresh = var_per_bond.quantile(quantile)
    high = var_per_bond[var_per_bond >= thresh].sort_values(ascending=False)
    return high
