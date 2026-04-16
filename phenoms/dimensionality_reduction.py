"""
PCA, t-SNE, Isomap on aggregated replicate H-bond data.

Use with one simulation set (replicates as points, colored by replicate or custom labels)
or two sets (ComparisonSet: two groups).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from phenoms.analysis import extract_residue_numbers

def aggregate_replicate_data(pivot_tables_list):
    """
    One row per replicate: mean occupancy per bond across frames.
    Shape (n_replicates, n_bonds). Fill missing bonds with 0.
    """
    all_bonds = set()
    for pt in pivot_tables_list:
        all_bonds.update(pt.index.tolist())
    bond_order = sorted(all_bonds, key=extract_residue_numbers)
    n_rep = len(pivot_tables_list)
    n_bonds = len(bond_order)
    out = np.zeros((n_rep, n_bonds))
    for i, pt in enumerate(pivot_tables_list):
        mean_occ = (pt > 0).astype(float).mean(axis=1)
        for j, b in enumerate(bond_order):
            out[i, j] = mean_occ.get(b, 0.0)
    return out


def run_pca(aggregated_data, n_components=2):
    """Standardize and run PCA. Returns (scores, explained_variance_ratio)."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    X = StandardScaler().fit_transform(aggregated_data)
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(X)
    return scores, pca.explained_variance_ratio_


def run_tsne(aggregated_data, perplexity=2, random_state=42):
    """Standardize and run t-SNE. Returns 2D scores."""
    from sklearn.manifold import TSNE
    from sklearn.preprocessing import StandardScaler
    X = StandardScaler().fit_transform(aggregated_data)
    tsne = TSNE(n_components=2, perplexity=min(perplexity, max(1, aggregated_data.shape[0] - 1)), random_state=random_state)
    return tsne.fit_transform(X)


def run_isomap(aggregated_data, n_neighbors=5, n_components=2):
    """Standardize and run Isomap. Returns 2D scores."""
    from sklearn.manifold import Isomap
    from sklearn.preprocessing import StandardScaler
    X = StandardScaler().fit_transform(aggregated_data)
    n_n = min(n_neighbors, aggregated_data.shape[0] - 1)
    if n_n < 2:
        return np.zeros((aggregated_data.shape[0], 2))
    iso = Isomap(n_components=n_components, n_neighbors=n_n)
    return iso.fit_transform(X)


def plot_manifold(
    scores,
    group_labels,
    replicate_labels=None,
    title="PCA",
    palette=None,
    ax=None,
):
    """
    Scatter plot of 2D manifold (PCA/t-SNE/Isomap) with points colored by group.

    Parameters
    ----------
    scores : np.ndarray
        (n_replicates, 2).
    group_labels : list or array
        One label per replicate (e.g. True/False for ligand, or "no_lig"/"lig").
    replicate_labels : list or None
        Optional text label per point.
    title : str
    palette : dict or None
        Mapping group_label -> color.
    ax : matplotlib axes or None
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    scores = np.asarray(scores)
    group_labels = list(group_labels)
    if palette is None:
        uniq = sorted(set(group_labels), key=str)
        palette = {g: c for g, c in zip(uniq, plt.cm.Set1.colors)}
    for g in set(group_labels):
        mask = [gl == g for gl in group_labels]
        ax.scatter(scores[mask, 0], scores[mask, 1], label=str(g), color=palette.get(g, "gray"), s=100, edgecolor="black")
    ax.set_xlabel("Component 1", fontsize=14, weight="bold")
    ax.set_ylabel("Component 2", fontsize=14, weight="bold")
    ax.set_title(title, fontsize=16, weight="bold")
    ax.legend(title="Group", frameon=False)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    return ax


def run_manifold_suite(
    pivot_tables,
    group_labels,
    perplexity=2,
    n_neighbors=5,
    random_state=42,
    plot=True,
):
    """
    Run PCA, t-SNE, and Isomap on the same aggregated matrix (glassfrog-style workflow).

    Parameters
    ----------
    pivot_tables : list of pd.DataFrame
        Bond × frame pivots (one per replicate).
    group_labels : list
        One label per replicate (e.g. replicate name, or True/False for two conditions).
    perplexity, n_neighbors, random_state
        Passed to t-SNE / Isomap.
    plot : bool
        If True, show three scatter plots (PCA, t-SNE, Isomap).

    Returns
    -------
    dict with keys 'pca' (scores, var_ratio), 'tsne' (scores), 'isomap' (scores),
    and 'aggregated_data' (np.ndarray).
    """
    data = aggregate_replicate_data(pivot_tables)
    scores_pca, var_ratio = run_pca(data, n_components=2)
    scores_tsne = run_tsne(data, perplexity=perplexity, random_state=random_state)
    scores_iso = run_isomap(data, n_neighbors=n_neighbors, n_components=2)
    if plot:
        plot_manifold(scores_pca, group_labels, title="PCA")
        plot_manifold(scores_tsne, group_labels, title="t-SNE")
        plot_manifold(scores_iso, group_labels, title="Isomap")
    return {
        "pca": (scores_pca, var_ratio),
        "tsne": scores_tsne,
        "isomap": scores_iso,
        "aggregated_data": data,
    }
