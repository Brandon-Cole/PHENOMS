import pandas as pd
import numpy as np

from phenoms.connectivity import (
    bond_occupancy_from_pivot_tables,
    bond_delta_from_pivot_tables,
    build_residue_graph_from_bond_deltas,
    build_residue_graph_with_communities,
    community_tables,
    export_residue_graph_html,
)


def _pivot_for_labels(labels, frames, data_rows):
    # data_rows: list of lists length=len(labels); each inner list length=len(frames)
    return pd.DataFrame(data_rows, index=labels, columns=frames)


def test_bond_occupancy_mean_across_reps():
    bond_labels = ["Ala1 -- Gly2", "Ala1 -- Ser3"]
    frames = [0, 1, 2]

    # Rep 1: bond1 present in frames 0,1; bond2 present in frame 2 only
    pt1 = _pivot_for_labels(bond_labels, frames, [[1, 1, 0], [0, 0, 1]])
    # Rep 2: bond1 present in frames 0 only; bond2 present in frames 0,1,2
    pt2 = _pivot_for_labels(bond_labels, frames, [[1, 0, 0], [1, 1, 1]])

    occ = bond_occupancy_from_pivot_tables([pt1, pt2], bond_labels)
    # bond1 rep occ: rep1=2/3, rep2=1/3 => mean=1/2
    assert np.isclose(occ.loc["Ala1 -- Gly2"], 0.5)
    # bond2 rep occ: rep1=1/3, rep2=1 => mean=2/3
    assert np.isclose(occ.loc["Ala1 -- Ser3"], 2 / 3)


def test_bond_delta_flip_difference_sign():
    bond_labels = ["Ala1 -- Gly2", "Ala1 -- Ser3"]
    frames = [0, 1, 2]

    # Set A: bond1 present in 0,1; bond2 present in 2 only
    pt_a = _pivot_for_labels(bond_labels, frames, [[1, 1, 0], [0, 0, 1]])
    # Set B: swap occupancies
    pt_b = _pivot_for_labels(bond_labels, frames, [[0, 0, 1], [1, 1, 0]])

    delta_no_flip = bond_delta_from_pivot_tables(
        [pt_a],
        [pt_b],
        bond_labels,
        flip_difference=False,
        impute_threshold=None,
    )
    # Without flip: delta = occ_a - occ_b
    assert delta_no_flip.loc["Ala1 -- Gly2"] > 0
    assert delta_no_flip.loc["Ala1 -- Ser3"] < 0

    delta_flip = bond_delta_from_pivot_tables(
        [pt_a],
        [pt_b],
        bond_labels,
        flip_difference=True,
        impute_threshold=None,
    )
    # With flip: sign inverted
    assert np.isclose(delta_flip.loc["Ala1 -- Gly2"], -delta_no_flip.loc["Ala1 -- Gly2"])


def test_build_residue_graph_top_k():
    bond_deltas = pd.Series(
        {
            "Ala1 -- Gly2": 0.9,
            "Ala1 -- Ser3": -0.1,
            "Thr4 -- Val5": 0.2,
            "Thr4 -- Met6": -0.4,
        }
    )

    graph = build_residue_graph_from_bond_deltas(bond_deltas, top_k_edges=2, directed=False)
    # Only 2 edges should be present.
    assert len(graph.edges) == 2


def test_build_residue_graph_with_communities_detects_multiple_modules():
    # Two disconnected modules should yield at least two communities.
    bond_values = pd.Series(
        {
            "ALA1 -- GLY2": 0.9,
            "GLY2 -- SER3": 0.8,
            "THR10 -- VAL11": 0.85,
            "VAL11 -- LYS12": 0.75,
        }
    )
    result = build_residue_graph_with_communities(
        bond_values,
        top_k_edges=10,
        directed=False,
    )
    community_ids = set(result.community_by_node.values())
    assert len(result.graph.nodes) >= 6
    assert len(community_ids) >= 2
    assert np.isfinite(result.modularity)


def test_community_tables_have_required_columns():
    bond_values = pd.Series(
        {
            "ALA1 -- GLY2": 0.9,
            "GLY2 -- SER3": 0.8,
            "THR10 -- VAL11": 0.85,
        }
    )
    result = build_residue_graph_with_communities(bond_values, top_k_edges=10, directed=False)
    nodes_df, summary_df = community_tables(result)
    assert {"residue", "community_id", "node_size"}.issubset(set(nodes_df.columns))
    assert {"community_id", "num_nodes", "intra_edge_weight_sum", "modularity"}.issubset(
        set(summary_df.columns)
    )


def test_export_residue_graph_html_includes_summary_and_toggle(tmp_path):
    bond_values = pd.Series({"ALA1 -- GLY2": 0.5, "GLY2 -- SER3": -0.4})
    result = build_residue_graph_with_communities(bond_values, top_k_edges=5, directed=False)
    out = tmp_path / "graph.html"
    export_residue_graph_html(
        result.graph,
        out,
        summary_stats={
            "n_nodes": len(result.graph.nodes),
            "n_edges": len(result.graph.edges),
            "n_communities": len(set(result.community_by_node.values())),
            "modularity": result.modularity,
        },
    )
    text = out.read_text(encoding="utf-8")
    assert "Show inter-community edges only" in text
    assert "Modularity:" in text


def test_export_residue_graph_html_can_hide_delta_legend(tmp_path):
    bond_values = pd.Series({"ALA1 -- GLY2": 0.2})
    result = build_residue_graph_with_communities(bond_values, top_k_edges=5, directed=False)
    out = tmp_path / "graph_no_legend.html"
    export_residue_graph_html(result.graph, out, show_delta_legend=False)
    text = out.read_text(encoding="utf-8")
    assert "Edge Δ > 0" not in text

