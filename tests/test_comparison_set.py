import pandas as pd

from phenoms.comparison_set import ComparisonSet


class _DummySet:
    def __init__(self, pivot_tables, resid_range=None):
        self._pivot_tables = pivot_tables
        self.resid_range = resid_range
        self.pdb_files = ["dummy.pdb"]

    def get_pivot_tables(self):
        return self._pivot_tables

    def get_plot_pivot_tables(self):
        return self._pivot_tables

    def get_plot_bond_labels_sorted(self):
        labels = set()
        for pt in self._pivot_tables:
            labels |= set(pt.index.tolist())
        return sorted(labels)

    def _plot_region_str(self):
        return "Entire protein"


def test_export_connectivity_community_graph_html_writes_outputs(tmp_path):
    labels = ["ALA1 -- GLY2", "GLY2 -- SER3", "THR10 -- VAL11"]
    frames = [0, 1, 2]
    pt_a = pd.DataFrame(
        [[1, 1, 0], [1, 0, 0], [0, 0, 1]],
        index=labels,
        columns=frames,
    )
    pt_b = pd.DataFrame(
        [[0, 1, 1], [0, 1, 1], [1, 1, 0]],
        index=labels,
        columns=frames,
    )
    set_a = _DummySet([pt_a])
    set_b = _DummySet([pt_b])
    comp = ComparisonSet(set_a, set_b, label_a="a", label_b="b")
    comp.compare(flip_difference=True)

    html_out = tmp_path / "community_diff.html"
    nodes_csv = tmp_path / "community_nodes_diff.csv"
    summary_csv = tmp_path / "community_summary_diff.csv"
    comp.export_connectivity_community_graph_html(
        str(html_out),
        graph_mode="diff",
        community_nodes_csv_path=str(nodes_csv),
        community_summary_csv_path=str(summary_csv),
    )

    assert html_out.exists()
    assert nodes_csv.exists()
    assert summary_csv.exists()
    assert not pd.read_csv(nodes_csv).empty
    assert not pd.read_csv(summary_csv).empty
