"""compare_two_sets and clipped difference column."""

import pandas as pd

from phenoms.comparison import compare_two_sets, differentially_protected, prefer_difference_clipped_column


def test_compare_two_sets_difference_clipped():
    p1 = pd.DataFrame({"b1": [1.0, 0.0]}, index=["ALA1 -- GLY2", "GLY2 -- ALA3"])
    p2 = pd.DataFrame({"b1": [0.0, 1.0]}, index=["ALA1 -- GLY2", "GLY2 -- ALA3"])
    out = compare_two_sets([p1], [p2], label_a="a", label_b="b", donor_aggregation="mean")
    assert "Difference_clipped" in out.columns
    assert out["Difference_clipped"].max() <= 1.0
    assert out["Difference_clipped"].min() >= -1.0


def test_prefer_difference_clipped_column():
    df = pd.DataFrame({"Difference": [1.0], "Difference_clipped": [0.5]})
    assert prefer_difference_clipped_column(df) == "Difference_clipped"
    df2 = pd.DataFrame({"Difference": [1.0]})
    assert prefer_difference_clipped_column(df2) == "Difference"


def test_differentially_protected_uses_clipped_by_default():
    df = pd.DataFrame(
        {
            "Difference": [2.0],
            "Difference_clipped": [1.0],
            "Residue Number": [1],
        }
    )
    hit = differentially_protected(df, threshold=0.5)
    assert len(hit) == 1
    miss = differentially_protected(df, threshold=1.5)
    assert len(miss) == 0
