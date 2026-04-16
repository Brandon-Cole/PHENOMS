"""
Unit tests for phenoms.analysis.
"""

import numpy as np
import pandas as pd
import pytest

from phenoms.analysis import (
    extract_residue_numbers,
    create_pivot_table,
    calculate_bond_statistics,
)


class TestExtractResidueNumbers:
    def test_simple(self):
        assert extract_residue_numbers("Ala1 -- Gly2") == (1, 2)

    def test_three_digit_residue(self):
        assert extract_residue_numbers("Ser100 -- Asp101") == (100, 101)

    def test_order_preserved(self):
        assert extract_residue_numbers("Gly5 -- Ala3") == (5, 3)


class TestCreatePivotTable:
    def test_basic(self):
        df = pd.DataFrame({
            "Bond Label": ["Ala1 -- Gly2", "Ala1 -- Gly2", "Gly2 -- Ser3"],
            "Frame": [0, 1, 0],
        })
        labels = ["Ala1 -- Gly2", "Gly2 -- Ser3"]
        pivot = create_pivot_table(df, labels)
        assert pivot.index.tolist() == labels
        assert pivot.shape[1] == 2  # frames 0 and 1
        assert pivot.loc["Ala1 -- Gly2", 0] == 1
        assert pivot.loc["Ala1 -- Gly2", 1] == 1
        assert pivot.loc["Gly2 -- Ser3", 0] == 1
        assert pivot.loc["Gly2 -- Ser3", 1] == 0

    def test_reindex_fill_zero(self):
        df = pd.DataFrame({
            "Bond Label": ["Ala1 -- Gly2"],
            "Frame": [0],
        })
        labels = ["Ala1 -- Gly2", "Gly2 -- Ser3"]
        pivot = create_pivot_table(df, labels)
        assert pivot.loc["Gly2 -- Ser3", 0] == 0


class TestCalculateBondStatistics:
    def test_empty_after_threshold(self):
        pivot = pd.DataFrame([[0, 0], [0, 0]], index=["A", "B"], columns=[0, 1])
        life, breaks = calculate_bond_statistics(pivot, threshold=0.5)
        assert life.empty and breaks.empty

    def test_one_bond_above_threshold(self):
        # one bond present in both frames -> lifetime 2, breaks 0
        pivot = pd.DataFrame(
            [[1, 1], [0, 0]],
            index=["Ala1 -- Gly2", "Gly2 -- Ser3"],
            columns=[0, 1],
        )
        life, breaks = calculate_bond_statistics(pivot, threshold=0.5)
        assert "Ala1 -- Gly2" in life.index
        assert life["Ala1 -- Gly2"] == 2.0
        assert breaks["Ala1 -- Gly2"] == 0

    def test_consecutive_runs(self):
        # present in 0,1 then 4,5 -> mean run length = 2
        pivot = pd.DataFrame(
            [[1, 1, 0, 1, 1]],
            index=["Ala1 -- Gly2"],
            columns=[0, 1, 2, 3, 4],
        )
        life, breaks = calculate_bond_statistics(pivot, threshold=0.1)
        assert life["Ala1 -- Gly2"] == 2.0
        assert breaks["Ala1 -- Gly2"] == 1
