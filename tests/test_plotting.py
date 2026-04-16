"""plot_difference threshold / autocorr helpers."""

import pandas as pd
import pytest

from phenoms.plotting import suggest_difference_threshold_autocorr


def test_suggest_difference_threshold_autocorr_bounds():
    y = [0.1 * i + 0.05 * ((-1) ** i) for i in range(20)]
    df = pd.DataFrame(
        {
            "Difference": y,
            "Difference_clipped": y,
            "Residue Number": range(20),
            "Donor Residue": [f"ALA{i}" for i in range(20)],
        }
    )
    t = suggest_difference_threshold_autocorr(df, kappa=1.96)
    assert 0.05 <= t <= 0.35


def test_suggest_difference_threshold_short_series():
    df = pd.DataFrame(
        {
            "Difference": [0.0, 0.5],
            "Difference_clipped": [0.0, 0.5],
            "Residue Number": [1, 2],
        }
    )
    assert suggest_difference_threshold_autocorr(df) == 0.2


@pytest.mark.parametrize("backend", ["agg"])
def test_plot_difference_return_meta(backend):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use(backend)
    import matplotlib.pyplot as plt

    from phenoms.plotting import plot_difference

    raw = [0.01, 0.5, -0.01, -0.4, 0.0]
    df = pd.DataFrame(
        {
            "Residue Number": [1, 2, 3, 4, 5],
            "Difference": raw,
            "Difference_clipped": raw,
            "Donor Residue": ["A1", "A2", "A3", "A4", "A5"],
        }
    )
    plt.close("all")
    meta = plot_difference(
        df,
        diff_threshold=0.2,
        diff_threshold_mode="manual",
        impute_small_differences=True,
        save_path=None,
        return_meta=True,
    )
    assert meta["diff_threshold_mode"] == "manual"
    assert meta["effective_threshold"] == 0.2
    assert meta["impute_small_differences"] is True
    plt.close("all")

    meta2 = plot_difference(
        df,
        diff_threshold_mode="autocorr",
        impute_small_differences=False,
        save_path=None,
        return_meta=True,
    )
    assert meta2["diff_threshold_mode"] == "autocorr"
    assert isinstance(meta2["effective_threshold"], float)
    plt.close("all")
