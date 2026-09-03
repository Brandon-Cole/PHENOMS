"""
Shared pytest fixtures.
"""

import pytest


@pytest.fixture(autouse=True)
def _isolate_phenoms_output_dir(tmp_path, monkeypatch):
    """
    SimulationSet now writes a standard output bundle by default. Point
    PHENOMS_OUTPUT_DIR at a per-test tmp dir so that default write never
    touches the repo working directory or leaks between tests.
    """
    monkeypatch.setenv("PHENOMS_OUTPUT_DIR", str(tmp_path / "phenom_outputs"))
