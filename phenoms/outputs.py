"""
User-visible output locations.

Set ``PHENOMS_OUTPUT_DIR`` to a writable directory to keep all default artifacts
outside the installed package tree. If unset, defaults to ``./phenom_outputs``
(relative to the process current working directory).
"""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path


def default_output_root() -> Path:
    """
    Root directory for PHENOMS artifacts (benchmarks, test-plot runs, exports).

    Resolution order:
    1. Environment variable ``PHENOMS_OUTPUT_DIR`` (expanded user path).
    2. Otherwise ``Path.cwd() / "phenom_outputs"``.
    """
    env = os.environ.get("PHENOMS_OUTPUT_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (Path.cwd() / "phenom_outputs").resolve()


def timestamped_run_dir(prefix: str = "run") -> Path:
    """
    A fresh directory under :func:`default_output_root`, named
    ``<prefix>_<YYYYmmdd_HHMMSS_ffffff>``.

    Used as the default ``output_dir`` for :class:`~phenoms.SimulationSet` so a
    run leaves the standard artifact bundle on disk without the caller having
    to name a path, while never colliding with a previous run's output.
    """
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    return default_output_root() / f"{prefix}_{stamp}"
