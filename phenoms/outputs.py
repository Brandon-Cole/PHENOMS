"""
User-visible output locations.

Set ``PHENOMS_OUTPUT_DIR`` to a writable directory to keep all default artifacts
outside the installed package tree. If unset, defaults to ``./phenom_outputs``
(relative to the process current working directory).
"""

from __future__ import annotations

import os
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
