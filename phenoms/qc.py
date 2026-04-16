"""
Quality control utilities for replicate trajectory analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import mdtraj as md
import numpy as np


DEFAULT_MDP_KEYS = (
    "integrator",
    "dt",
    "nsteps",
    "tcoupl",
    "pcoupl",
    "constraints",
    "constraint-algorithm",
    "coulombtype",
    "rcoulomb",
    "vdwtype",
    "rvdw",
    "tc-grps",
    "tau-t",
    "ref-t",
    "tau-p",
    "ref-p",
)


def parse_mdp(path: str | Path) -> dict[str, str]:
    """Parse a GROMACS .mdp file into a normalized key/value dictionary."""
    out: dict[str, str] = {}
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, val = line.split("=", 1)
        out[key.strip().lower()] = val.strip()
    return out


def check_mdp_key_consistency(
    mdp_files: Iterable[str | Path],
    keys: Iterable[str] = DEFAULT_MDP_KEYS,
) -> dict:
    """
    Check that selected MDP keys are consistent across all provided files.
    Returns a report dictionary with per-key mismatches.
    """
    files = [str(Path(p)) for p in mdp_files]
    parsed = {f: parse_mdp(f) for f in files}
    norm_keys = [k.lower() for k in keys]
    mismatches: dict[str, dict[str, str | None]] = {}
    for k in norm_keys:
        vals = {f: parsed[f].get(k) for f in files}
        uniq = set(vals.values())
        if len(uniq) > 1:
            mismatches[k] = vals
    return {
        "ok": len(mismatches) == 0,
        "files": files,
        "checked_keys": norm_keys,
        "mismatches": mismatches,
    }


def assess_series_convergence(
    series: np.ndarray,
    *,
    last_fraction: float = 0.2,
    mean_tolerance: float = 0.05,
    slope_tolerance: float = 1e-4,
) -> dict:
    """
    Simple convergence check:
    - compare mean of penultimate vs final window
    - fit slope on final window
    """
    x = np.asarray(series, dtype=float)
    n = len(x)
    if n < 10:
        return {
            "converged": False,
            "reason": "too_few_points",
            "n_points": int(n),
        }
    w = max(5, int(round(n * last_fraction)))
    prev = x[-2 * w : -w] if n >= 2 * w else x[:w]
    last = x[-w:]
    prev_mean = float(np.mean(prev))
    last_mean = float(np.mean(last))
    denom = max(abs(prev_mean), 1e-8)
    rel_diff = abs(last_mean - prev_mean) / denom
    t = np.arange(w, dtype=float)
    slope = float(np.polyfit(t, last, 1)[0]) if w >= 2 else 0.0
    converged = rel_diff <= mean_tolerance and abs(slope) <= slope_tolerance
    reason = "ok" if converged else "window_shift_or_drift"
    return {
        "converged": bool(converged),
        "reason": reason,
        "n_points": int(n),
        "window_size": int(w),
        "prev_mean": prev_mean,
        "last_mean": last_mean,
        "relative_mean_diff": float(rel_diff),
        "last_window_slope": slope,
        "mean_tolerance": float(mean_tolerance),
        "slope_tolerance": float(slope_tolerance),
    }


def rmsd_convergence_report(
    traj: md.Trajectory,
    *,
    atom_selection: str = "protein and backbone",
    last_fraction: float = 0.2,
    mean_tolerance: float = 0.05,
    slope_tolerance: float = 1e-4,
) -> dict:
    """Compute RMSD-to-first-frame series and run a simple convergence check."""
    sel = traj.top.select(atom_selection)
    if len(sel) == 0:
        sel = traj.top.select("protein")
    if len(sel) == 0:
        sel = np.arange(traj.n_atoms, dtype=int)
    rmsd = md.rmsd(traj, traj, 0, atom_indices=sel)
    conv = assess_series_convergence(
        rmsd,
        last_fraction=last_fraction,
        mean_tolerance=mean_tolerance,
        slope_tolerance=slope_tolerance,
    )
    conv.update(
        {
            "metric": "rmsd_to_first_frame_nm",
            "atom_selection": atom_selection,
            "n_selected_atoms": int(len(sel)),
        }
    )
    return conv
