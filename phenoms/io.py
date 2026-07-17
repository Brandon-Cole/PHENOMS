"""
I/O and trajectory loading for PHENOMS.

Supports multi-frame PDBs (simple path) and native MD formats (trajectory + topology).
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import mdtraj as md


def load_trajectory(path, top=None):
    """
    Load a trajectory from a PDB or from a native trajectory + topology pair.

    Parameters
    ----------
    path : str or path-like
        Trajectory file (``.pdb``, ``.xtc``, ``.trr``, ``.dcd``, ``.nc``, …).
    top : str or path-like or None
        Topology for non-PDB trajectories (``.pdb``, ``.prmtop``, ``.gro``, ``.tpr``, …).
        Ignored when ``path`` is already a multi-frame PDB with topology embedded.

    Returns
    -------
    mdtraj.Trajectory
    """
    path = str(path)
    if top is None:
        return md.load(path)
    return md.load(path, top=str(top))


def load_and_select_residues(path, resid_range=None, top=None):
    """
    Load trajectory and optionally select a residue range.

    If ``resid_range`` is None, select all protein atoms (whole protein).

    Parameters
    ----------
    path : str or path-like
        Trajectory / PDB path.
    resid_range : tuple (int, int) or None
        (start, end) residue numbers (1-based). None = whole protein.
    top : str or path-like or None
        Topology for native trajectory formats.

    Returns
    -------
    mdtraj.Trajectory
        Loaded (and optionally sliced) trajectory.
    """
    trajectory = load_trajectory(path, top=top)
    if resid_range is None:
        selected_atoms = trajectory.top.select("protein")
    else:
        selected_atoms = trajectory.top.select(
            f"resid {resid_range[0]} to {resid_range[1]}"
        )
    return trajectory.atom_slice(selected_atoms)


def normalize_topology_list(
    n_replicates: int,
    *,
    topology=None,
    topologies: Sequence | None = None,
) -> list | None:
    """
    Resolve a single shared topology or per-replicate topologies.

    Returns ``None`` when no topology is needed (PDB-only inputs).
    """
    if topologies is not None and topology is not None:
        raise ValueError("Pass only one of topology= or topologies=, not both.")
    if topologies is not None:
        tops = [str(Path(t).expanduser()) if t is not None else None for t in topologies]
        if len(tops) != n_replicates:
            raise ValueError(
                f"topologies length ({len(tops)}) must match number of "
                f"replicates ({n_replicates})."
            )
        return tops
    if topology is not None:
        top = str(Path(topology).expanduser())
        return [top] * n_replicates
    return None
