"""
I/O and trajectory loading for PHENOMS.

Supports multi-frame PDBs (simple path) and native MD formats (trajectory + topology).
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path
from typing import Sequence

import mdtraj as md

_TPR_HELP = (
    "GROMACS .tpr topology detected ({tpr}), but MDTraj cannot parse .tpr directly.\n"
    "Install the optional MDAnalysis dependency to read it automatically:\n"
    "  pip install \"phenoms[gromacs]\"\n"
    "or convert it yourself first:\n"
    "  gmx trjconv -s {tpr} -o topol.gro -dump 0"
)


def tpr_to_mdtraj_topology(tpr_path) -> Path:
    """
    Convert a GROMACS .tpr into a throwaway PDB that MDTraj can use as a topology.

    Uses MDAnalysis's pure-Python .tpr parser (no GROMACS binary required) to read
    atoms/bonds/resids, then writes them out as PDB. Caller is responsible for
    deleting the returned path.
    """
    tpr_path = Path(tpr_path)
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise ValueError(_TPR_HELP.format(tpr=tpr_path)) from exc

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        universe = mda.Universe(str(tpr_path))
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
        tmp.close()
        universe.atoms.write(tmp.name)
    return Path(tmp.name)


def load_trajectory(path, top=None):
    """
    Load a trajectory from a PDB or from a native trajectory + topology pair.

    Parameters
    ----------
    path : str or path-like
        Trajectory file (``.pdb``, ``.xtc``, ``.trr``, ``.dcd``, ``.nc``, …).
    top : str or path-like or None
        Topology for non-PDB trajectories (``.pdb``, ``.prmtop``, ``.gro``, ``.tpr``, …).
        ``.tpr`` is converted on the fly via MDAnalysis (optional dependency).
        Ignored when ``path`` is already a multi-frame PDB with topology embedded.

    Returns
    -------
    mdtraj.Trajectory
    """
    path = str(path)
    if top is None:
        return md.load(path)

    top = Path(top)
    if top.suffix.lower() == ".tpr":
        converted = tpr_to_mdtraj_topology(top)
        try:
            return md.load(path, top=str(converted))
        finally:
            converted.unlink(missing_ok=True)

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
