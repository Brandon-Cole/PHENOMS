"""
I/O and trajectory loading for PHENOMS.
Supports whole-protein or residue-range selection.
"""

import mdtraj as md


def load_trajectory(pdb_file):
    """
    Load a trajectory from a PDB file (whole protein, no selection).

    Parameters
    ----------
    pdb_file : str or path-like
        Path to the PDB file.

    Returns
    -------
    mdtraj.Trajectory
        Full trajectory.
    """
    return md.load(pdb_file)


def load_and_select_residues(pdb_file, resid_range=None):
    """
    Load trajectory and optionally select a residue range.
    If resid_range is None, select all protein atoms (whole protein).

    Parameters
    ----------
    pdb_file : str or path-like
        Path to the PDB file.
    resid_range : tuple (int, int) or None
        (start, end) residue numbers (1-based). None = whole protein.

    Returns
    -------
    mdtraj.Trajectory
        Loaded (and optionally sliced) trajectory.
    """
    trajectory = md.load(pdb_file)
    if resid_range is None:
        selected_atoms = trajectory.top.select("protein")
    else:
        selected_atoms = trajectory.top.select(
            f"resid {resid_range[0]} to {resid_range[1]}"
        )
    return trajectory.atom_slice(selected_atoms)
