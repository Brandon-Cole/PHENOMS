"""
Hydrogen bond detection (Baker–Hubbard).

- **Rust** (``phenoms_hbond_rs``): fast batched geometry when the extension is built;
  falls back to MDTraj otherwise.
- **Polars** (required): per-bond occupancy summaries in :func:`hbond_occupancy_table`
  and any export path that uses it (Rust does not call Polars; detection stays in Rust).
- Backbone mode (default): donor/acceptor filtered to backbone-relevant N–O pairs.
- All-bonds mode: include all donor/acceptor atom classes supported by topology.
"""

import numpy as np
import pandas as pd
import mdtraj as md
from joblib import Parallel, delayed
from tqdm import tqdm

try:
    from phenoms_hbond_rs import run_baker_hubbard as _run_baker_hubbard_rs
except ImportError:
    _run_baker_hubbard_rs = None
try:
    from phenoms_hbond_rs import run_baker_hubbard_with_threads as _run_baker_hubbard_rs_threads
except ImportError:
    _run_baker_hubbard_rs_threads = None
try:
    from phenoms_hbond_rs import run_baker_hubbard_brute as _run_baker_hubbard_rs_brute
except ImportError:
    _run_baker_hubbard_rs_brute = None

import polars as pl


def label_hbond(hbond, trajectory):
    """
    Label a hydrogen bond as 'residue1 -- residue2' only if it is N–O (backbone-relevant).
    Returns None otherwise.

    Parameters
    ----------
    hbond : tuple (donor_idx, hydrogen_idx, acceptor_idx)
    trajectory : mdtraj.Trajectory (single frame or full; only topology is used)

    Returns
    -------
    str or None
        e.g. 'ALA1 -- GLY2' or None
    """
    atom1 = trajectory.topology.atom(hbond[0])
    atom2 = trajectory.topology.atom(hbond[2])
    atom1_name = atom1.name
    atom2_name = atom2.name

    if atom1_name == "N" and atom2_name == "O":
        donor_residue = atom1.residue
        acceptor_residue = atom2.residue
        return f"{donor_residue} -- {acceptor_residue}"
    if atom1_name == "O" and atom2_name == "N":
        donor_residue = atom2.residue
        acceptor_residue = atom1.residue
        return f"{donor_residue} -- {acceptor_residue}"
    return None


def _topology_donor_acceptor_lists(trajectory, backbone_only=True):
    """
    Build (donor_h_pairs, acceptor_indices) from MDTraj topology for Rust.
    Returns (None, None) if topology has no bonds (e.g. minimal PDB).
    """
    top = trajectory.topology
    n_atoms = top.n_atoms
    if backbone_only:
        acceptor_indices = [i for i in range(n_atoms) if top.atom(i).name in ("N", "O")]
    else:
        acceptor_indices = []
        for i in range(n_atoms):
            atom = top.atom(i)
            sym = atom.element.symbol if atom.element is not None else ""
            if sym in ("N", "O", "S"):
                acceptor_indices.append(i)
    bonds = list(top.bonds)
    if not bonds:
        return None, None
    donor_h_pairs = []
    for a1, a2 in bonds:
        i, j = a1.index, a2.index
        name_i, name_j = top.atom(i).name, top.atom(j).name
        try:
            elem_i = top.atom(i).element.symbol
            elem_j = top.atom(j).element.symbol
        except Exception:
            elem_i, elem_j = "", ""
        is_h_i = name_i == "H" or elem_i == "H"
        is_h_j = name_j == "H" or elem_j == "H"
        if backbone_only:
            donor_ok_i = name_i in ("N", "O")
            donor_ok_j = name_j in ("N", "O")
        else:
            donor_ok_i = elem_i in ("N", "O", "S")
            donor_ok_j = elem_j in ("N", "O", "S")

        if is_h_i and donor_ok_j:
            donor_h_pairs.append((j, i))
        elif is_h_j and donor_ok_i:
            donor_h_pairs.append((i, j))
    if not donor_h_pairs:
        return None, None
    return donor_h_pairs, acceptor_indices


def _process_frames_rust(trajectory, sub_frames=None, backbone_only=True, n_jobs=1):
    """Run Baker-Hubbard in Rust for all frames."""
    total_frames = trajectory.n_frames
    n_process = min(sub_frames, total_frames) if sub_frames else total_frames
    xyz = trajectory.xyz[:n_process].astype(np.float64)
    n_frames = xyz.shape[0]
    n_atoms = xyz.shape[1]
    xyz_flat = np.ascontiguousarray(xyz.reshape(-1))

    donor_h_pairs, acceptor_indices = _topology_donor_acceptor_lists(
        trajectory,
        backbone_only=backbone_only,
    )
    if donor_h_pairs is None or acceptor_indices is None:
        return None

    if _run_baker_hubbard_rs_threads is not None and int(n_jobs) > 1:
        raw = _run_baker_hubbard_rs_threads(
            xyz_flat,
            n_frames,
            n_atoms,
            donor_h_pairs,
            acceptor_indices,
            int(n_jobs),
        )
    else:
        raw = _run_baker_hubbard_rs(
            xyz_flat,
            n_frames,
            n_atoms,
            donor_h_pairs,
            acceptor_indices,
        )
    if raw.size == 0:
        df = pd.DataFrame(
            columns=["Frame", "Donor Index", "Hydrogen Index", "Acceptor Index", "Bond Label"],
        )
        return df
    df = pd.DataFrame(
        raw,
        columns=["Frame", "Donor Index", "Hydrogen Index", "Acceptor Index"],
    )
    donors = df["Donor Index"].to_numpy(dtype=np.int64, copy=False)
    acceptors = df["Acceptor Index"].to_numpy(dtype=np.int64, copy=False)
    mult = np.int64(n_atoms + 1)
    pair_ids = donors * mult + acceptors
    unique_ids, inverse = np.unique(pair_ids, return_inverse=True)
    unique_pairs = [(int(pid // mult), int(pid % mult)) for pid in unique_ids]
    pair_to_label = _pair_label_map(trajectory, unique_pairs, backbone_only=backbone_only)
    labels_for_unique = np.empty(len(unique_ids), dtype=object)
    for i, (d, a) in enumerate(unique_pairs):
        labels_for_unique[i] = pair_to_label.get((d, a))
    df["Bond Label"] = labels_for_unique[inverse]
    if backbone_only:
        df = df.dropna(subset=["Bond Label"])
    return df


def _pair_label_map(trajectory, donor_acceptor_pairs, *, backbone_only=True):
    """
    Build a cached label map for (donor_atom_idx, acceptor_atom_idx) pairs.
    """
    top = trajectory.topology
    out = {}
    for donor_idx, acceptor_idx in donor_acceptor_pairs:
        donor_atom = top.atom(int(donor_idx))
        acceptor_atom = top.atom(int(acceptor_idx))
        if backbone_only:
            donor_name = donor_atom.name
            acceptor_name = acceptor_atom.name
            if donor_name == "N" and acceptor_name == "O":
                out[(int(donor_idx), int(acceptor_idx))] = (
                    f"{donor_atom.residue} -- {acceptor_atom.residue}"
                )
            elif donor_name == "O" and acceptor_name == "N":
                out[(int(donor_idx), int(acceptor_idx))] = (
                    f"{acceptor_atom.residue} -- {donor_atom.residue}"
                )
            else:
                out[(int(donor_idx), int(acceptor_idx))] = None
        else:
            out[(int(donor_idx), int(acceptor_idx))] = (
                f"{donor_atom.residue}:{donor_atom.name} -- {acceptor_atom.residue}:{acceptor_atom.name}"
            )
    return out


def process_frame(frame, trajectory):
    """
    Run Baker–Hubbard on one frame and label N–O bonds (MDTraj path).
    """
    hbonds = md.baker_hubbard(trajectory[frame], periodic=False, freq=0)
    frame_df = pd.DataFrame(
        hbonds,
        columns=["Donor Index", "Hydrogen Index", "Acceptor Index"],
    )
    bond_labels = [
        label_hbond(
            (int(bond["Donor Index"]), int(bond["Hydrogen Index"]), int(bond["Acceptor Index"])),
            trajectory[frame],
        )
        for _, bond in frame_df.iterrows()
    ]
    frame_df["Bond Label"] = bond_labels
    frame_df["Frame"] = frame
    return frame_df


def label_hbond_all(hbond, trajectory):
    """
    Label any detected H-bond as 'DONORRES:ATOM -- ACCEPTORRES:ATOM'.
    """
    donor_atom = trajectory.topology.atom(hbond[0])
    acceptor_atom = trajectory.topology.atom(hbond[2])
    donor_residue = donor_atom.residue
    acceptor_residue = acceptor_atom.residue
    return f"{donor_residue}:{donor_atom.name} -- {acceptor_residue}:{acceptor_atom.name}"


def process_frame_all(frame, trajectory):
    """
    Run Baker-Hubbard on one frame and keep all detected donor/acceptor types.
    """
    hbonds = md.baker_hubbard(trajectory[frame], periodic=False, freq=0)
    frame_df = pd.DataFrame(
        hbonds,
        columns=["Donor Index", "Hydrogen Index", "Acceptor Index"],
    )
    bond_labels = [
        label_hbond_all(
            (int(bond["Donor Index"]), int(bond["Hydrogen Index"]), int(bond["Acceptor Index"])),
            trajectory[frame],
        )
        for _, bond in frame_df.iterrows()
    ]
    frame_df["Bond Label"] = bond_labels
    frame_df["Frame"] = frame
    return frame_df


def process_frames(trajectory, sub_frames=None, n_jobs=1, use_rust=True):
    """
    Process all (or sub_frames) frames; return N–O H-bonds only.
    Uses Rust when available and topology has bonds; otherwise MDTraj (parallel).

    The Rust path uses the same geometric criteria as MDTraj Baker–Hubbard:
    r(H···A) < 0.25 nm and ∠D–H–A > 120° (angle at H). Use ``use_rust=False`` only
    for debugging or if the extension is not built.
    """
    total_frames = trajectory.n_frames
    if sub_frames is None:
        sub_frames = total_frames
    n_process = min(sub_frames, total_frames)

    if use_rust and _run_baker_hubbard_rs is not None:
        df = _process_frames_rust(
            trajectory,
            sub_frames=n_process,
            backbone_only=True,
            n_jobs=n_jobs,
        )
        if df is not None:
            return df
    hbond_frames = Parallel(n_jobs=n_jobs)(
        delayed(process_frame)(frame, trajectory)
        for frame in tqdm(
            range(n_process),
            desc="Processing Frames",
            unit="frame",
        )
    )
    hbonds_df = pd.concat(hbond_frames, ignore_index=True)
    hbonds_df = hbonds_df.dropna(subset=["Bond Label"])
    return hbonds_df


def process_frames_all(trajectory, sub_frames=None, n_jobs=1, use_rust=True):
    """
    Process all (or sub_frames) frames; return all detected H-bonds (not backbone-only).
    Uses Rust when available and topology has bonds; otherwise MDTraj (parallel).
    """
    total_frames = trajectory.n_frames
    if sub_frames is None:
        sub_frames = total_frames
    n_process = min(sub_frames, total_frames)

    if use_rust and _run_baker_hubbard_rs is not None:
        df = _process_frames_rust(
            trajectory,
            sub_frames=n_process,
            backbone_only=False,
            n_jobs=n_jobs,
        )
        if df is not None:
            return df

    hbond_frames = Parallel(n_jobs=n_jobs)(
        delayed(process_frame_all)(frame, trajectory)
        for frame in tqdm(
            range(n_process),
            desc="Processing Frames",
            unit="frame",
        )
    )
    hbonds_df = pd.concat(hbond_frames, ignore_index=True)
    hbonds_df = hbonds_df.dropna(subset=["Bond Label"])
    return hbonds_df


def hbond_occupancy_table(hbonds_df):
    """
    Build a per-bond occupancy summary from a frame-level H-bond DataFrame.

    Returns columns:
      - Bond Label
      - Present Frames
      - Total Frames
      - Occupancy

    Uses Polars for the aggregation (required dependency).
    """
    if hbonds_df is None or hbonds_df.empty:
        return pd.DataFrame(columns=["Bond Label", "Present Frames", "Total Frames", "Occupancy"])

    total_frames = int(hbonds_df["Frame"].nunique())
    sub = hbonds_df[["Frame", "Bond Label"]]
    lf = pl.from_pandas(sub, nan_to_null=True)
    out_pl = (
        lf.group_by("Bond Label")
        .agg(pl.col("Frame").n_unique().alias("Present Frames"))
        .with_columns(
            pl.lit(total_frames).alias("Total Frames"),
            (pl.col("Present Frames") / float(max(total_frames, 1))).alias("Occupancy"),
        )
        .sort("Occupancy", descending=True)
    )
    return pd.DataFrame(out_pl.to_dicts())


def export_hbond_occupancy_csv(hbonds_df, output_csv_path):
    """
    Write per-bond occupancy summary CSV.
    """
    table = hbond_occupancy_table(hbonds_df)
    table.to_csv(output_csv_path, index=False)
    return table
