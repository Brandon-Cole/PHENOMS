"""
Engine-aware input discovery and trajectory normalization for PHENOMS.

This module converts replicate folders from GROMACS/OpenMM/AMBER layouts into
normalized multi-frame protein-only PDBs consumable by SimulationSet.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mdtraj as md


_TRAJ_EXTS = {".xtc", ".trr", ".dcd", ".nc", ".mdcrd", ".pdb"}


@dataclass(frozen=True)
class ReplicateInput:
    replicate_dir: Path
    engine: str
    trajectory_path: Path
    topology_path: Path | None


def _files_with_exts(directory: Path, exts: set[str]) -> list[Path]:
    return sorted(
        [p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in exts],
        key=lambda p: p.name.lower(),
    )


def _pick_first(paths: Iterable[Path], suffixes: tuple[str, ...]) -> Path | None:
    suffixes = tuple(s.lower() for s in suffixes)
    for p in paths:
        if p.suffix.lower() in suffixes:
            return p
    return None


def detect_replicate_input(replicate_dir: str | Path) -> ReplicateInput:
    """
    Detect simulation engine and enforce required files for one replicate folder.

    Required files:
    - gromacs: trajectory (.xtc/.trr) + topology (.tpr preferred, else .gro/.pdb)
    - openmm: trajectory (.dcd/.xtc) + topology (.pdb/.prmtop)
    - amber: trajectory (.nc/.mdcrd) + topology (.prmtop/.parm7)
    - pdb: one multi-model .pdb trajectory (already normalized input)
    """
    rep = Path(replicate_dir).expanduser().resolve()
    if not rep.is_dir():
        raise ValueError(f"Replicate path is not a directory: {rep}")

    files = [p for p in rep.iterdir() if p.is_file()]
    if not files:
        raise ValueError(f"No files found in replicate directory: {rep}")

    all_traj = _files_with_exts(rep, _TRAJ_EXTS)

    # GROMACS
    gmx_traj = _pick_first(all_traj, (".xtc", ".trr"))
    if gmx_traj is not None:
        gmx_top = _pick_first(files, (".tpr", ".gro", ".pdb"))
        if gmx_top is None:
            raise ValueError(
                f"GROMACS replicate {rep} missing topology file (.tpr/.gro/.pdb) for {gmx_traj.name}."
            )
        return ReplicateInput(rep, "gromacs", gmx_traj, gmx_top)

    # OpenMM
    openmm_traj = _pick_first(all_traj, (".dcd", ".xtc"))
    if openmm_traj is not None:
        openmm_top = _pick_first(files, (".pdb", ".prmtop"))
        if openmm_top is None:
            raise ValueError(
                f"OpenMM replicate {rep} missing topology file (.pdb/.prmtop) for {openmm_traj.name}."
            )
        return ReplicateInput(rep, "openmm", openmm_traj, openmm_top)

    # AMBER
    amber_traj = _pick_first(all_traj, (".nc", ".mdcrd"))
    if amber_traj is not None:
        amber_top = _pick_first(files, (".prmtop", ".parm7"))
        if amber_top is None:
            raise ValueError(
                f"AMBER replicate {rep} missing topology file (.prmtop/.parm7) for {amber_traj.name}."
            )
        return ReplicateInput(rep, "amber", amber_traj, amber_top)

    # PDB-only pathway (already normalized)
    pdbs = [p for p in files if p.suffix.lower() == ".pdb"]
    if len(pdbs) == 1:
        return ReplicateInput(rep, "pdb", pdbs[0], None)

    raise ValueError(
        "Could not detect replicate engine/input files in "
        f"{rep}. Expected GROMACS/OpenMM/AMBER files or a single .pdb."
    )


def discover_replicate_dirs(input_dir: str | Path) -> list[Path]:
    """
    Directory policy:
    - if input_dir itself looks like one replicate, use it
    - otherwise treat each direct child directory as a replicate
    """
    root = Path(input_dir).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"Input path is not a directory: {root}")

    try:
        detect_replicate_input(root)
        return [root]
    except ValueError:
        pass

    reps = sorted([p for p in root.iterdir() if p.is_dir()], key=lambda p: p.name.lower())
    if not reps:
        raise ValueError(
            f"No replicate directories found in {root}. "
            "Pass either one replicate dir, or a set dir containing replicate subdirs."
        )
    return reps


def _load_rep(rep: ReplicateInput):
    if rep.topology_path is None:
        return md.load(str(rep.trajectory_path))
    return md.load(str(rep.trajectory_path), top=str(rep.topology_path))


def _select_frames_by_time_ps(traj, start_ps: float | None, end_ps: float | None, frame_dt_ps: float | None):
    mask = None
    if start_ps is not None:
        m = traj.time >= float(start_ps)
        mask = m if mask is None else (mask & m)
    if end_ps is not None:
        m = traj.time <= float(end_ps)
        mask = m if mask is None else (mask & m)

    if mask is None:
        idx = list(range(traj.n_frames))
    else:
        idx = [i for i, keep in enumerate(mask) if bool(keep)]

    if not idx:
        raise ValueError("No frames selected after start/end filtering.")

    if frame_dt_ps is not None:
        t0 = traj.time[idx[0]]
        step = float(frame_dt_ps)
        picked = [idx[0]]
        last_t = t0
        for i in idx[1:]:
            ti = traj.time[i]
            if ti >= last_t + step - 1e-9:
                picked.append(i)
                last_t = ti
        idx = picked
    return idx


def normalize_replicate_to_pdb(
    replicate: ReplicateInput,
    out_pdb_path: str | Path,
    *,
    frame_dt_ps: float | None = 1000.0,
    start_ps: float | None = None,
    end_ps: float | None = None,
    apply_imaging: bool = True,
    center: bool = True,
    fit: bool = True,
) -> Path:
    """
    Normalize one replicate:
    image/unwrap (best-effort) -> center -> fit -> protein-only -> save PDB.
    """
    traj = _load_rep(replicate)
    protein_idx = traj.topology.select("protein")
    if len(protein_idx) > 0:
        traj = traj.atom_slice(protein_idx)

    if apply_imaging:
        try:
            traj.image_molecules(inplace=True)
        except Exception:
            # Best effort only; fallback keeps trajectory usable.
            pass

    if center:
        traj.center_coordinates()

    if fit and traj.n_frames > 1:
        atom_idx = list(range(traj.n_atoms))
        traj.superpose(traj, frame=0, atom_indices=atom_idx)

    frame_idx = _select_frames_by_time_ps(traj, start_ps=start_ps, end_ps=end_ps, frame_dt_ps=frame_dt_ps)
    traj = traj.slice(frame_idx)

    out = Path(out_pdb_path).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    traj.save_pdb(str(out))
    return out


def prepare_set_from_dir(
    input_dir: str | Path,
    prepared_dir: str | Path,
    *,
    frame_dt_ps: float | None = 1000.0,
    start_ps: float | None = None,
    end_ps: float | None = None,
    apply_imaging: bool = True,
    center: bool = True,
    fit: bool = True,
) -> list[str]:
    """
    Prepare one set directory into normalized PDB replicates.
    Returns sorted list of prepared PDB paths.
    """
    reps = discover_replicate_dirs(input_dir)
    out_root = Path(prepared_dir).expanduser().resolve()
    out_root.mkdir(parents=True, exist_ok=True)

    out_paths: list[str] = []
    for rep_dir in reps:
        rep = detect_replicate_input(rep_dir)
        out_path = out_root / f"{rep_dir.name}_protein_centered_fit_nopbc.pdb"
        normalize_replicate_to_pdb(
            rep,
            out_path,
            frame_dt_ps=frame_dt_ps,
            start_ps=start_ps,
            end_ps=end_ps,
            apply_imaging=apply_imaging,
            center=center,
            fit=fit,
        )
        out_paths.append(str(out_path))
    return sorted(out_paths)
