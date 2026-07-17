#!/usr/bin/env python3
import gc
import os
import sys
import tempfile
import time

import mdtraj as md
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

if "/app" not in sys.path:
    sys.path.insert(0, "/app")
from phenoms.hbond import _run_baker_hubbard_rs_threads, _topology_donor_acceptor_lists


def main():
    traj = os.environ.get("DBG_TRAJ", "/prep/plcg2_preprocessed_1000f.xtc")
    top = os.environ.get("DBG_TOP", "/prep/plcg2_preprocessed_top.pdb")
    frames = int(os.environ.get("DBG_FRAMES", "1000"))
    workers = int(os.environ.get("DBG_WORKERS", "4"))

    mtop = md.load_topology(top)
    residues = [r for r in mtop.residues]
    windows = []
    for i in range(10):
        keep = len(residues) - i * 100
        if keep <= 0:
            break
        windows.append((f"W{i}", [r.index for r in residues[:keep]]))
    print(f"windows={len(windows)} frames={frames} workers={workers}", flush=True)

    for wname, ridx in windows:
        t0 = time.perf_counter()
        ridx_set = set(ridx)
        atom_idx = [a.index for a in mtop.atoms if a.residue.index in ridx_set]
        sub = md.load(traj, top=top, atom_indices=atom_idx)[:frames]
        print(f"{wname} load atoms={sub.n_atoms} time={time.perf_counter()-t0:.2f}s", flush=True)

        # Rust stage
        donor_h_pairs, acceptors = _topology_donor_acceptor_lists(sub, backbone_only=False)
        xyz_nm = sub.xyz.astype(np.float64)
        xyz_flat = np.ascontiguousarray(xyz_nm.reshape(-1))
        t1 = time.perf_counter()
        raw = _run_baker_hubbard_rs_threads(
            xyz_flat,
            xyz_nm.shape[0],
            xyz_nm.shape[1],
            donor_h_pairs,
            acceptors,
            workers,
        )
        print(f"{wname} rust rows={0 if raw is None else int(raw.shape[0])} time={time.perf_counter()-t1:.2f}s", flush=True)

        # MDTraj stage
        t2 = time.perf_counter()
        hb_md = md.baker_hubbard(sub, periodic=False)
        print(f"{wname} mdtraj rows={len(hb_md)} time={time.perf_counter()-t2:.2f}s", flush=True)

        # MDA stage
        td = tempfile.mkdtemp(prefix="dbg_w_")
        pdb_path = os.path.join(td, f"{wname}.pdb")
        sub.save_pdb(pdb_path)
        u = mda.Universe(pdb_path)
        hb = HBA(
            universe=u,
            donors_sel="protein and element N",
            hydrogens_sel="protein and element H",
            acceptors_sel="protein and (element O or element N)",
            d_h_cutoff=1.2,
        )
        t3 = time.perf_counter()
        hb.run(start=0, stop=frames, step=1, verbose=False, n_workers=workers, backend="multiprocessing")
        rows = 0 if hb.results.hbonds is None else len(hb.results.hbonds)
        print(f"{wname} mda rows={rows} time={time.perf_counter()-t3:.2f}s", flush=True)

        del hb, u, hb_md, raw, xyz_nm, xyz_flat, sub
        gc.collect()

    print("done", flush=True)


if __name__ == "__main__":
    main()
