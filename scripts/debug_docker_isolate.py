#!/usr/bin/env python3
import gc
import os
import tempfile
import time

import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA


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
        t1 = time.perf_counter()
        hb.run(start=0, stop=frames, step=1, verbose=False, n_workers=workers, backend="multiprocessing")
        rows = 0 if hb.results.hbonds is None else len(hb.results.hbonds)
        print(f"{wname} mda rows={rows} time={time.perf_counter()-t1:.2f}s", flush=True)
        del hb, u, sub
        gc.collect()

    print("done", flush=True)


if __name__ == "__main__":
    main()
