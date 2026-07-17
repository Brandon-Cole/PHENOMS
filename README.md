# PHENOMS

Python-based Hydrogen-Deuterium Exchange of Molecular Dynamics Simulations

![PHENOMS workflow](./assets/phenoms.png)

**Documentation:** [https://brandon-cole.github.io/PHENOMS/](https://brandon-cole.github.io/PHENOMS/)

PHENOMS analyzes backbone H-bond networks from MD trajectories for HDX-MS–style interpretation: Rust-accelerated Baker–Hubbard detection, replicate/`ComparisonSet` workflows, occupancy heatmaps, differential protection, and connectivity exports. Backbone N–O is the default; all-bond mode and native traj+topology inputs are optional.

## Install

```bash
conda env create -f environment.yml
conda activate phenoms
pip install -e .
```

Rust kernel (recommended):

```bash
pip install maturin && maturin develop --release
```

## Quick start

```python
from phenoms import SimulationSet, ComparisonSet, default_output_root

sim = SimulationSet(
    pdb_files=["rep1.pdb", "rep2.pdb"],
    resid_range=(50, 70),
    sub_frames=100,
    output_dir=default_output_root() / "my_run",
)
sim.run()

a = SimulationSet(["a1.pdb", "a2.pdb"], sub_frames=100)
b = SimulationSet(["b1.pdb", "b2.pdb"], sub_frames=100)
cmp = ComparisonSet(a, b, label_a="apo", label_b="holo")
a.run(); b.run()
cmp.compare()
cmp.export_comparison_artifacts(default_output_root() / "comparison")
```

CLI (optional): `phenoms prep`, `phenoms run`, `phenoms compare` — see the [docs](https://brandon-cole.github.io/PHENOMS/).

Outputs default to `./phenom_outputs/` (or `$PHENOMS_OUTPUT_DIR`).

## Benchmarks

Kernel timing (Docker, 4 CPUs, 3×250 frames, 2,722 atoms): Rust **0.047 s** vs MDTraj **0.803 s** (~17×) vs MDAnalysis **1.126 s** (~24×). Reproduce via [docker/README.md](./docker/README.md).

## License

[MIT](./LICENSE)
