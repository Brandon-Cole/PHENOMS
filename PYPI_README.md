# PHENOMS

Python-based Hydrogen-Deuterium Exchange of Molecular Dynamics Simulations

PHENOMS analyzes backbone H-bond networks from MD trajectories for HDX-MS–style
interpretation: Rust-accelerated Baker–Hubbard detection, replicate/`ComparisonSet`
workflows, occupancy heatmaps, differential protection, and connectivity exports.

## Install

```bash
pip install phenoms
```

## Quick start

```python
from phenoms import SimulationSet, default_output_root

sim = SimulationSet(
    pdb_files=["rep1.pdb", "rep2.pdb"],
    sub_frames=100,
    output_dir=default_output_root() / "my_run",
)
sim.run()
```

A CLI is also available: `phenoms prep`, `phenoms run`, `phenoms compare`.

**Docs:** https://brandon-cole.github.io/PHENOMS/
**Source:** https://github.com/Brandon-Cole/PHENOMS

## License

MIT
