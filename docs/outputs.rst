Outputs & benchmarks
====================

Core analysis outputs
---------------------

Written under ``output_dir/raw_data/`` when ``output_dir`` is set (or via the CLI):

* **``*_hbonds.csv``** — frame-level donor/hydrogen/acceptor indices, frame, bond label
* **``*_occupancy.csv``** — present-frame counts and occupancy fractions
* **``*_pivot.csv``** — bond-label × frame matrices for heatmaps and comparisons
* **``manifest.json``** — run metadata (inputs, ``sub_frames``, ``resid_range``,
  ``backbone_only``, …)
* **``qc_report.json``** — present when QC was enabled

Differential / comparison outputs
---------------------------------

Via :meth:`phenoms.ComparisonSet.export_comparison_artifacts` and related methods:

* **``comparison.csv``** — per-donor occupancy in each condition, raw and clipped
  difference (clipped to ``[-1, 1]`` for visualization / B-factors)
* **Difference plots** — :meth:`~phenoms.ComparisonSet.plot_difference` with manual
  or autocorrelation-based thresholds
* **Connectivity HTML** — residue graphs from occupancy or Δ occupancy
  (:meth:`~phenoms.ComparisonSet.export_connectivity_graph_html`,
  community-aware variant available)
* **PDB B-factors** — map differential or single-set metrics for PyMOL/Chimera

Default location
----------------

``./phenom_outputs/`` or ``$PHENOMS_OUTPUT_DIR``. Prefer explicit ``output_dir`` /
``--output-dir`` for project-organized runs.

Kernel benchmarks
-----------------

Local:

.. code-block:: bash

   python scripts/benchmark_kernel.py --sub-frames 250

Docker reproduction: see ``docker/README.md`` in the repository.

Reference (4-CPU Docker, 3 replicates × 250 frames, 2,722 atoms each):

======= ============================ ============
Method  Total time (s)               Rust speedup
======= ============================ ============
Rust    0.047                        1.00×
MDTraj  0.803                        17.19×
MDA     1.126                        24.12×
======= ============================ ============
