Outputs
=======

Core analysis outputs
---------------------

``SimulationSet.run()`` writes this bundle by default — see :ref:`default-outputs`
below for exactly when and where. Under ``output_dir/raw_data/``:

* **``*_hbonds.csv``** — frame-level donor/hydrogen/acceptor indices, frame, bond label
* **``*_occupancy.csv``** — present-frame counts and occupancy fractions
* **``*_pivot.csv``** — bond-label × frame matrices for heatmaps and comparisons
* **``manifest.json``** — run metadata (inputs, ``sub_frames``, ``resid_range``,
  ``backbone_only``, …)
* **``qc_report.json``** — present when QC was enabled

And directly under ``output_dir``:

* **``plots/``** — per-replicate heatmaps, an aggregated heatmap, and (when
  ``bond_statistics_threshold`` is set) lifetime/break-frequency bar plots
* **``structure_bfactors.pdb``** — frame-0 reference structure colored by
  per-residue H-bond variance across replicates (see :ref:`default-outputs` for
  where "frame 0" comes from on each input path)

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

A fresh, timestamped directory under ``./phenom_outputs/`` (or
``$PHENOMS_OUTPUT_DIR``) per run — see :ref:`default-outputs` for the exact
naming and how to override it with ``output_dir=`` / ``--output-dir``, or
disable default writing entirely with ``output_dir=False``.
