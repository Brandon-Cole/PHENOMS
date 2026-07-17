Usage
=====

Output location
---------------

Artifacts default to ``./phenom_outputs/`` next to your working directory.
Override globally:

.. code-block:: bash

   export PHENOMS_OUTPUT_DIR=/path/to/my_phenom_runs

Or pass ``output_dir=`` / ``--output-dir`` per run.

Simple form: multi-frame PDBs
-----------------------------

This is the recommended path when you already have protein-centered,
fitted multi-frame PDBs. Backbone N–O analysis is the default
(``backbone_only=True``):

.. code-block:: python

   from phenoms import SimulationSet, default_output_root

   sim = SimulationSet(
       pdb_files=["rep1.pdb", "rep2.pdb"],
       resid_range=(50, 70),       # plot filter only; None = whole protein
       sub_frames=100,             # None = all frames
       backbone_only=True,
       output_dir=default_output_root() / "my_run" / "set_a",
   )
   sim.run()
   # raw_data/*_hbonds.csv, *_occupancy.csv, *_pivot.csv, manifest.json

One-shot helper:

.. code-block:: python

   from phenoms import run_backbone_hbond_analysis

   sim = run_backbone_hbond_analysis(
       ["r1.pdb", "r2.pdb"],
       sub_frames=100,
       plot_heatmaps=False,
       output_dir=default_output_root() / "quick",
   )

Native trajectories
-------------------

Load engine-native files without writing intermediate PDBs:

.. code-block:: python

   from phenoms import SimulationSet

   sim = SimulationSet.from_trajectories(
       ["rep1.xtc", "rep2.xtc"],
       topology="system.pdb",          # or topologies=[...]
       sub_frames=200,
       backbone_only=True,
   )
   sim.run()

Equivalent constructor form:

.. code-block:: python

   SimulationSet(
       trajectories=["rep1.xtc", "rep2.xtc"],
       topology="system.pdb",
   )

For PBC imaging / centering / fitting first, use
:func:`phenoms.prepare_set_from_dir` or ``phenoms prep``.

All-bond mode
-------------

Backbone N–O remains the HDX-style default. Opt into all donor/acceptor classes:

.. code-block:: python

   sim = SimulationSet(["rep1.pdb"], backbone_only=False, sub_frames=100)
   sim.run()

Comparison sets
---------------

:class:`~phenoms.ComparisonSet` may be constructed before ``.run()``; methods
that need pivots validate at call time:

.. code-block:: python

   from phenoms import SimulationSet, ComparisonSet, default_output_root

   base = default_output_root() / "my_run"
   a = SimulationSet(["a1.pdb", "a2.pdb"], resid_range=(50, 70), sub_frames=100)
   b = SimulationSet(["b1.pdb", "b2.pdb"], resid_range=(50, 70), sub_frames=100)
   cmp = ComparisonSet(a, b, label_a="apo", label_b="holo")
   a.run()
   b.run()
   cmp.compare()
   cmp.export_comparison_artifacts(base / "comparison")
   cmp.plot_difference()
   # cmp.export_connectivity_graph_html("network.html")

Engine folders (GROMACS / OpenMM / AMBER)
-----------------------------------------

Directory helpers discover replicates, normalize to protein-only multi-frame
PDBs, then return ready ``SimulationSet`` objects:

.. code-block:: python

   from phenoms import simulation_set_from_dir, comparison_sets_from_dirs

   sim = simulation_set_from_dir(
       input_dir="./wt_set",
       prepared_dir="./prepared/wt",
       frame_dt_ps=1000,
       start_ps=1000,
       end_ps=500000,
   ).run()

   set_a, set_b, comp = comparison_sets_from_dirs(
       dir_a="./wt_set",
       dir_b="./mut_set",
       prepared_dir_a="./prepared/wt",
       prepared_dir_b="./prepared/mut",
       label_a="wt",
       label_b="mut",
   )
   set_a.run()
   set_b.run()
   comp.compare()

Required files per replicate directory:

* **GROMACS**: ``.xtc``/``.trr`` + topology ``.tpr`` (preferred) or ``.gro``/``.pdb``
* **OpenMM**: ``.dcd``/``.xtc`` + ``.pdb``/``.prmtop``
* **AMBER**: ``.nc``/``.mdcrd`` + ``.prmtop``/``.parm7``

Direct detection helpers
------------------------

Skip the ``SimulationSet`` wrapper when you only need bond tables:

.. code-block:: python

   from phenoms import detect_hbonds, detect_hbonds_with_occupancy

   df_bb = detect_hbonds("trajectory.pdb", backbone_only=True)
   df = detect_hbonds("rep.xtc", top="system.pdb")
   hbonds_df, occ_df = detect_hbonds_with_occupancy(
       "trajectory.pdb",
       output_csv_path="hbond_occupancy.csv",
   )

Quality control
---------------

Optional fail-fast QC inside ``SimulationSet.run``:

.. code-block:: python

   sim = SimulationSet(["rep1.pdb", "rep2.pdb"], sub_frames=250)
   sim.run(
       qc=True,
       mdp_files=["rep1.mdp", "rep2.mdp"],  # optional GROMACS consistency
       qc_fail_on_nonconverged=True,
   )
   report = sim.get_qc_report()

Checks include:

* RMSD window-shift / drift convergence per replicate
* MDP key consistency when ``mdp_files`` are provided

Renumbering & GROMACS export scripts
------------------------------------

For mismatched residue numbering across mutants/constructs:

.. code-block:: bash

   python scripts/renumber_many_to_reference.py \
     --reference ref.pdb --mobile-dir ./exports \
     --mobile-glob "*.pdb" --output-dir ./exports_renumbered

Standalone GROMACS → PDB helper:

.. code-block:: bash

   python scripts/preprocess_gromacs_to_pdb.py \
     --run-dir /path/to/run --out-dir /path/to/exports --frame-dt-ps 1000
