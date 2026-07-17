Usage
=====

Outputs
-------

By default, artifacts go under ``./phenom_outputs/``. Override with:

.. code-block:: bash

   export PHENOMS_OUTPUT_DIR=/path/to/my_phenom_runs

Simple form: multi-frame PDBs
-----------------------------

Backbone N–O analysis is the default (``backbone_only=True``):

.. code-block:: python

   from phenoms import SimulationSet, default_output_root

   sim = SimulationSet(
       pdb_files=["rep1.pdb", "rep2.pdb"],
       resid_range=(50, 70),
       sub_frames=100,
       backbone_only=True,
       output_dir=default_output_root() / "my_run" / "set_a",
   )
   sim.run()

Native trajectories
-------------------

Load engine-native files without writing intermediate PDBs:

.. code-block:: python

   from phenoms import SimulationSet

   sim = SimulationSet.from_trajectories(
       ["rep1.xtc", "rep2.xtc"],
       topology="system.pdb",
       sub_frames=200,
   )
   sim.run()

For PBC imaging / centering / fitting first, use
:func:`phenoms.prepare_set_from_dir` or ``phenoms prep``.

All-bond mode
-------------

Opt in when you need all donor/acceptor classes:

.. code-block:: python

   sim = SimulationSet(["rep1.pdb"], backbone_only=False, sub_frames=100)
   sim.run()

Comparison sets
---------------

:class:`~phenoms.ComparisonSet` can be constructed before ``.run()``:

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

Engine folders (GROMACS / OpenMM / AMBER)
-----------------------------------------

.. code-block:: python

   from phenoms import simulation_set_from_dir, comparison_sets_from_dirs

   sim = simulation_set_from_dir(
       input_dir="./wt_set",
       prepared_dir="./prepared/wt",
       frame_dt_ps=1000,
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

Direct detection helpers
------------------------

.. code-block:: python

   from phenoms import detect_hbonds, detect_hbonds_with_occupancy

   df_bb = detect_hbonds("trajectory.pdb", backbone_only=True)
   df = detect_hbonds("rep.xtc", top="system.pdb")
   hbonds_df, occ_df = detect_hbonds_with_occupancy(
       "trajectory.pdb",
       output_csv_path="hbond_occupancy.csv",
   )
