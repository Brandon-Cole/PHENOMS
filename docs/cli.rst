Command-line interface
======================

The Python API is primary. After ``pip install -e .``, a Click CLI is also available:

.. code-block:: bash

   phenoms --help
   python -m phenoms --help

``phenoms prep``
----------------

Normalize GROMACS / OpenMM / AMBER replicate folders to multi-frame PDBs:

.. code-block:: bash

   phenoms prep --input-dir ./wt_set --prepared-dir ./prepared/wt

``phenoms run``
---------------

Analyze PDB replicates or native trajectories:

.. code-block:: bash

   phenoms run --pdb rep1.pdb --pdb rep2.pdb --sub-frames 100 --output-dir ./out
   phenoms run --traj rep1.xtc --traj rep2.xtc --topology top.pdb
   phenoms run --pdb rep1.pdb --all-bonds   # opt into all-bond mode

``phenoms compare``
-------------------

Prep two class directories, run analysis, and export comparison artifacts:

.. code-block:: bash

   phenoms compare \
     --dir-a ./wt --dir-b ./mut \
     --prepared-dir-a ./prepared/wt --prepared-dir-b ./prepared/mut \
     --label-a wt --label-b mut
