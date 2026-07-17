PHENOMS documentation
=====================

**PHENOMS** (Python-based Hydrogen-Deuterium Exchange of Molecular Dynamics Simulations)
analyzes backbone hydrogen-bond networks from MD trajectories for HDX-MS–style interpretation.

.. image:: ../assets/phenoms.png
   :alt: PHENOMS workflow
   :align: center
   :width: 85%

Quick links
-----------

* :doc:`installation` — conda / pip setup and optional Rust acceleration
* :doc:`usage` — Python API for PDB and native trajectories
* :doc:`cli` — optional Click CLI (``phenoms``)
* :doc:`api/index` — API reference

At a glance
-----------

.. code-block:: python

   from phenoms import SimulationSet

   sim = SimulationSet(
       pdb_files=["rep1.pdb", "rep2.pdb"],
       sub_frames=100,
       backbone_only=True,  # default (HDX-style N–O)
   )
   sim.run()

.. toctree::
   :maxdepth: 2
   :hidden:

   installation
   usage
   cli
   api/index
