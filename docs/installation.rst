Installation
============

Conda environment
-----------------

Use the project ``environment.yml``:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate phenoms
   pip install -e .

Rust extension (recommended)
----------------------------

For faster Baker–Hubbard detection:

.. code-block:: bash

   pip install maturin
   maturin develop --release

Detection uses the Rust backend when available; Polars accelerates occupancy aggregation only.

Optional extras
---------------

.. code-block:: bash

   pip install -e ".[benchmark]"   # MDAnalysis for kernel benchmarks
   pip install -e ".[structure]"   # BioPython for PDB B-factor writes
   pip install -e ".[docs]"        # Sphinx + Shibuya for local docs builds

Build these docs locally
------------------------

.. code-block:: bash

   pip install -e ".[docs]"
   cd docs
   make html
   # open _build/html/index.html
