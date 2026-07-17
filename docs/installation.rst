Installation
============

Conda environment
-----------------

.. code-block:: bash

   conda env create -f environment.yml
   conda activate phenoms
   pip install -e .

The editable install registers the ``phenoms`` CLI entry point.

Rust extension (recommended)
----------------------------

Install `Rust <https://rustup.rs/>`_, then from the repo root:

.. code-block:: bash

   pip install maturin
   maturin develop --release

This builds ``phenoms_hbond_rs``. Detection uses the Rust geometry kernel when
available; **Polars** speeds occupancy aggregation only (not the Baker–Hubbard
kernel itself). Without Rust, PHENOMS falls back to MDTraj.

Optional extras
---------------

.. code-block:: bash

   pip install -e ".[benchmark]"   # MDAnalysis for kernel benchmarks
   pip install -e ".[structure]"   # BioPython for PDB B-factor writes

Verify
------

.. code-block:: bash

   python -c "import phenoms; print(phenoms.__version__)"
   phenoms --help
