Installation
============

Conda environment
-----------------

Install `Rust <https://rustup.rs/>`_ first — the build is driven by
`maturin <https://www.maturin.rs/>`_, which compiles the Rust geometry kernel
(``phenoms.phenoms_hbond_rs``) into the same wheel as the Python package.

.. code-block:: bash

   conda env create -f environment.yml
   conda activate phenoms
   pip install -e .

The editable install builds the Rust extension and registers the ``phenoms``
CLI entry point in one step. Detection uses the Rust geometry kernel when it's
built; **Polars** speeds occupancy aggregation only (not the Baker–Hubbard
kernel itself). Without a Rust toolchain at install time, PHENOMS falls back
to MDTraj for detection.

For local iteration on the Rust source, ``maturin develop --release`` from the
repo root rebuilds just the extension without reinstalling dependencies.

Using uv
--------

.. code-block:: bash

   uv venv
   source .venv/bin/activate
   uv pip install -e .

Optional extras
---------------

.. code-block:: bash

   pip install -e ".[gromacs]"     # MDAnalysis to read GROMACS .tpr topologies
   pip install -e ".[benchmark]"   # MDAnalysis for kernel benchmarks

BioPython (PDB B-factor writes via :func:`phenoms.write_pdb_bfactors`) is a
core dependency and needs no extra.

Verify
------

.. code-block:: bash

   python -c "import phenoms; print(phenoms.__version__)"
   phenoms --help
