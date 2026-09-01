Benchmarks
==========

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
