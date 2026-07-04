Getting Started
===============

Run the cross-plane example
---------------------------

From the repository root:

.. code-block:: bash

   ./build/PhonoMC example/input_cross_100nm.toml

PhonoMC creates an indexed directory such as
``example/results/Cross_100nm_1`` so that an existing result is not
overwritten.

Select the OpenMP thread count
------------------------------

.. code-block:: bash

   OMP_NUM_THREADS=16 ./build/PhonoMC example/input_cross_100nm.toml

Start from a copied input
-------------------------

.. code-block:: bash

   cp example/input_cross_100nm.toml input.toml
   ./build/PhonoMC input.toml

If no argument is supplied, PhonoMC searches for ``input.toml`` near the
current working directory and executable.

Inspect the first results
-------------------------

The result directory contains at least:

- ``summary.txt``: normalized input, geometry, boundary, and runtime details
- ``grid_centers.csv``: grid-center coordinates and cell volumes
- ``convergence.txt``: temperature and transport observables versus time

Plot a one-dimensional run
--------------------------

.. code-block:: bash

   python3 tools/plot_convergence.py example/results/Cross_100nm_1

The command writes temperature, heat-flux, and thermal-conductivity figures
under ``plots_1d`` in the selected result directory.

Validate before a long run
--------------------------

Before increasing particle count or iteration count:

1. run with a small particle count
2. confirm every boundary region matches the intended facet
3. inspect ``summary.txt`` and ``grid_centers.csv``
4. check particle balance in ``convergence.txt``
5. then scale threads, particles, grid size, and iterations
