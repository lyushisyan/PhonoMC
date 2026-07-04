Reference
=========

Command-line interface
----------------------

.. code-block:: text

   PhonoMC [input.toml]

With no argument, the executable searches for ``input.toml`` in the current
directory and near the executable. Only sectioned TOML input is supported.

Boundary codes
--------------

``T``
   Thermal reservoir. ``boundary_values`` is temperature in K. Particles that
   reach the facet are absorbed, and reservoir particles are injected using
   the one-to-one refill bookkeeping.

``P``
   Periodic boundary. The value must be ``0``. The paired facet determines the
   translation applied to the particle position.

``R``
   Rough boundary. The value is RMS roughness in nm. Reflection uses the
   precomputed mode-selection data and specularity model.

Output files
------------

``convergence.txt``
   Space-separated time history. Columns include ``timestep``, ``time_ps``,
   grid temperatures, ``heatflux``, ``kappa_int``, ``kappa_eff``, ``absorbed``,
   ``injected``, ``recovered``, and ``net``.

``summary.txt``
   Normalized input, geometry bounds, volume, facet counts, boundary counts,
   grid information, runtime, OpenMP settings, and optional profiling data.

``grid_centers.csv``
   ``x_nm,y_nm,z_nm,volume_nm3`` for every retained control-volume center.

``rough_boundary_mode_map.csv``
   Written when rough-facet preprocessing data is available. Contains incoming
   and reflected mode mapping diagnostics.

Plotting commands
-----------------

.. code-block:: bash

   python3 tools/plot_convergence.py RESULT_DIRECTORY

.. code-block:: bash

   python3 tools/plot_temperature_3d.py \
     --input INPUT.toml --results RESULT_DIRECTORY

Developer modules
-----------------

- ``src/SimulationConfig.cpp``: strict TOML parsing and validation
- ``src/SurfaceMesh.cpp``: geometry topology and spatial queries
- ``src/SimulationDomain.cpp``: boundary and grid construction
- ``src/PhononMaterial.cpp``: POSCAR/HDF5 material loading
- ``src/MonteCarloSolver.cpp``: particle transport and observables

Tests
-----

.. code-block:: bash

   ctest --test-dir build --output-on-failure

The CTest targets cover configuration, domain/boundary behavior, the provided
geometry examples, and strict material loading.
