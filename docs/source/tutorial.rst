Tutorial
========

This tutorial uses the repository-provided examples. Build PhonoMC before
starting.

Cross-plane conductivity
------------------------

Run:

.. code-block:: bash

   OMP_NUM_THREADS=8 ./build/PhonoMC example/input_cross_100nm.toml

This example uses two thermal reservoirs on the x faces and periodic
boundaries on y and z. ``initial_temperature = "linear"`` initializes a
profile between 299 K and 301 K, and ``compute_kappa = true`` enables the
endpoint and fitted conductivity estimates.

After the run:

.. code-block:: bash

   python3 tools/plot_convergence.py example/results/Cross_100nm_1

Inspect ``kappa_vs_time.png`` together with particle balance and temperature
convergence; a flat conductivity curve alone is not sufficient validation.

In-plane transport with rough surfaces
---------------------------------------

Run:

.. code-block:: bash

   OMP_NUM_THREADS=8 ./build/PhonoMC \
     example/input_inplane_x1000nm_z100nm_r1nm.toml

The x faces are thermal reservoirs, the y faces are periodic, and the z faces
have 1 nm roughness. Compare this result with the cross-plane case to study
geometry and surface-scattering effects.

FinFET heating
--------------

The FinFET example loads ``model/finfet.stl`` and injects a Gaussian
volumetric heat source using absolute nm coordinates in the mesh coordinate
system:

.. code-block:: bash

   OMP_NUM_THREADS=8 ./build/PhonoMC \
     example/input_finfet_stl_heat1e20.toml

The repository input is a lightweight tutorial/smoke-test case
(``particle_count = 200000``, ``time_step = 0.05`` ps, and
``iterations = 2000``). For production calculations, copy the input and then
increase ``particle_count``, ``iterations``, and, if needed, ``grid_xyz`` after
checking the generated ``summary.txt``.

Plot an existing or completed result:

.. code-block:: bash

   python3 tools/plot_temperature_3d.py \
     --input example/input_finfet_stl_heat1e20.toml \
     --results example/results/FinFET_heat1e20_1

The command writes a 3D temperature view and x/y slice images under
``plots_3d``.

Create a new case
-----------------

1. Copy the closest example TOML file.
2. Select ``model = "box"`` or a valid STL/OBJ path.
3. Define all boundary arrays with identical lengths.
4. Use ``0`` for every periodic boundary value.
5. Pair every periodic region exactly once.
6. Point ``material_folder`` to a directory containing POSCAR and HDF5 data.
7. Start with a small grid and particle count.
8. Inspect the generated summary before scaling the calculation.
