Troubleshooting
===============

Unknown TOML key or section
---------------------------

PhonoMC intentionally rejects aliases and misspelled keys. Move the key into
one of the five documented sections and compare it with :doc:`input_files`.

Boundary array length error
---------------------------

``boundary_position``, ``boundary_conditions``, and ``boundary_values`` must
have the same number of entries. Periodic entries still require a placeholder
value of ``0``.

Boundary region did not match a facet
-------------------------------------

Selectors are relative to the geometry bounding box and match facet
centroids. Inspect the mesh bounds and widen the thin selector interval, for
example from ``0.0`` to ``[-0.01,0.01]`` around a boundary plane.

Periodic boundary validation failed
-----------------------------------

Check that:

- both regions are marked ``P``
- both boundary values are ``0``
- every periodic region occurs exactly once in ``periodic_pair``
- the pair matches the same number of facets
- normals are opposite and areas/edge spectra are compatible

Strict POSCAR loading failed
----------------------------

The selected ``material_folder`` must contain a readable POSCAR. PhonoMC no
longer substitutes a unit-cell volume. Fix the path or POSCAR contents rather
than bypassing the error.

HDF5 loading failed
-------------------

Verify that ``kappa-fbz.hdf5`` or ``kappa.hdf5`` exists and contains the
required datasets with compatible dimensions. Use an absolute
``material_folder`` temporarily when diagnosing relative paths.

``qdelaunay`` is missing
------------------------

Install Qhull command-line tools:

.. code-block:: bash

   sudo apt install qhull-bin

OpenMP was not found
--------------------

The project still builds serially. Install the compiler OpenMP runtime and
reconfigure from a clean build directory. On Ubuntu:

.. code-block:: bash

   sudo apt install libomp-dev

Particle balance drifts
-----------------------

Inspect ``absorbed``, ``injected``, ``recovered``, and ``net`` columns. Confirm
periodic pairing, reservoir selection, time step, and mesh closure before
increasing particle count.

The 3D run is too slow
----------------------

First reduce ``particle_count``, ``iterations``, and ``grid_xyz``. Confirm the
geometry and source region using a smoke test, then scale one dimension at a
time and enable OpenMP.
