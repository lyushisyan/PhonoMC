Performance
===========

Release builds
--------------

Use ``CMAKE_BUILD_TYPE=Release`` for production simulations. Debug builds are
valuable for diagnosis but can be dramatically slower in particle and geometry
loops.

OpenMP
------

Enable OpenMP at configure time and select the thread count at runtime:

.. code-block:: bash

   cmake -S . -B build \
     -DCMAKE_BUILD_TYPE=Release \
     -DPHONOMC_ENABLE_OPENMP=ON
   cmake --build build -j
   OMP_NUM_THREADS=32 ./build/PhonoMC input.toml

Do not assume that the largest available thread count is fastest. Compare a
short representative run at 1, 2, 4, 8, and higher thread counts. Memory
bandwidth, dynamic collision scheduling, grid reductions, and the active-mode
table size can limit scaling.

Built-in profiling
------------------

Enable:

.. code-block:: toml

   [simulation]
   profile_timers = true

The runtime summary separates major stages such as main particle advance,
particle removal, injection construction, collision-cache updates, injected
particle advance, temperature update, lifetime scattering, and statistics.
``update_temp`` includes every temperature refresh in the timestep (up to three
with a source), but excludes initialization. ``source`` includes both source
halves, ``boundary_checks`` covers boundary-counter merging and periodic domain
validation, and ``energy_ledger`` covers final energy bookkeeping. ``other``
reports remaining orchestration time; stage times are exclusive, not nested.
Progress printing is outside ``profile_total_seconds``.

``cell_index_rebuilds`` and ``temperature_refreshes`` include initialization.
Normally the former is ``iterations + 1``: occupation-only changes reuse counts
and stable particle buckets, while energy and temperature are still refreshed.
Fixed-background Bose occupations and material energy weights are cached.
There is no new duplicate angular-frequency table: material frequencies are
already stored in a lookup array. Speedups must be measured on a representative
production case, not inferred from regression-test durations.

Scaling parameters
------------------

``particle_count``
   Usually the largest direct memory and transport-work multiplier. Increase
   it only after a low-count case is geometrically correct.

``iterations``
   Controls total simulated time and total transport work.

``grid_xyz``
   Increases temperature/flux reduction buffers and output width. For mesh
   geometries, initialization also builds a full bounding-box voxel map. The
   runtime position-to-grid lookup remains constant time.

``temperature_lookup_dt``
   Smaller values create larger temperature-energy tables and increase
   initialization work, while improving interpolation resolution.

``convergence_write_interval``
   Larger values reduce file I/O and conductivity-statistics frequency.

Memory considerations
---------------------

Particle state is stored in parallel arrays for modes, positions, velocities,
temperatures, occupations, energy, grid IDs, collision locations, collision
facets, boundary conditions, and alive flags. Ten million particles therefore
require substantially more memory than the position array alone suggests.

OpenMP grid reductions allocate per-thread buffers proportional to
``thread count × grid count``. Very fine grids combined with many threads can
consume significant additional memory.

Recommended scaling workflow
----------------------------

1. Validate input with a very small case.
2. Choose a grid using a grid-sensitivity study.
3. Select a stable time step.
4. Increase particles until noise is acceptable.
5. Increase iterations until observables are stationary.
6. Benchmark thread count on the final problem shape.
7. Enable profiling for at least one representative production run.

Known performance hotspots
--------------------------

- boundary intersection and collision processing
- rough-boundary mode selection
- particle-array compaction after absorption
- temperature and heat-flux grid reductions
- initialization of detailed mesh-volume sampling data

Surface winding repair and volume-consistency checks are initialization-only,
linear passes over existing faces or tetrahedra. They do not add work to the
per-particle time-stepping path.

Axis-aligned ``box`` domains sample initial positions directly from three
uniform coordinate distributions. This avoids general volume decomposition
for high-aspect-ratio films and ribbons. The position distribution is unchanged,
but a fixed seed produces different particle trajectories than the older
tetrahedral sampler.

Reproducibility
---------------

Set ``simulation.random_seed`` explicitly for production comparisons. OpenMP
threads receive independently derived streams; therefore exact replay also
requires the same thread count and scheduling environment.
