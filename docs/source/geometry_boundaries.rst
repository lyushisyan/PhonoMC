Geometry and Boundaries
=======================

Box geometry
------------

For a rectangular domain:

.. code-block:: toml

   [geometry]
   model = "box"
   sizes = [1000, 100, 100]

The three sizes are positive lengths in nm. The box is triangulated internally
and coplanar triangles are merged into six facets by default.

Mesh geometry
-------------

``model`` may point to an ASCII/binary STL or OBJ file. STL coordinates are
interpreted as nm and converted internally to angstroms. The loaded mesh is
shifted so its bounding-box minimum is at the origin.

A production mesh should be:

- closed and manifold
- free of zero-area triangles
- consistently scaled in nm
- simple enough that facet merging and ray intersection remain robust

``merge_coplanar_facets = true`` groups connected coplanar triangles. Boundary
selectors operate on the resulting facet centroids, not on every triangle
vertex.

Relative boundary selectors
---------------------------

Each boundary region is a six-value relative box:

.. code-block:: text

   [xmin, ymin, zmin, xmax, ymax, zmax]

``0`` and ``1`` correspond to the geometry bounding-box minimum and maximum.
Small extensions such as ``-0.01`` and ``1.01`` make boundary-plane selectors
robust to floating-point tolerances.

Every declared selector must match at least one facet. The three arrays
``boundary_position``, ``boundary_conditions``, and ``boundary_values`` must
have identical lengths. When selectors overlap, the later selector determines
the final condition.

Thermal reservoirs
------------------

``T`` values are non-negative temperatures in K. A particle that reaches a
thermal facet is marked absorbed, and the facet's leaving count contributes to
the next reservoir-refill step.

Use at least two reservoirs with different temperatures for a conductivity
calculation. Linear initialization chooses the coldest and hottest valid
reservoirs and interpolates along their centroid-to-centroid direction.

Periodic boundaries
-------------------

``P`` values must be ``0``. Every periodic selector must appear exactly once
in ``periodic_pair`` and the total number of regions must be even.

For each pair, PhonoMC validates:

- both selected facets are marked periodic
- the two regions match the same number of facets
- normals point in opposite directions
- areas agree within tolerance
- boundary-edge length spectra agree
- no facet is assigned to two different partners

The centroid difference defines the translation from one periodic facet to
the other.

Rough boundaries
----------------

``R`` values are non-negative RMS roughness values in nm. Unselected facets
default to ``R`` with zero roughness. Explicitly list important surfaces rather
than relying on this default in production inputs.

Grid construction
-----------------

``grid_xyz`` divides the bounding box into regular candidate cells. For a box,
all centers are retained and position-to-grid indexing is direct. For a mesh,
only centers inside the surface are retained. Initialization fills a complete
voxel-to-active-grid table from those cells, so particle positions use an
``O(1)`` lookup for both boxes and arbitrary closed meshes.

Every 100 time steps, mesh particles are checked against the actual closed
surface rather than only its axis-aligned bounding box. A failed boundary trace
is retried from a nearby interior point without changing velocity; a persistent
failure stops the run with the particle position and an STL validation hint.
Likewise, more than 64 boundary collisions within one time step is reported as
an input/time-step error instead of silently dropping the remaining travel
time.

Current limitation: retained mesh cells receive the same volume
``domain volume / retained cell count``. Cells cut by a curved or oblique
surface are not assigned exact intersection volumes. Perform a grid-sensitivity
study for strongly non-rectangular devices.

The domain volume itself is evaluated from the oriented surface using the
signed tetrahedral-volume identity. For topologically closed meshes,
inconsistent triangle winding is repaired once during loading. Cached volume
tetrahedra are accepted only when their summed volume agrees with the surface
volume; otherwise sampling falls back to the robust surface test.

Geometry validation workflow
----------------------------

Before a long mesh run:

1. use a coarse grid and small particle count
2. verify the reported face and facet counts
3. inspect ``grid_centers.csv``
4. check thermal-reservoir, rough-facet, and periodic-pair counts in
   ``summary.txt``
5. render the mesh and temperature slices with ``plot_temperature_3d.py``
