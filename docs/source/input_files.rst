Input Files
===========

PhonoMC accepts only sectioned TOML files. Top-level key aliases, legacy text
files, unknown sections, unknown keys, duplicate keys, and invalid scalar
types are rejected.

Canonical sections
------------------

The supported sections are:

- ``[geometry]``
- ``[simulation]``
- ``[boundary]``
- ``[heat_source]``
- ``[io]``

``[geometry]``
--------------

``model``
   ``"box"`` or a path to an STL/OBJ mesh. Relative paths are resolved from
   the working directory first and then from the input-file directory.

``sizes``
   Required for ``model = "box"``. Exactly three positive values ``[x,y,z]``
   in nm.

``merge_coplanar_facets``
   Optional boolean, default ``true``. Connected coplanar triangles are
   treated as one facet when enabled.

Example:

.. code-block:: toml

   [geometry]
   model = "box"
   sizes = [100, 100, 100]
   merge_coplanar_facets = true

``[simulation]``
----------------

``particle_count``
   Positive particle count. Default: ``10000``.

``time_step``
   Positive time step in ps. Default: ``1.0``.

``iterations``
   Positive integer number of time steps. Default: ``10000``.

``random_seed``
   Unsigned 64-bit seed used to derive the initial, injection, and OpenMP
   random streams. Default: ``12345``. Reusing the seed reproduces a run when
   the executable, input, and OpenMP thread configuration are unchanged.

``grid_xyz``
   Required positive integer grid dimensions ``[nx,ny,nz]``.

``initial_temperature``
   Either an unquoted numeric temperature in K or ``"linear"``. Linear mode
   interpolates between the coldest and hottest valid thermal reservoirs.

``compute_kappa``
   Enable conductivity estimates. Automatically disabled when a volumetric
   heat source is enabled.

``convergence_write_interval``
   Positive integer output stride. Default: ``10``.

``profile_timers``
   Write solver timing diagnostics when ``true``.

``progress_temperature_summary_only``
   Print only Tmin/Tavg/Tmax instead of every grid temperature.

``temperature_lookup_dt``
   Positive lookup-table temperature step in K. Default: ``0.1``.

``background_temperature``
   Unquoted, finite, non-negative temperature in K. It is the invariant global
   reference for deviational energy. ``"local"`` is deliberately rejected
   because a moving reference can create gauge energy during carrier motion.
   Use ``0`` for full-phonon sampling or a fixed value such as ``300`` for
   deviational sampling.

``lifetime_temperature``
   Either an unquoted, finite, non-negative fixed temperature in K or
   ``"local"`` to evaluate each mode lifetime at the local particle
   temperature.

``[boundary]``
--------------

Boundary regions use relative coordinates in the domain bounding box:

``[xmin,ymin,zmin,xmax,ymax,zmax]``

Three arrays must have identical lengths:

``boundary_position``
   Region selectors. Every declared region must match at least one mesh facet.

``boundary_conditions``
   One code per region: ``"T"`` thermal reservoir, ``"P"`` periodic, or
   ``"R"`` rough boundary.

``boundary_values``
   Temperature in K for ``T``, roughness in nm for ``R``, and exactly ``0``
   for ``P``.

``periodic_pair``
   Every periodic region must appear exactly once. Entries are grouped into
   pairs, and paired facets must have opposite normals and compatible areas
   and boundary-edge lengths.

Later boundary regions take precedence if selectors overlap. Prefer
non-overlapping selectors unless deliberate precedence is required.

Example:

.. code-block:: toml

   [boundary]
   boundary_position = [
     [-0.01, 0, 0, 0.01, 1, 1],
     [ 0.99, 0, 0, 1.01, 1, 1],
     [0, -0.01, 0, 1, 0.01, 1],
     [0,  0.99, 0, 1, 1.01, 1]
   ]
   boundary_conditions = ["T", "T", "P", "P"]
   boundary_values = [301, 299, 0, 0]
   periodic_pair = [
     [0, -0.01, 0, 1, 0.01, 1],
     [0,  0.99, 0, 1, 1.01, 1]
   ]

``[heat_source]``
-----------------

``enabled``
   Enable a local volumetric energy source. The source changes the occupation
   carried by the existing fixed-weight Monte Carlo carriers; it does not
   create additional computational particles. Energy is distributed over the
   sampled modes as a Bose--Einstein thermal increment, with a cell-wise
   pseudo-temperature solve enforcing the prescribed energy exactly.

``profile``
   ``"uniform"`` or ``"gaussian"``.

``power_density``
   W/m³. For Gaussian sources this is the peak value.

Uniform profile keys
~~~~~~~~~~~~~~~~~~~~

``min`` and ``max`` define opposite corners of an absolute region in nm, the
same unit as ``[geometry].sizes`` for box geometries. For STL/mesh geometries,
use the mesh coordinate system expressed in nm.

.. code-block:: toml

   [heat_source]
   enabled = true
   profile = "uniform"
   min = [6.6, 21.6, 21.0]
   max = [15.4, 28.8, 30.0]
   power_density = 1.0e20

Gaussian profile keys
~~~~~~~~~~~~~~~~~~~~~

``center`` and ``sigma`` are absolute lengths in nm, the same unit as
``[geometry].sizes`` for box geometries. For STL/mesh geometries, use the
mesh coordinate system expressed in nm. A non-positive sigma makes that axis
uniform.

.. code-block:: toml

   [heat_source]
   enabled = true
   profile = "gaussian"
   center = [11.0, 25.2, 30.0]
   sigma = [0.0, 1.44, 3.0]
   power_density = 1.0e20

``[io]``
--------

``material_folder``
   Path to a directory containing a valid POSCAR and HDF5 material file.

``output_folder``
   Base result path. PhonoMC appends ``_0``, ``_1``, and so on.

.. code-block:: toml

   [io]
   material_folder = "../material/Si/"
   output_folder = "results/Cross_100nm/"

Strict validation summary
-------------------------

- use only literal TOML booleans ``true`` and ``false``
- do not quote numeric ``initial_temperature`` values
- use integer, positive ``grid_xyz`` dimensions
- keep all boundary arrays aligned
- use finite, non-negative thermal and roughness values
- assign every periodic facet to exactly one valid pair
- provide valid POSCAR and HDF5 material data
