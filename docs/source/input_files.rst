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
- ``[scattering]``
- ``[driving]``
- ``[sampling]``
- ``[convergence]``
- ``[boundary]``
- ``[materials]``
- ``[heat_source]``
- ``[mode_excitation]``
- ``[io]``

``[scattering]``
----------------

``model``
   ``"rta"`` (default) or ``"callaway"``. Both collision paths update existing
   carrier occupations without changing carrier count. Callaway uses
   fixed-reference N/R dual relaxation (R = U plus supplied isotope, impurity
   and defect linewidths), with discrete energy and N momentum conservation;
   see :doc:`callaway`. ``rta_linearized`` and ``full_matrix`` are no longer
   public input selections.

``matrix_files``
   Forbidden for the public RTA and Callaway models. The matrix adapter is
   retained as an internal experimental API; see :doc:`full_scattering`.

``[driving]``
-------------

``temperature_gradient``
   Optional signed nonzero imposed temperature gradient in K/m along
   ``simulation.transport_axis``. Enables periodic linear-response transport
   with RTA or Callaway. Requires paired periodic faces in the drive direction,
   no thermal reservoirs, and ``initial_temperature = lifetime_temperature =
   background_temperature > 0``. See :doc:`temperature_gradient` and
   ``example/input_inplane_gradient.toml`` for the complete film setup.

``[sampling]`` and ``[convergence]``
------------------------------------

Existing-carrier RTA/Callaway requires ``sampling.resample_interval = 0``.
Optional gradient warm starts (``driving.warm_start_ps``) and stationarity
screens (``convergence.enabled`` and window/tolerance settings) are documented
in :doc:`temperature_gradient`. They are not available for ordinary
reservoir-driven or locally heated runs. A passed stationarity screen does not
establish particle-count or time-step convergence.

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

``transport_axis``
   Optional explicit heat-flow axis, ``"x"``, ``"y"`` or ``"z"``. When omitted,
   legacy runs infer it from the largest grid dimension. Required with
   ``[driving] temperature_gradient`` (signed, nonzero K/m). Gradient driving
   uses periodic boundaries and fixed-reference linearized RTA or Callaway;
   see :doc:`temperature_gradient` for constraints and the output definition.

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

Select ``scattering.model = "rta"`` or ``"callaway"``. Ordinary RTA uses its
existing Bose pseudo-temperature collision update and optional local lifetime
lookup. Callaway is linearized at ``background_temperature``. Use small
temperature perturbations and matched fixed rates when comparing their N
scattering treatment. Gradient-driven RTA is linearized automatically.

``enabled``
   Enable a local volumetric source that changes existing carrier occupations.
   RTA thermal injection uses a Bose--Einstein increment; Callaway uses the
   reference heat-capacity weights. Both conserve the deposited energy.
   Empty cells, or cells without eligible sampled modes, defer energy locally.
   No source carriers are created.

``profile``
   ``"uniform"`` or ``"gaussian"``.

``power_density``
   W/m³. For Gaussian sources this is the peak value.

``total_power``
   Optional total lattice-heating power in W, mutually exclusive with
   ``power_density`` and ``power_densities``. The spatial profile is normalized
   over the selected simulation grid using cell volumes. Multiple Gaussian
   centers share this one total budget; this is not the power of each center.
   The time envelope and ``amplitude`` multiply this base power. Region selection
   still uses grid centers, so a partially intersected cell is not integrated
   geometrically. An enabled source selecting no grid centers is an error.

``spectrum``
   ``"thermal"`` (default): a Bose--Einstein occupation increment over eligible
   sampled modes. ``"weighted"``: distribute cell energy in proportion to the
   configured branch weights of eligible carriers, then divide each carrier's
   allocated energy by its phonon energy and statistical weight. Both are
   prescribed, phenomenological spectra, not electron-phonon calculations.

``branches``
   Optional unique **zero-based** branch indices, intersected with the frequency
   window. Indices follow the loaded HDF5 branch ordering; no acoustic/optical
   classification or branch tracking through crossings is inferred.

``spectrum = "gaussian"``
   Gaussian **energy weight per sampled carrier**, supported with ``rta`` and
   ``callaway``. Set ``frequency_center`` (finite, nonnegative) and
   ``frequency_sigma`` (finite, positive standard deviation) in THz,
   using f = omega/(2 pi). They are independent of the spatial ``profile``,
   ``center`` and ``sigma``. Branch/frequency filters intersect the Gaussian;
   ``branch_weights`` cannot be combined with it.

   For eligible carriers in each cell, use
   ``w_i = exp(-0.5*((f_i-frequency_center)/frequency_sigma)^2)`` and
   ``dE_i = P_cell * dt * w_i / sum_i(w_i)``. Repeated modes count as repeated
   quadrature samples. The prescribed cell energy is exact, but the realized
   spectrum has sampling error; it is no longer independent of carrier counts.
   Missing modes are not inserted. Cells with no eligible carriers retain a
   pending source budget. Check spectral sampling and pending energy when
   using narrow bands. With uniform mode sampling the histogram approaches
   DOS times Gaussian, not a prescribed Gaussian density per THz.

   Log-space normalization prevents all-zero underflow. Explicit
   ``frequency_min/max`` can truncate the spectrum; there is no implicit cutoff.
   Filters selecting no eligible material modes are rejected at initialization.
   This is a phenomenological lattice source, not an electron-phonon calculation.
   ``summary.txt`` records centre, width and existing-carrier normalization.
   See ``example/input_gaussian_spectrum.toml``.

``branch_weights``
   Required for ``spectrum = "weighted"``; otherwise omitted. Supply one finite
   nonnegative value per branch, with at least one positive selected weight.
   The array length must match every loaded material's branch count. Values are
   relative **energy weights per sampled mode/carrier**, not occupation weights
   and not prescribed total fractions for entire branches. The realized branch
   totals therefore depend on the number of eligible carriers. Zero excludes a
   branch. Material-specific and fully q-resolved tables are not implemented.

Source energy accounting
~~~~~~~~~~~~~~~~~~~~~~~~

This input prescribes energy delivered to the lattice; it does not solve for
electrons or infer actual device dissipation. Positive heating only is supported.

Every cell separately tracks requested, deposited and pending energy:
``requested = deposited + pending`` within floating-point accuracy. If a cell
has no eligible carriers, its budget remains there; it is never redistributed
to another cell. Pending energy can be deposited after the pulse ends when
eligible carriers arrive. It is **not** included in phonon energy or temperature
before deposition. This deferred-carrier policy applies to both RTA and
Callaway. A significant final pending budget indicates insufficient
sampling or an unsupported modal selection, not a completed heating calculation.

The latest per-cell ledger is written to ``heat_source_ledger.txt`` at convergence
output and final summary. ``convergence.txt`` retains the original columns and
appends ``hs_prescribed_total_ev``, ``hs_pending_total_ev`` and
``hs_ledger_residual_ev``. Existing ``hs_injected_energy_ev`` remains the energy
actually deposited in that step. ``summary.txt`` records the input type, spectrum,
filters, branch weights and the assumption ``electron_transport = false``.

See ``example/input_box_total_power.toml`` for a complete illustrative input.

``frequency_min`` and ``frequency_max``
   Optional modal filter for targeted phonon excitation. Values are specified
   in THz for :math:`\omega/2\pi`. If either cutoff is omitted, that side of
   the frequency window is left open. If both are omitted, the heat source
   uses all sampled modes as in the original all-mode occupation update.

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

Targeted modal excitation example
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   [heat_source]
   enabled = true
   profile = "uniform"
   min = [0, 0, 0]
   max = [100, 100, 100]
   power_density = 1.0e18
   frequency_min = 2.0
   frequency_max = 5.0

``[mode_excitation]``
---------------------

This optional section applies an initial over-occupation to modes in a
specified frequency window while leaving the Monte Carlo mode-sampling
normalization unchanged. It is useful for transient tests of targeted phonon
excitation under ordinary boundary-driven transport.

``enabled``
   Enable the initial modal occupation multiplier.

``frequency_min`` and ``frequency_max``
   Frequency window in THz for :math:`\omega/2\pi`.

``occupation_multiplier``
   Multiplicative factor applied to the initial Bose--Einstein occupation of
   eligible particles.

.. code-block:: toml

   [mode_excitation]
   enabled = true
   frequency_min = 2.0
   frequency_max = 2.5
   occupation_multiplier = 2.0

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

``[materials]``
---------------

This optional section enables an axis-aligned layered material stack for box
or STL geometry. If it is omitted, ``[io].material_folder`` keeps the original
single-material behavior.

``folders``
   Two or more material directories, ordered from the low-coordinate side to
   the high-coordinate side. Each directory must contain POSCAR and HDF5 data.

``interface_axis``
   ``"x"``, ``"y"``, or ``"z"``.

``interface_positions``
   Exactly ``folders.size()-1`` strictly increasing absolute positions in nm.
   Every interface must coincide with a grid-cell boundary and lie strictly
   inside the loaded geometry bounds. Each material layer must contain active
   cells. For STL geometry, these checks use the actual mesh bounds rather
   than ``geometry.sizes``; translated absolute coordinates are supported.

``interface_model``
   ``"diffuse_mismatch"`` (alias ``"dmm"``). This is currently the only model.

.. code-block:: toml

   [materials]
   folders = ["../material/Si/", "../material/Ge/"]
   interface_axis = "x"
   interface_positions = [50]
   interface_model = "diffuse_mismatch"
   interface_frequency_bin_thz = 0.1

``interface_frequency_bin_thz`` sets the shared DMM frequency-band width in THz
(default 0.1; finite and positive). Both materials use the same half-open bands.
It is a spectral discretization parameter, not interface roughness. Check
band-width convergence for the supplied q meshes. Multimaterial runs sample
equal physical state weights and use fixed incoming reservoir fluxes; the
material-specific carrier volumes are recorded in ``summary.txt``.

Thermal-reservoir facets must be normal to the interface axis. Periodic
boundaries normal to that axis are rejected because they would connect the two
end materials.

Strict validation summary
-------------------------

- use only literal TOML booleans ``true`` and ``false``
- do not quote numeric ``initial_temperature`` values
- use integer, positive ``grid_xyz`` dimensions
- keep all boundary arrays aligned
- use finite, non-negative thermal and roughness values
- assign every periodic facet to exactly one valid pair
- provide valid POSCAR and HDF5 material data
- align every layered-material interface with a grid-cell boundary

Finite-domain nonlinear Callaway and ribbon fields
--------------------------------------------------

``scattering.callaway_formulation = "nonlinear"`` is an opt-in for
``scattering.model = "callaway"``. It uses complete displaced Bose targets and
finite-step energy/momentum constraints. The default ``"linearized"`` keeps
historical fixed-reference behavior. In the nonlinear formulation,
``simulation.lifetime_temperature = "local"`` selects N/R rates from the old
cell temperature; a number freezes both channels at that temperature.
See :doc:`callaway` for the discrete model and its limits.

``simulation.conservative_boundaries = true`` lets RTA use the same elastic
frequency-shell reflections and fixed incoming reservoir flux as Callaway.
It does not linearize Bose occupations or the energy-temperature relation.

``geometry.extruded_z = true`` enables validated constant-height prism sampling,
projected-triangle ray intersection, and exact Cartesian cut-cell volumes for
polygonal ribbons specified as STL. Input STL coordinates remain in nm.

``simulation.write_cell_heat_flux = true`` writes ``grid_heat_flux.csv`` at the
convergence output interval for finite-domain simulations. Columns contain
cell index, time, temperature and all three heat-flux components. These fields
are written automatically for imposed-gradient runs. The usual convergence
file format is retained.

Acoustic/mixed interfaces
-------------------------

Layered materials may select ``interface_model = "mmm"`` with an explicit
``interface_amm_fraction`` in [0, 1]. Each event chooses the acoustic channel
with probability p and DMM with probability 1-p. Resistances are measured from
the temperature jump and actual flux; they are not interpolated. The p=0
endpoint retains the original DMM draw sequence. ``interface_model = "amm"``
selects p=1 directly.

For a nonzero acoustic fraction, provide one positive ``mass_densities_kg_m3``
value per material. The acoustic implementation uses geometric overlaps of
finite q cells, with dispersion linearized using the input group velocities.
``interface_parallel_bin_inv_a`` and, within AMM, ``interface_frequency_bin_thz``
only control the spatial search index. They do not discretize transmission
channels. Within DMM, the frequency parameter retains its original physical
band-sampling role. For example::

    [materials]
    folders = ["material/Si", "material/SiC"]
    interface_axis = "z"
    interface_positions = [10]
    interface_model = "mmm"
    interface_amm_fraction = 0.5
    mass_densities_kg_m3 = [2288.69038952758, 3173.18491635673]
    interface_frequency_bin_thz = 0.1
    interface_parallel_bin_inv_a = 0.05

These example densities were calculated from the relaxed input cells and are
not universal constants. A full uniform q mesh is required; symmetry-reduced
point sets must first be expanded. Acoustic reflection partners are found at
the mirrored reciprocal-grid address, resolving degenerate branches by their
frequencies and velocities. Reflection-only optical states may use a strict
frequency/velocity partner search when eigenvectors at the mirrored address
are not symmetry covariant. Unsupported banks fail explicitly.

Output identifies the new kernel as
``reciprocal_linearized_qcell_overlap_v2`` and records incident, transmitted,
no-overlap and capacity-correction flux weights. The old point-bin AMM kernel
artificially removed channels; its p>0 results should not be reused.

This remains a scalar dispersive acoustic approximation: the three lowest
branches transmit without polarization conversion, and optical modes reflect.
It does not solve full anisotropic elastic boundary conditions. A constant p
is a mixing probability, not a roughness height. Validate the phonon mesh,
finite-cell reconstruction, and non-reference equilibrium bias separately
from search-index independence and Monte Carlo sampling convergence.
