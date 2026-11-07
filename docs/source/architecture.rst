Solver architecture
===================

PhonoMC separates input and material loading, in-memory material physics,
particle evolution, conservation accounting and result formatting. The public
entry point reads a TOML configuration, constructs the domain and materials,
and advances ``MonteCarloSolver``.

Library targets
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Target
     - Responsibility
   * - ``phonomc_material_data``
     - Owned raw material arrays and collision-operator representations.
   * - ``phonomc_material_io``
     - POSCAR/HDF5 loading and schema validation; depends on HDF5.
   * - ``phonomc_material``
     - In-memory modes, linewidth interpolation, thermodynamics and Callaway data.
   * - ``phonomc_accounting``
     - Boundary and energy ledgers, independent of particle mutation.
   * - ``phonomc_output``
     - Serialization of output snapshots; no solver or HDF5 dependency.
   * - ``phonomc_core``
     - Configuration, geometry, material factory and transport orchestration.
   * - ``PhonoMC``
     - Command-line executable linked to the shared core.

Material construction
---------------------

``load_phonon_material(config, index)`` handles paths and strict file loading.
``PhononMaterial(data, temperature_lookup_dt)`` instead constructs an in-memory
model from owned ``MaterialData`` without filesystem or configuration access.
Production and explicit synthetic fixtures share the same validation and
mode-table initialization. Synthetic fallback is opt-in and intended only for
local validation. ``PhononMaterial.cpp`` constructs mode tables and active-mode
maps; ``material/Thermodynamics.cpp`` implements occupations, linewidth/lifetime
queries and temperature-energy lookup construction.
``PhononMaterialCollisionData.cpp`` supplies RTA and Callaway collision-data views.
These implementation units share the same owned model and preserve active-mode
ordering and reduction order.

``build_synthetic_test_modes()`` constructs the explicit synthetic fixture.
The public ``mode_frequency_window()`` returns a mesh-derived angular-frequency
tolerance in rad/ps, not ordinary frequency in THz or a scattering linewidth.
It is used when matching discrete modes at boundaries.

Configuration and geometry
--------------------------

``SimulationConfig.cpp`` is the public configuration facade. ``config/ConfigParsing.cpp``
handles assignment parsing, supported keys, aliases, defaults and unit conversion;
``config/ConfigValidation.cpp`` owns physical constraints and public simulation-control,
heat-source and scattering validators. Material and heat-source parsing are separate
functions; source lengths are converted together from nm to Angstrom. Numeric
conversion errors are reported at parsing, while physical constraints are checked
after conversion. Direct C++ solver construction validates timestep and population
before allocation; it does not clamp invalid controls to arbitrary defaults.
``config/OutputDirectory.cpp`` owns indexed result-directory creation. The small
configuration helper headers stay private under ``src/config``.

``geometry/MeshReader.cpp`` reads OBJ and ASCII/binary STL into owned raw vertex
and face arrays, preserving file order and STL vertex deduplication. It neither
scales coordinates nor samples particles. ``SimulationDomain`` resolves paths,
applies the existing STL nm-to-Angstrom conversion and origin shift, and assigns
boundary conditions and grid cells. ``geometry/DomainOutput.cpp`` serializes its
summary and grid centers; the private selector parser is shared with boundary
assignment so aliases and region counts remain consistent.

``SurfaceMesh`` owns topology, geometric caches and intersection queries. Edge IDs
and face incidence are built in one traversal, preserving first-encounter order.
Facet normals use the already oriented faces. The unused ``interfaces()`` accessor
and internal nonmanifold-interface cache have been removed; physical material
interfaces are represented by ``LayeredMaterialLayout``.
``geometry/SurfaceSampling.cpp`` implements its existing surface, rejection and
tetrahedral sampling methods. Distribution lifetimes, face/simplex ordering,
coordinate draws and fallback behavior are unchanged. No new RNG is introduced.

Particle and cell state
-----------------------

``ParticleStorage`` owns structure-of-arrays particle state. ``CellParticleIndex``
builds stable cell membership, reused across occupation-only updates.
``GridStatistics`` reconstructs energy, temperature and flux; ``BackgroundCache``
stores fixed-reference equilibrium values. Transport, absorption, insertion and
compaction invalidate membership as needed. Cell volumes and material-specific
carrier normalization enter all physical energy and flux sums consistently.

Time stepping
-------------

``MonteCarloSolver`` coordinates the following operations:

#. Apply the first half of any prescribed heat or gradient source.
#. Drift carriers and process surface, periodic, reservoir and material-interface events.
#. Reconstruct cell energy and temperature.
#. Apply the selected conservative collision update on existing carriers.
#. Apply the second half of the source and refresh observables.
#. Update energy ledgers and write results at the requested cadence.

``run_timestep()`` spells out this physical sequence. ``transport_particles()``
handles the existing/injected batches and compaction; ``Diagnostics.cpp`` owns
observation cadence and console formatting. Timers use the existing scope timer.

``refresh_thermodynamic_state()`` updates particle energies, cell energy densities
and temperatures, then particle temperatures. ``apply_phonon_collisions()``
dispatches the configured collision model. ``record_observables_if_due()``
updates heat flux and conductivity, samples convergence and appends output at
step 1 and at the configured ``convergence_write_interval``.

Solver resampling is unsupported: no validated configuration allowed the old
full-matrix/gradient-only integration. The standalone experimental modal resampler
remains available. Its unused solver RNG and counters have been removed, along with
the five resampling fields formerly written in ``[gradient_acceleration]``.
Warm-start and convergence fields remain.

Source half steps do not make the transport/collision composition second order.
See :doc:`methodology`, :doc:`callaway` and :doc:`temperature_gradient` for the
different collision targets and linear-response restrictions.

Implementation components
-------------------------

``RtaCollisionOperator`` implements finite-step Bose relaxation with an
energy-conserving auxiliary temperature. ``CarrierCollisionOperator`` implements
the linear-response RTA/Callaway update on the sampled carrier quadrature;
``NonlinearCallawayCollision`` handles the explicitly selected nonlinear path.
The full-matrix operator remains experimental (:doc:`full_scattering`).

``CallawayOperator::apply(energy)`` evaluates the linear collision generator
``L * energy`` (eV/ps); ``evolve(energy, dt_ps)`` evaluates the finite-time update
``exp(dt_ps * L) * energy`` (eV). Both use modal energy deviations from the fixed
reference equilibrium. ``CarrierCollisionOperator`` applies this evolution to
the existing equal-weight carrier samples, including repeated modes, rather than
requiring the full material mode grid in every cell.

``LayeredMaterialLayout`` owns the layer axis and interface planes. It locates
materials, computes the next plane crossing, and validates cell alignment and
boundary compatibility against ``SimulationDomain``. It has no particle state,
material physics or random-number generator. A point exactly on a plane belongs
to the upper layer; crossing times at or below the existing tolerance are skipped.

``interact_boundary`` evaluates absorption, periodic translation and wall
reflection from a ``BoundaryState``, ``BoundaryHit`` and narrow ``BoundaryPolicy``.
``InterfaceInteraction`` owns DMM/AMM/MMM sampling tables, borrows immutable
materials and applies the configured reference temperature and material energy
weights. Both return a value ``BoundaryOutcome`` without changing particle
storage, geometry, ledgers or global diagnostic counters. The outcome includes
the outgoing state, exchange energy, diagnostic counts and any mapping failure.
Failed mappings retain the incident modal state.

``TransportOperator`` applies those outcomes and owns the single-particle event
loop, position pushes, material/grid indices, ray-cache refresh, absorbed-particle
marking and transport diagnostics. It borrows immutable geometry, material layout,
material tables and samplers; each operation receives particle storage, the
boundary ledger and random streams explicitly. It has no solver back-reference,
source/collision operators, output writer or complete simulation configuration.

``advance_existing`` and ``advance_injected`` execute the two transport batches.
The latter advances each injected carrier only for the part of the timestep
remaining after its emission. The solver retains the original sequence: advance
existing carriers, compact, inject, build new caches, advance injected carriers,
compact again, then reduce the ledger. Structural storage changes occur outside
the operator's parallel loops, and the solver invalidates cell membership before
transport. Initialization and carrier-creating source/collision paths use the
same operator to build new boundary caches.

Physical flight times, numerical corrections, cache units and event ordering
are unchanged. Cached collision times still store remaining physical time divided
by the configured timestep. An exact surface/interface time tie selects the
surface. The event loop retains the 64-event guard and reports ray/mapping/guard
failures after leaving the parallel region. Reported transport or containment
failures latch the operator as unusable; the solver also retains its failed-step
guard. The operator exposes value diagnostic snapshots for the existing output
columns. Initialization of interface and rough-wall sampling tables remains in
``InterfaceInitialization.cpp`` and ``RoughBoundaryInitialization.cpp``.

``InterfaceModeSampler`` and ``AcousticInterfaceSampler`` construct interface
channels; ``RoughBoundarySampler`` and ``ReservoirSampler`` construct boundary
sampling tables. Random-number generators are passed explicitly into event
evaluation. DMM samples the reflected bank, then the transmitted bank, then
the transmission decision. MMM consumes a model-selection draw only for an
interior mixing fraction. The legacy direction fallback creates a fresh normal
distribution for each event, preserving its random-number consumption.
``TransportRandomStreams`` borrows the original serial and per-thread generators;
the two parallel batches retain ``schedule(dynamic, 64)``. An OpenMP batch checks
that every potential worker has a stream before it starts, preventing accidental
sharing of the serial fallback when the thread count exceeds the supplied pool.
The solver is explicitly non-movable because its operators borrow its services.

``include/solver/BoundaryFluxSampling.h`` contains internal residual diffuse-flux,
prefix-sum and weighted-window helpers used by ``RoughBoundarySampler`` and its
regression tests.

``PrescribedHeatSource`` tracks requested, deposited and pending cell energy.
``PhononSourceOperator`` allocates that energy to eligible carrier occupations.
A cell without eligible carriers retains pending source energy rather than
creating carriers or depositing heat in another cell.

``HeatSourceDriving.cpp`` contains solver-level source timing and deposition
orchestration. ``PhononSourceOperatorFullMatrix.cpp`` implements the experimental
full-matrix source path, including creation of carriers for unrepresented modes.

``BoundaryLedger`` and ``EnergyLedger`` record physical exchanges.
``BoundaryLedger`` also retains the previous step's reservoir absorption counts
for the legacy one-to-one injection policy. ``ReservoirInjector`` borrows the
configured ``ReservoirSampler``, geometry and material tables, and stores material
energy weights plus narrow injection options. It receives particle storage, the
ledger and one RNG explicitly. The solver selects fixed incoming flux for
conservative transport or multiple materials, and one-to-one otherwise. The
first one-to-one step uses rounded initial flux counts; subsequent steps use the
previous committed absorption counts, never the current step's counts.

Injection runs outside parallel regions, using the original serial RNG or OpenMP
worker-zero stream. Each reservoir draws its count (fixed-flux fractional rounding
only), all mode/time pairs, then all surface positions before processing any
particle. The injector preserves inward pushes, nonlinear/linearized occupations,
legacy velocity reversal and conservative outward-mode rejection. It appends each
particle, records its index/emission time, then records signed deviational energy
with that material's weight. Collision caches retain insertion defaults until the
solver calls ``TransportOperator``. It neither advances particles nor reduces or
commits the ledger. An injection exception aborts the step, without rollback or
retry. ``ReservoirInitialization.cpp`` only configures the sampler, initial counts
and injector; configuration and borrowed services stay fixed throughout a run.
``ResultSnapshots`` and ``ResultWriter`` format these values separately from the
solver. Output column names remain stable for existing analysis tools.

Conservation and failure contracts
----------------------------------

Collisions and heat sources preserve carrier count. In closed domains the
population remains fixed; open reservoirs can change it. Material interfaces
preserve represented deviational energy per event and do not branch carriers.
Mode discretization and approximate interface spectra still require convergence
and equilibrium checks beyond global energy conservation.

Unresolved transport events or invalid state abort the time step rather than
silently relocating or thermalizing a carrier. A failed solver instance must
not be resumed. Stable cell indexing and independent ledgers support the local
regression suite, which is optional and excluded from the public source tree.

Reservoir configuration, initialization and timestep orchestration still use
private ``MonteCarloSolver`` members. Material layout, single-event physics,
transport and reservoir injection can be tested without constructing a solver.
Local operator tests
exercise multiple reflections, periodic-cache refresh, interface/wall ordering,
injected flight times, compaction between batches, numerical ray corrections,
failure guards, per-thread RNG consumption and parallel absorption accounting.
Local injector tests check count policies across committed steps, exact RNG state
and batch ordering, hot/cold material-weighted energy, occupation policies, zero
emissions, uniform fallback, outward-mode behavior and insertion cache defaults.
The integration suite checks complete trajectories and energy ledgers. Identical
fixed-seed trajectories are compared using the same compiler/backend in serial
or one-thread OpenMP runs; different OpenMP thread counts are not expected to give
identical trajectories. The local MMM regression accepts an optional second
executable to compare pre/post-refactor output files for all seven scenarios.

Code organization and formatting
--------------------------------

``tools/`` contains reusable utilities such as ``format_cpp.py``,
``plot_convergence.py`` and ``export_phono3py_material.py``, alongside research
campaign scripts such as ``prepare_fet_comparison_campaign.py`` and
``prepare_gnr_size_scan.py``. Campaign scripts encode specific model directories
and may copy sources and companion scripts into frozen run directories. Their
existing locations and command-line workflows are retained.

``MonteCarloSolver::run_timestep`` retains the physical stage order; console
progress reporting lives in ``Diagnostics.cpp``. Material energy weights are
collected through one helper for interface, transport and reservoir setup.
``GridLocation.h`` provides the shared position-to-cell lookup, so initialization
and injection do not depend on the transport class for geometry queries.

The repository ``.clang-format`` defines four-space indentation, a 100-column
target, explicit control-flow braces and consistent include/function spacing.
All production C++ headers and sources under ``include/`` and ``src/`` can be
formatted and checked with::

    python3 tools/format_cpp.py
    python3 tools/format_cpp.py --check

The previous ``tools/format_solver.py`` command remains a compatibility entry
point for solver-only formatting, equivalent to ``format_cpp.py --scope solver``.

Use ``--clang-format /path/to/clang-format`` or ``CLANG_FORMAT`` when the executable
is not on ``PATH``. The current formatting was checked with clang-format 23.1.0;
use a consistent formatter version across contributors. The command discovers
new implementation subdirectories automatically; local tests, campaign tools
and generated data are outside its production-source scope.

Use descriptive names for policies and physical quantities, retain unit suffixes,
and prefer named option assignments when positional flags obscure intent.
Comments should explain physical constraints, numerical tolerances, ownership
and sampling order. Preserve arithmetic evaluation order, RNG calls and ledger
commit boundaries during readability changes; validate them with the existing
serial/OpenMP regressions and fixed-seed comparisons.

The local ``module_compatibility`` test checks raw mesh-reader units, complete
domain construction through OBJ and both STL encodings, quoted configuration
comments and validation precedence. Its optional ``--snapshot file`` output
records material mode/thermal queries, geometry, sampled points and full RNG
states for pre/post comparisons. Together with energy snapshots and the seven
MMM/DMM/AMM scenarios, these checks distinguish implementation organization from
changes in physical behavior or random-stream consumption.
