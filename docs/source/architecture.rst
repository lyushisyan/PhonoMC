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
local validation.

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

``Transport`` and ``Interfaces`` handle geometry events. ``InterfaceModeSampler``
and ``AcousticInterfaceSampler`` construct interface channels;
``RoughBoundarySampler`` and ``ReservoirSampler`` construct boundary sampling
tables. Random-number generators are supplied by the solver.

``PrescribedHeatSource`` tracks requested, deposited and pending cell energy.
``PhononSourceOperator`` allocates that energy to eligible carrier occupations.
A cell without eligible carriers retains pending source energy rather than
creating carriers or depositing heat in another cell.

``BoundaryLedger`` and ``EnergyLedger`` record physical exchanges.
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

Some transport and boundary files remain private ``MonteCarloSolver`` member
implementations. Separate compilation units do not imply fully independent
ownership of all simulation state.
