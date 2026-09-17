Methodology
===========

This chapter describes the equations implemented by the current code. It is a
code-level description, not a claim that every approximation is appropriate
for every material, geometry, or transport regime.

Phonon modes
------------

A mode is identified by a q-point and branch index. Modes with negligible
group-velocity magnitude are excluded from the active particle-mode bank.
For angular frequency :math:`\omega` and temperature :math:`T`, PhonoMC uses
the Bose--Einstein occupation

.. math::

   n(T,\omega) = \frac{1}{\exp(\hbar\omega/k_\mathrm{B}T)-1}.

The corresponding modal energy, including the zero-point term where required,
is

.. math::

   E_\mathrm{mode} = \hbar\omega\left(n+\frac{1}{2}\right).

The crystal energy density is normalized by the q-point count and POSCAR
unit-cell volume. A monotonic temperature-energy table is precomputed at the
``temperature_lookup_dt`` interval and is used in both directions:

- temperature to equilibrium energy density
- energy density to grid temperature

Particle representation
-----------------------

Each simulated particle stores mode, position, group velocity, occupation,
energy deviation, temperature, grid identifier, and cached information for
the next surface collision. The implementation uses structure-of-arrays
storage so the main OpenMP loops can process each state component efficiently.

Grid temperature update
-----------------------

Within each grid cell, PhonoMC compares particle occupation with an equilibrium
occupation at the selected background temperature. The deviation energy is

.. math::

   \Delta E_i = \hbar\omega_i\left(n_i-n_\mathrm{eq}\right).

Each carrier has the fixed spatial weight ``domain_volume / initial_count``.
The sampled deviation is converted to an energy density with that weight, the
active-mode count, the grid volume, and the crystal normalization, then added
to the exact equilibrium energy at ``background_temperature``. The temperature
is recovered from the inverse energy table.

The background reference is always a fixed numeric temperature. A moving local
reference is not supported because changing the reference as a carrier crosses
a grid boundary changes its represented deviational energy. The supported
combinations are:

``background_temperature = 0`` and ``lifetime_temperature = "local"``
   Full-phonon sampling with temperature-dependent local lifetimes.

``background_temperature = 300`` and ``lifetime_temperature = 300``
   Fixed-reference deviational sampling with linewidths held at 300 K.
   This choice does not linearize the Bose equilibrium target.

``background_temperature = 300`` and ``lifetime_temperature = "local"``
   Fixed-reference deviational sampling with local temperature-dependent
   lifetimes.

Lifetime scattering
-------------------

This section describes the default RTA path. Selecting
``scattering.model = "callaway"`` uses :doc:`callaway` for conserving N/U
dual relaxation on existing carrier quadrature. Both collision paths preserve
carrier count. The full-matrix adapter remains an internal experimental API,
documented in :doc:`full_scattering`, and is not a public TOML mode.

The HDF5 ``gamma`` value is converted to a lifetime in ps using

.. math::

   \tau = \frac{1}{4\pi\gamma}.

The total linewidth (``gamma`` plus independent isotope, impurity and defect
contributions when supplied) is interpolated before taking this reciprocal.
For ``gamma = 0``, the lifetime is infinite and the relaxation weight is
exactly zero. Finite-step weights are evaluated directly as
``-expm1(-4*pi*gamma*dt)`` rather than by dividing by an infinite lifetime.

Linewidths are linearly interpolated between tabulated temperatures and
clamped outside the available range for RTA. During one time step, each
occupation relaxes toward an auxiliary equilibrium target:

.. math::

   n(t+\Delta t) = n_0 + \left[n(t)-n_0\right]
   \exp\left(-\frac{\Delta t}{\tau}\right).

Here ``n_0 = n_BE(T_R)``. The cell auxiliary temperature ``T_R`` is chosen so
that the sum of modal energy changes, weighted by the finite-step relaxation
factors, vanishes on the actual carrier sample. It is distinct from the fixed
reference temperature used to represent deviations. ``lifetime_temperature``
selects rates evaluated at the reconstructed local temperature or at a fixed
numeric temperature. Periodic-gradient RTA instead uses the linear-response
operator described in :doc:`temperature_gradient`.

Boundary scattering
-------------------

Thermal boundaries absorb incident particles. Ordinary single-material RTA
uses one-to-one population refill. Callaway, conservative-boundary RTA and
multimaterial transport use prescribed incoming modal flux, so open-domain
particle count can fluctuate independently of collision count preservation.

Periodic boundaries translate a particle from one paired facet to the other.
The pairing is validated for opposite normals, area, and boundary-edge
spectrum.

For rough boundaries, the fallback specularity expression is

.. math::

   p_\mathrm{spec} = \exp\left[-\left(2\eta k|\cos\theta|\right)^2\right],

where :math:`\eta` is RMS roughness, :math:`k` is estimated from
:math:`\omega/|v_g|`, and :math:`\theta` is the incidence angle. The main
rough-boundary path additionally uses precomputed frequency-compatible mode
maps. Diffuse candidates inside the incoming mode's elastic frequency window
are sampled in proportion to their inward normal flux
:math:`|\mathbf{v}_g\!\cdot\!\mathbf{n}|`; this is required to avoid
over-populating nearly grazing modes.

Material interfaces
-------------------

For a layered box or STL domain, every carrier carries a material identifier in addition to
its mode. When a free flight reaches an internal plane, PhonoMC constructs
outgoing mode sets in shared, half-open frequency bands on both sides.
Each candidate is weighted by its normal phase-space flux
:math:`|v_{g,n}|/(N_q V_\mathrm{cell})`. The diffuse-mismatch transmission probability
for the common frequency band is

.. math::

   \mathcal{T}_{1\rightarrow2}(\omega) =
   \frac{\Phi_2(\omega)}{\Phi_1(\omega)+\Phi_2(\omega)}.

A transmitted or reflected outgoing mode is then sampled from its flux-weighted
bank. Multimaterial carrier positions are sampled in proportion to the active
state density :math:`M_m/(N_{q,m}V_{\mathrm{cell},m})`. Consequently all materials
use the same physical state weight per carrier. Their fixed statistical spatial
volumes differ; collisions, heat sources and reservoir injection use these same
material-specific volumes. Interface events do not split or create carriers.
Thermal reservoirs inject a prescribed incoming carrier flux independently of
absorbed carriers, including for nonlinear RTA.

The occupation deviation is rescaled by the frequency ratio so that each event
preserves represented deviational energy. A finite band still approximates
elastic frequency matching: shared bands restore reciprocity of the phase-space
transition channels but do not make different discrete frequencies identical.
Check band-width convergence and non-reference isothermal tests; exact energy
conservation or a zero-deviation test alone does not prove equilibrium preservation.
The band width, quadrature volumes, and interface event counts are written to
``summary.txt``. Single-material sampling remains unchanged.

The current implementation is deliberately restricted to axis-aligned layers
in box or STL geometry. Interfaces must coincide with grid-cell boundaries, so each
temperature cell has exactly one material energy-temperature relation.

Heat flux and conductivity
--------------------------

Grid heat flux is formed from the sampled sum of particle velocity times
deviation energy, followed by mode-count, crystal-volume, and unit
normalization.

When ``compute_kappa = true``, PhonoMC reports two Fourier-law estimates:

.. math::

   \kappa = -\frac{q}{\mathrm{d}T/\mathrm{d}x}.

``kappa_eff``
   Uses the two reservoir temperatures and full domain length.

``kappa_int``
   Uses a least-squares fit to grid temperature versus position.

The transport axis can be specified with ``simulation.transport_axis``;
otherwise it is inferred from the largest grid dimension. Conductivity
estimation is disabled when a local volumetric heat source is active.

With :doc:`temperature_gradient` driving, the periodic cell represents an
infinite-length film or bulk material. ``kappa_eff`` uses the imposed gradient,
not reservoir temperatures, and ``kappa_int`` is not applicable.

Prescribed lattice heat sources can use a power density or a total power
normalized over the sampled spatial profile. These are effective heating inputs,
not calculated electronic dissipation or electron-phonon emission rates.

The default modal allocation is a thermal occupation increment. Alternatively,
``spectrum = "weighted"`` distributes a cell's energy in proportion to each
eligible carrier's branch weight. Branch and frequency filters may restrict
either model. These weights describe energy allocation, not equal occupation
increments or exact integrated branch fractions.

The opt-in ``spectrum = "gaussian"`` assigns energy weights
``exp(-0.5*((f-f_center)/sigma_f)^2)`` to existing eligible carriers, with
normalization within each cell. RTA and Callaway share this sampled source
allocation. Repeated modes carry repeated quadrature weights; the realized
spectrum depends on finite mode sampling. It is a phenomenological model of
frequency-selective heating, not an electronic emission calculation. The
uniform-mode sampling limit includes the material density of states.
Callaway thermal injection instead uses reference heat-capacity weights.

Each cell retains requested, deposited and pending energy ledgers. Empty cells
and cells with no eligible samples retain their prescribed energy locally,
including after a heating pulse ends. Neither collision nor heating creates
carriers in RTA/Callaway. Only deposited energy contributes to the phonon
energy balance; a pending budget is a numerical backlog, not deposited heat.
It must be checked when converging spatial and spectral sampling.

For transient boundary-driven tests, an initial modal over-occupation can also
be applied through ``[mode_excitation]``. This multiplies the initial
occupation of modes in a selected frequency window while retaining uniform
full-band mode sampling, so the change represents a physical spectral
over-population rather than an importance-sampling change in the modal
quadrature.

Statistical interpretation
--------------------------

PhonoMC results contain Monte Carlo noise. A credible result should be checked
for convergence with respect to particles, time step, total simulated time,
grid resolution, and random sampling. Agreement between ``kappa_int`` and
``kappa_eff`` is useful evidence but is not, by itself, a convergence proof.

MMM interface operator
-----------------------

The event kernel is ``K_MMM = p K_acoustic + (1-p) K_DMM``. Existing carriers
retain their weighted deviational energy through an occupation remapping.
The following describes the scalar acoustic approximation, not a complete
elastic scattering matrix.

Each uniform reciprocal-grid point represents a finite parallelepiped.
Linearize ``omega(k)`` within it using the group velocity, then project the
cell into ``(omega/(2*pi), k_parallel_1, k_parallel_2)``. Intersections of the
left and right projected cells provide channels sharing frequency and parallel
wavevector. Convex-polyhedron clipping computes their geometry; positive
second-order tetrahedral quadrature integrates the scalar transmission factor
``4 Z1 Z2 cos(theta1) cos(theta2)/(Z1 cos(theta1)+Z2 cos(theta2))^2``.
The impedances use density times ``omega/|k|``, reconstructed at the common
quadrature point. Channels outside the overlap have no represented propagating
partner. Search bins only accelerate candidate lookup: sorted candidate pairs
produce the same transition table for all valid search-bin widths.

For incoming mode i, let ``w_i = |v_normal,i|/(Nq*Vcell)``. Divide this measure
by the projected cell volume to obtain its channel-space density. An overlap
between i and j receives the smaller of the two densities. In a consistent
uniform crystal grid these densities coincide. Neighboring linearized cells
can overlap on curved dispersion surfaces. A shared edge correction based on
both total overlap capacities bounds the two marginals before transmission
sampling. The resulting reciprocal edge weight ``g_ij`` gives probabilities
``g_ij/w_i`` and ``g_ij/w_j``; the remaining probability reflects. Consequently,
row probabilities remain bounded and the equilibrium carrier-flux measure is
stationary. Capacity corrections are reported explicitly.

Continuous overlap points share a frequency, but the stored carrier modes
remain centered on finite q cells. Therefore mapping back to those modes
introduces a finite-mesh spectral approximation. It preserves event-level
deviational energy and exactly preserves the reference equilibrium; it does
not guarantee exact Bose equilibrium at every temperature on a finite mesh.
Quantify this separately using the signed isothermal channel current for
``E_i(T)-E_i(T_ref)`` and mesh refinement. Conservation tests alone do not
establish physical AMM accuracy. Optical transmission, polarization conversion
and interface-specific force constants remain outside this approximation.
