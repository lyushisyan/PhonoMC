Callaway dual relaxation
========================

Select ``scattering.model = "callaway"`` for a fixed-reference, linearized
Callaway model with separate normal (N) and resistive (R) relaxation. It retains
the material's actual mode frequencies, group velocities and Cartesian
wavevectors. This is an approximation to the collision integral, not an import
of phono3py's full LBTE collision matrix.

Configuration and material data
-------------------------------

.. code-block:: toml

   [simulation]
   background_temperature = 300
   initial_temperature = 301

   [scattering]
   model = "callaway"

All other geometry, time-step, particle and boundary settings remain required
as usual. ``matrix_files`` is forbidden. ``lifetime_temperature`` applies only
to RTA: Callaway fixes both N/R rates and equilibrium derivatives at
``background_temperature``. Positive-frequency stationary modes are included.

The material's ``kappa-fbz.hdf5`` must include ``gamma_N`` and ``gamma_U`` with
the same ``(temperature, qpoint, branch)`` shape and THz linewidth units as
``gamma``. Both arrays must be finite, nonnegative and satisfy
``gamma_N + gamma_U = gamma`` (relative tolerance 1e-7, absolute 1e-14 THz).
The N and U inverse lifetimes, before adding the optional resistive channels
below, are ``4*pi*gamma_N`` and ``4*pi*gamma_U`` in inverse ps. Rates are interpolated linearly at the reference
temperature; extrapolation is rejected. Zero rates remain zero.

Additional resistive scattering is enabled automatically when the HDF5 contains
``gamma_isotope``, ``gamma_impurity`` or ``gamma_defect``. The isotope name follows
phono3py; impurity and defect are PhonoMC extension field names. Each table must
be a finite, nonnegative THz HWHM linewidth with shape ``(qpoint, branch)`` or
``(temperature, qpoint, branch)``. Static tables are shared across temperatures;
temperature-dependent tables use the same temperature axis as ``gamma``.

The effective rates are

.. math::

   r_R = 4\pi(\gamma_U+\gamma_\mathrm{isotope}
                 +\gamma_\mathrm{impurity}+\gamma_\mathrm{defect}),
   \qquad r_\mathrm{RTA}=4\pi(\gamma+\gamma_\mathrm{isotope}
                 +\gamma_\mathrm{impurity}+\gamma_\mathrm{defect}).

Missing tables contribute zero. These fields must contain **independent**
contributions: ``gamma`` remains the three-phonon N+U linewidth and ``gamma_U``
remains U only. Do not supply a combined impurity table alongside its isotope or
defect components. No dataset is guessed from its name, and inverse-lifetime
``rate_*`` exports are not added a second time. Explicit boundary events are
not included in this sum. The full-matrix solver uses its supplied matrix;
this automatic rate addition applies to Callaway and RTA.

All resistive mechanisms share one energy-conserving ordinary equilibrium
projection constructed from the **total** R rate. This is a Callaway relaxation
approximation, not a frequency-resolved elastic disorder collision matrix.
The export tool preserves these fields when expanding to the full BZ and
recomputes total-rate RTA conductivity; original phono3py conductivity is kept
as ``kappa_phono3py``.

Generate N/U linewidths from the same phono3py run, for example::

   phono3py --mesh 11 11 11 --br --nu --ts 300

Use the force constants, primitive cell and calculation settings appropriate
to the material. Expand all scalar linewidth tables and velocities consistently
onto a complete Gamma-centred grid with unit q weights and first-Brillouin-zone
wavevector representatives. First-BZ membership is checked in the actual
reciprocal-lattice metric, not by wrapping each fractional component into a cube.
Duplicate and non-mesh qpoints are rejected. Match the POSCAR primitive cell to
the HDF5 data. See the `phono3py N/U documentation
<https://phonopy.github.io/phono3py/command-options.html#nu-n-u-true>`_ and
`linewidth units <https://phonopy.github.io/phono3py/input-output-files.html#gamma>`_.

The bundled Si/Ge files lack N/U tables. The bundled 3C-SiC NEP table includes
consistent ``gamma_N``/``gamma_U`` plus independent ``gamma_isotope``. Select a
reference temperature within its 100–1000 K range; the remaining mode-grid
and boundary compatibility checks still apply.
No automatic rescaling, inferred N/U ratio, or synthetic material fallback is
used to conceal missing/inconsistent Callaway input.

A runnable synthetic demonstration can be generated from the repository root::

   python3 tools/create_callaway_demo.py
   ./build/PhonoMC example/input_callaway.toml

This demo is an artificial mode bank for checking the executable workflow; its
results are not material predictions. Replace ``io.material_folder`` with
validated first-principles data for scientific use.

Conserving collision operator
-----------------------------

The model is

.. math::

   C_\lambda=-r_{N,\lambda}(n_\lambda-n^d_\lambda(T_N,\mathbf u))
             -r_{R,\lambda}(n_\lambda-n^0_\lambda(T_R)).

The displaced Bose distribution uses energy
:math:`\hbar\omega_\lambda-\hbar\mathbf q_\lambda\cdot\mathbf u`.
Both targets are linearized at :math:`T_0`. Their auxiliary temperatures and
drift satisfy the **rate-weighted** constraints

.. math::

   \sum_\lambda\hbar\omega_\lambda C_{N,\lambda}=0,\qquad
   \sum_\lambda\hbar\mathbf q_\lambda C_{N,\lambda}=0,\qquad
   \sum_\lambda\hbar\omega_\lambda C_{R,\lambda}=0.

For mode-dependent rates these auxiliary temperatures generally differ from
each other and from the energy-derived output temperature. Momentum means
crystal momentum; N scattering does not generally conserve heat current or
phonon number. Related formulations are discussed by
`Guo and Wang <https://doi.org/10.1103/PhysRevB.96.134312>`_.

The implementation eliminates the auxiliary variables analytically on the
**existing carriers in each cell**. In the formulas below, the index denotes
a carrier sample, not every entry of the material mode bank. Each carrier has
the same fixed statistical volume; duplicate modes remain separate samples.
Let
:math:`c_\lambda=\hbar\omega_\lambda\partial_Tn^0_\lambda`,
:math:`x_\lambda=f_\lambda/\sqrt{c_\lambda}`, and form the four columns
:math:`B_N=[\sqrt c,\sqrt c\,q_x/\omega,\sqrt c\,q_y/\omega,
\sqrt c\,q_z/\omega]`. For each mechanism,

.. math::

   S_j=R_j-R_jB_j(B_j^{\mathsf T}R_jB_j)^+B_j^{\mathsf T}R_j,
   \qquad B_R=\sqrt c,\qquad \dot x=-(S_N+S_R)x.

Here :math:`R_j` is the diagonal rate matrix and + denotes a pseudoinverse.
Rate-weighted, twice-orthogonalized basis columns implement this projection
without constructing a dense matrix. Rank-deficient directions, including
absent momentum components in 2D, are supported; relative residual columns
below 1e-12 are discarded. Each :math:`S_j` is symmetric positive semidefinite.
Both preserve energy; :math:`S_N` also annihilates the drift subspace.

The complete coupled exponential acts on individual carrier deviation energies. The code
does not freeze a target distribution and independently decay each mode, which
would fail finite-step conservation with unequal rates. Each matrix-vector
action and stored kernel requires O(N_cell) space/work up to a fixed small number
of constraint columns. A scaled Taylor action controls its truncation error
using the bound ``max(rate_N + rate_R)`` and a scaled norm at most four;
excessive work is rejected with a request to reduce the time step.

Transport, sources and boundaries
---------------------------------

``CarrierCollisionOperator`` constructs the rate-weighted conservation Gram
matrix from the actual carriers present in each cell. It updates occupations
only: count, mode, position, velocity and boundary cache remain unchanged.
Empty cells stay empty. No target modes are inserted. Restricting only the
state while retaining a full-grid conservation projector would be incorrect;
both are built on the same sampled quadrature here.

The finite carrier sample approximates the full mode integral. The collision
step preserves its discrete energy (and pure-N crystal momentum), but a sparse
sample can underresolve relaxation. A one-carrier cell, for example, cannot
relax while conserving that carrier's energy. Converge carrier count per cell,
spatial grid and timestep. Conservation alone does not establish transport
accuracy. This is a numerical-discretization change from earlier full-mode
carrier-completion runs; old and new trajectories must not be joined.

Initialization and thermal reservoirs use the same linearized equilibrium.
Reservoirs inject fixed incoming modal flux, so an open system's global count
can fluctuate from boundary injection/absorption even though collisions keep
it fixed. Temperature output uses ``T0 + delta_u/C(T0)``. Thermal sources weight
existing carriers by their reference-temperature heat capacities. Gaussian
spectral sources normalize their energy weights over eligible existing samples.
Empty cells and cells with no eligible samples retain their prescribed heat
in the local pending-energy ledger; sources do not create carriers.

The default ``rta`` retains its nonlinear Bose pseudo-temperature update on
existing carriers, with fixed or local lifetimes. In imposed-gradient runs,
``rta`` uses the linear-response energy projection on the same carrier sample
as Callaway. The public input has only ``rta`` and ``callaway``; no separate
``rta_linearized`` selection is needed. At small temperature perturbations the
ordinary RTA update approaches its linear-response limit; arbitrary finite
steps with unequal rates need not equal the coupled exponential.

Callaway reflecting boundaries use disjoint equal-frequency shells (relative
frequency tolerance 1e-8) and preserve equilibrium energy flux within each shell.
An imbalance exceeding relative 1e-7 plus a roundoff-sized absolute flux floor
is rejected during setup. The floor is 64 machine epsilons times the largest
shell flux and accommodates velocity noise on nominally stationary modes.
Discrete nearest
mirror-velocity matches can be many-to-one; their specular probabilities are
capped by each outgoing mode's equilibrium flux, with the remainder going into
the diffuse channel. This preserves shell equilibrium but is a boundary
discretization approximation: nominal zero roughness need not be perfectly
specular when the mode grid cannot represent the reflection. The roughness
factor uses the wavevector's normal component. Cross-frequency fallback is
disabled. Explicit geometric boundaries are not also added as lifetime rates.

The transport-then-collision ordering remains first-order Lie splitting.
An accurate collision exponential does not remove time-splitting errors or make
large time steps physically accurate. Strong N scattering can make the current
explicit transport scheme expensive. Existing material-interface scattering
still uses its frequency-window DMM approximation and needs separate interface
and mesh convergence checks.

Scope and verification
-----------------------

This version includes N relaxation and combined resistive relaxation from U
and any supplied isotope, impurity and defect linewidths. Electron scattering
and four-phonon mechanisms are not automatically included.
The reference temperature is fixed, and positive heat capacities must remain
numerically representable. Large temperature changes and nonlinear displaced
Bose dynamics are outside this version's scope.

Tests cover N energy/momentum conservation, stationary displaced equilibrium,
analytic resistive drift decay (including disorder with U=0), nonuniform N/R
rates, static and temperature-dependent extra HDF5 linewidths, zero rates, 2D rank deficiency,
dissipativity, exponential composition, HDF5 validation and units, full solver
energy accounting, empty-cell sources, and boundary equilibrium. They establish
the numerical model's contract; they do not validate a material's conductivity,
hydrodynamic window or second-sound speed against a full collision calculation.

Opt-in nonlinear finite-domain formulation
------------------------------------------

``scattering.callaway_formulation = "nonlinear"`` selects complete ordinary
and displaced Bose targets on the existing cell carriers. The default remains
``"linearized"`` for compatibility with earlier calculations. Nonlinear runs
use full Bose initialization/reservoirs and the nonlinear energy-temperature
relation. ``simulation.lifetime_temperature = "local"`` interpolates both N and
R rates at the old cell temperature; a numeric value freezes the rates there.
The fixed deviational background is retained and does not fix the local physics.
Out-of-range N/R temperatures are rejected, never extrapolated or clamped.

A step applies N(dt/2), R(dt), N(dt/2), with rates held fixed across that collision
step. For each channel, a_i=1-exp(-r_i dt_channel) and
n_i'=(1-a_i)n_i+a_i n_target,i. Auxiliary inverse temperature and drift parameters
are solved against **a-weighted** energy and, for N, all independent crystal
momentum constraints on the actual carrier sample. The target is
1/[exp(beta*hbar*omega-hbar*q.dot(beta*u))-1]. A damped Newton solve of the convex
Bose dual enforces a positive exponent and verifies the moment residual. The
update is a convex combination of nonnegative populations. It changes occupations
only and stages all cells before committing; an invalid moment solve aborts.
In the dt->0 limit these become the rate-weighted Callaway constraints. Finite
steps are not the exact nonlinear flow: channel splitting and target matching
require timestep convergence. Transport/collision composition remains first order.

These finite-domain runs use fixed incoming modal reservoir flux and elastic
frequency-shell reflection. RTA comparisons can opt into those same boundary
rules using ``simulation.conservative_boundaries = true``. The imported mode
bank, mesh, and reflection discretization still require convergence. A finite
sample conserving moments is not a validation against full LBTE.

For constant-height polygonal ribbons, ``geometry.extruded_z = true`` validates
vertical sides and flat top/bottom, samples the volume using top-face triangles,
and clips those triangles against Cartesian cells to obtain exact cell volumes.
It also uses projected-triangle ray intersection to handle extreme aspect ratios
and very small but nonzero group velocities. This option leaves the default
geometry workflow unchanged.
