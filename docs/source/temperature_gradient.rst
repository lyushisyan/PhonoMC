Infinite-length in-plane conductivity
=====================================

A uniform imposed temperature gradient allows a finite periodic simulation
cell to represent an infinitely repeated film. For a film of thickness ``d``
normal to z, choose x/y periodic boundaries, reflecting z surfaces and a drive
along x. The output is the fully developed in-plane conductivity
``kappa_parallel(d)`` without contacts at the ends of the heat-flow direction.
Making z periodic as well gives the homogeneous bulk problem.

Configuration
-------------

See ``example/input_inplane_gradient.toml`` for a complete Si example.
The essential settings are:

.. code-block:: toml

   [simulation]
   transport_axis = "x"
   grid_xyz = [1, 1, 20]
   compute_kappa = true
   initial_temperature = 300
   background_temperature = 300
   lifetime_temperature = 300

   [driving]
   temperature_gradient = 100000 # signed K/m along x

   [scattering]
   model = "callaway" # or "rta"

The x/y geometry sizes are periodic repeat lengths. The z size is the physical
film thickness. Change the z size to scan thickness, keeping the surface
roughness specified independently in nm. Changing the repeat length along x
must not change the converged conductivity of a homogeneous film.

The imposed gradient must be finite and nonzero, and ``transport_axis`` must
be explicitly ``"x"``, ``"y"`` or ``"z"``. The resolved faces normal to that
axis must be paired periodic boundaries. Thermal reservoirs are forbidden.
The current implementation supports a single-material box with RTA or
Callaway, a positive fixed reference temperature within the material data,
and uniform initial equilibrium at that reference. Local lifetime lookup,
volumetric heating, initial mode excitation and imported full-matrix collision
data cannot be combined with this drive. The latter restriction avoids implying
that an imported steady odd-parity matrix is a validated general transient
collision operator.

Periodic reference and source
-----------------------------

All coefficients are evaluated at ``T0``. Write the linearized local reference
as ``n0(T0) + (dn0/dT) G x`` and let ``g_lambda`` be the modal energy departure
from that affine reference. Only this departure is periodic; a physical
temperature field with a nonzero mean gradient is not strictly periodic.
The BTE for the departure is

.. math::

   \partial_t g_\lambda + \boldsymbol v_\lambda\cdot\nabla g_\lambda
   = \mathcal C_\lambda[g] - c_\lambda v_{\lambda a}G,
   \qquad G = \partial_a T.

Here ``c_lambda`` is the modal heat capacity, not heat capacity per volume.
The gradient source has positive and negative modal contributions and, for a
balanced full-BZ bank, zero total energy. It creates directional heat flow,
not uniform heating. The code rejects a relative imbalance in ``sum(c*v_a)``
larger than ``1e-8`` and removes only the accepted residual in the energy mode.
The implementation follows the affine-reference gradient-source formulation
of `Peraud and Hadjiconstantinou (2012)
<https://doi.org/10.1063/1.4757607>`_.

The source is evaluated on each cell's existing carriers. With equal fixed
statistical weights, define ``s_i = -c_i v_i,a G`` and apply

.. math::

   \Delta E_i^{\rm raw}=\Delta t\left(s_i-c_i
       \frac{\sum_{j\in cell}s_j}{\sum_{j\in cell}c_j}\right).

This projects out the finite sample's energy component, so the gradient
supplies no net heat in each nonempty cell. The physical carrier energy weight
remains ``M*Vp/(Nq*Vuc)``. Source and collision steps never create carriers;
empty cells stay empty. Finite mode sampling introduces quadrature error,
which must be reduced by carrier-count and spatial-grid convergence. The
full material bank is still required and checked for balanced velocities.

Collisions and boundaries
-------------------------

Gradient-driven RTA uses a conservative linearized modal operator at T0,
instead of the nonlinear carrierwise pseudo-temperature update used by legacy
thermal-reservoir RTA runs. With ``R_lambda = 1/tau_lambda`` its collision term is

.. math::

   \mathcal C^{\rm RTA}_\lambda[g]
   = -R_\lambda\left(g_\lambda-c_\lambda
     \frac{\sum_\mu R_\mu g_\mu}{\sum_\mu R_\mu c_\mu}\right).

Zero rates give zero collision terms. This operator preserves energy and the
uniform temperature mode. Both its state and conservation projector are
constructed on the existing carrier sample and advanced with the same coupled
exponential used by Callaway. Positive-frequency stationary modes remain
eligible for initial sampling. RTA uses total ``gamma`` plus independently
supplied isotope/impurity/defect linewidths. Callaway keeps its N energy/momentum
projection and combined resistive energy projection unchanged.

Each carrier retains its own diagonal attenuation, including opposite signed
weights at different positions. The conservation gain is constructed with
sample multiplicities; it is not divided by the number of carriers of a mode.
Ordinary thermal-reservoir RTA retains its nonlinear Bose update.

Both models use signed-weight, physical-mode reflections and the same
equal-frequency, equilibrium-flux-balanced surface kernel at T0. A roughness
value of zero requests specular reflection, but exact specularity requires the
discrete mode bank to represent the mirror modes. Large roughness does not
mean every mode is diffuse: the model uses each mode's normal wavevector.
Surface scattering is not also added to the bulk linewidth.
Shell flux-balance validation includes a floating-point absolute floor
(``64*epsilon`` times the largest resolved shell flux), so velocity noise on
nominally stationary modes does not cause a spurious relative-error failure.
The relative tolerance remains ``1e-7``; resolvable shell asymmetry is rejected.

Output and convergence
----------------------

``heatflux`` is the volume-averaged heat flux along the explicit axis, even
when most grid cells lie along a different axis. With the drive enabled,

.. math::

   \kappa_{\rm eff}=-\langle q_a\rangle/G.

``kappa_eff`` therefore means the imposed-gradient conductivity, without contact
jumps. ``kappa_int`` is zero and not applicable. The summary records this
definition. ``T_1``, etc. represent ``T0 + periodic temperature correction``;
they do not include the imposed affine temperature ramp. The thermal energy
ledger likewise tracks the periodic departure plus the uniform T0 background.
Gradient occupation updates and actual net energy increments are recorded in
separate columns appended to ``convergence.txt``; ordinary heat-source columns
are unaffected. Legacy runs retain their existing column layout.

Gradient runs also write ``grid_heat_flux.csv`` at the convergence-output cadence.
Each row contains timestep, time in ps, zero-based cell index, periodic-reference
temperature in K, and signed x/y/z heat flux in W/m². Join the index with
``grid_centers.csv`` for positions and cell volumes. For a 1×1×Nz film this gives
the layer-averaged ``qx(z)`` profile; average over a stated time window after
checking convergence. These samples are not automatically steady-state results.
Older executables do not write this file; spatial flux cannot be reconstructed
from their volume-averaged ``convergence.txt`` alone.

The source is applied in two half steps around transport/collision. The
transport/collision split remains first order; source half steps do not make
the whole solver second order. Converge the time step, thickness grid,
statistical population and averaging duration, and verify the gradient sign
and amplitude do not change conductivity. A finite q mesh and the Callaway
approximation retain their physical limitations. Pure momentum-conserving
systems without any momentum relaxation need not have a finite steady bulk
conductivity; the drive does not impose artificial damping.

Tests compare paired-stream bulk transients to analytic RTA/Callaway solutions,
check SI units, signed source normalization and empty-cell count preservation, periodic repeat
length and gradient sign/amplitude invariance, independent disorder rates,
energy accounting, and a discrete-stream diffuse-wall thin-film RTA solution.
Specular films recover the corresponding bulk limit on the symmetric test bank.

Optional acceleration and stationarity stopping
-----------------------------------------------

The following settings are opt-in. Omitting them preserves the original
initialization, carrier population policy and fixed iteration count::

   [sampling]
   resample_interval = 0  # required for existing-carrier RTA/Callaway

   [driving]
   temperature_gradient = 100000
   warm_start_ps = 0      # 0: equilibrium; positive: bulk RTA response initial guess

   [convergence]
   enabled = true
   min_time_ps = 2000
   window_ps = 500
   windows = 4
   relative_tolerance = 0.005
   absolute_tolerance = 0.01  # W/(m K)

Nonzero ``resample_interval`` is rejected for the existing-carrier collision
scheme. The former modal merging operation preserves signed energy but changes
its quadrature weights, hence its relaxation target. Retaining it would require
explicit statistical-weight tracking, beyond a simple population cap. No
near-zero carriers are deleted. In closed periodic/reflecting domains the
carrier count remains exactly the initialized count.

A positive ``warm_start_ps = a`` initializes the modal response
``g_lambda = -c_lambda*v_lambda,a*G*(1-exp(-r_lambda*a))/r_lambda``,
using the total RTA rate including supplied disorder. Zero-rate modes use the
continuous limit ``a``. The uniform temperature component is removed and the
initial energy ledger is initialized afterwards. This is a bulk initial guess
for both RTA and Callaway, not a steady film solution, not a restart checkpoint,
and not an elapsed-time offset. The actual boundary and collision operators
remain unchanged. In films it can overshoot or converge from above; compare
equilibrium and warm starts before claiming a reduction in required time.
Warm-started curves must not be interpreted as physical switch-on transients.

Stopping discards samples before ``min_time_ps``, integrates the sampled global
and cell ``-qx/G`` values by piecewise-linear quadrature over non-overlapping
windows, and requires the range of **every** component across the latest
``windows`` complete windows to be no greater than
``absolute_tolerance + relative_tolerance*abs(mean_global_kappa)``.
The mean global conductivity must exceed the absolute tolerance. At least four
windows and four output intervals per window are required. ``iterations`` is
still a hard maximum. A passed screen is empirical stationarity, not a bound
on the remaining long-time tail, an independent-sample confidence interval,
or a substitute for timestep, spatial, modal and population convergence.

``summary.txt`` records whether this screen passed, the final time, the checked
window means (global followed by cells), thresholds, warm-start age, removed
carrier count, resampling wall time and its energy residual. An early stop is a
normal successful process exit; consumers must inspect the recorded final time
and termination reason rather than assume ``iterations*time_step`` was reached.
