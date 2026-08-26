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
   Fixed-reference deviational sampling with lifetimes linearized at 300 K.

``background_temperature = 300`` and ``lifetime_temperature = "local"``
   Fixed-reference deviational sampling with local temperature-dependent
   lifetimes.

Lifetime scattering
-------------------

The HDF5 ``gamma`` value is converted to a lifetime in ps using

.. math::

   \tau = \frac{1}{4\pi\gamma}.

Lifetimes are linearly interpolated between the tabulated temperatures and
clamped outside the available range. During one time step, occupation relaxes
toward equilibrium:

.. math::

   n(t+\Delta t) = n_0 + \left[n(t)-n_0\right]
   \exp\left(-\frac{\Delta t}{\tau}\right).

``lifetime_temperature`` independently selects a local particle temperature
with ``"local"`` or a fixed numeric temperature.

Boundary scattering
-------------------

Thermal boundaries absorb incident particles. Reservoir injection uses the
previous-step leaving count for one-to-one population refill.

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

The transport axis is inferred from the largest grid dimension. Conductivity
estimation is disabled when a local volumetric heat source is active.

Statistical interpretation
--------------------------

PhonoMC results contain Monte Carlo noise. A credible result should be checked
for convergence with respect to particles, time step, total simulated time,
grid resolution, and random sampling. Agreement between ``kappa_int`` and
``kappa_eff`` is useful evidence but is not, by itself, a convergence proof.
