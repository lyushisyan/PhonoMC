Results and Analysis
====================

Result-directory naming
-----------------------

``output_folder`` is a base name. PhonoMC scans sibling directories and creates
the next numeric suffix, for example ``Cross_100nm_0`` followed by
``Cross_100nm_1``. Record the exact printed result path in batch workflows.

``convergence.txt``
-------------------

The header begins with ``#`` and defines all columns. The current columns are:

``timestep``
   Completed solver step.

``time_ps``
   Simulated elapsed time in ps.

``T_1 ... T_n``
   Grid temperatures in K, in the same index order as ``grid_centers.csv``.

``heatflux``
   Volume-weighted average heat flux along the inferred transport axis in
   W/m².

``kappa_int``
   Conductivity from a least-squares grid-temperature gradient.

``kappa_eff``
   Conductivity from reservoir endpoint temperatures and full domain length.

``particle_count``, ``absorbed``, ``injected``, ``recovered``, ``net``
   Per-step reservoir and recovery bookkeeping. ``net`` is reservoir injection
   minus absorption. A local heat source does not change ``particle_count``.

``reservoir_absorbed_energy_ev`` / ``reservoir_injected_energy_ev``
   Deviational energy removed and supplied by thermal reservoirs in the step.

``hs_occupation_updates`` / ``hs_injected_energy_ev``
   Number of carrier-occupation updates made across the two split source
   half-steps and the prescribed physical source energy added in the step.

``lifetime_energy_residual_ev``
   Numerical energy residual of the finite-step lifetime scattering solve. It
   should remain close to floating-point roundoff.

``energy_balance_residual_ev``
   Remaining step-wise closure error after accounting for the heat source,
   reservoir exchange, and the recorded lifetime-scattering residual. A
   closed or correctly accounted domain should keep this near roundoff.

``total_thermal_energy_ev``
   Domain-integrated thermal energy with the zero-point contribution removed.

Temperature convergence
-----------------------

Look for a statistically stationary profile rather than a perfectly flat
trace. Increase simulation time if the profile is still drifting. Increase
particle count if the steady trace is dominated by noise.

For a nominally one-dimensional case, plot several representative ``T_i``
traces and inspect the final spatial profile using ``grid_centers.csv``.

Particle balance
----------------

Reservoir-driven steady cases should not exhibit a persistent unexplained
particle drift. Compare cumulative absorption and injection in ``summary.txt``
and inspect ``recovered`` events. Frequent recovery indicates a geometry or
collision-tolerance problem, not a normal convergence mechanism.

Conductivity interpretation
---------------------------

Compare ``kappa_int`` and ``kappa_eff`` after the temperature profile has
settled. Differences can arise from boundary jumps, nonlinearity, insufficient
grid resolution, or sampling noise.

Repeat a conductivity calculation with:

- more particles
- a smaller time step
- longer simulated time
- at least one finer grid

Report the convergence study alongside the selected conductivity value.

``summary.txt`` diagnostics
---------------------------

The summary includes normalized input, geometry, grid and boundary counts,
runtime, OpenMP settings, rough-scattering selection statistics, escaped
particle recovery, reservoir totals, and heat-source injection totals.

Large rough-boundary fallback counts deserve investigation. The rough-boundary
summary also reports normal-flux-weighted frequency-window selections and the
number of guarded uniform fallbacks; the latter should normally remain zero.
Nonzero or large fallback counts can indicate
that no frequency-compatible reflected mode was available for many events.

One-dimensional plots
---------------------

.. code-block:: bash

   python3 tools/plot_convergence.py RESULT_DIRECTORY \
     --max-temp-lines 40 --show-legend

Outputs:

- ``plots_1d/temperature_vs_time.png``
- ``plots_1d/heatflux_vs_time.png``
- ``plots_1d/kappa_vs_time.png``

Three-dimensional plots
-----------------------

.. code-block:: bash

   python3 tools/plot_temperature_3d.py \
     --input INPUT.toml \
     --results RESULT_DIRECTORY \
     --x-slice-rel 0.5 \
     --y-slice-rel 0.5

Use slice images to check whether the heat-source region, periodic directions,
and thermal contacts produce the intended spatial pattern.
