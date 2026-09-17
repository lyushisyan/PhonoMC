Material Data
=============

The RTA datasets below do not contain mode-to-mode collision couplings. The
internal matrix-import work is experimental and not a public simulation
workflow; its restrictions are described in :doc:`full_scattering`. The bundled Si/Ge/SiC files
do not contain those matrices.

Bundled materials
-----------------

The public data directories are ``material/Si``, ``material/Ge`` and
``material/SiC``. All three production tables use 31×31×31 full q grids.
Si and Ge provide total three-phonon linewidths at 0–1000 K in 10 K steps,
without separate N/U tables. The 3C-SiC table spans 100–1000 K in 50 K steps
and includes N/U and independent isotope linewidths. Its force constants
come from a SiC-specific NEP potential, rather than direct DFT calculations.
See ``material/SiC/README.md``, ``provenance.json`` and ``validation.json``
for the source and export checks. Graphene is not bundled in this release.

Directory layout
----------------

A material directory should contain one POSCAR and one production HDF5 file:

.. code-block:: text

   material/Si/
   |-- POSCAR
   `-- kappa-fbz.hdf5

PhonoMC first looks for ``kappa-fbz.hdf5`` and then ``kappa.hdf5``. For
reproducibility, keep only the intended production HDF5 file in a material
directory. POSCAR loading is always strict. If neither preferred filename
exists, exactly one other ``.hdf5`` file may be used; ambiguous directories
are rejected.

POSCAR use
----------

PhonoMC reads the scale and three lattice vectors from POSCAR. Their
determinant gives the unit-cell volume in Å³, and the inverse lattice defines
the reciprocal basis used to convert fractional q-points to wave vectors.

An unreadable, incomplete, or singular POSCAR aborts material loading. The
program does not substitute an arbitrary cell volume.

Required HDF5 datasets
----------------------

``frequency``
   Shape ``(nq, nb)``. Frequencies in THz. Values are converted to angular
   frequency in rad/ps.

``qpoint``
   Shape ``(nq, 3)``. Fractional reciprocal coordinates.

``mesh``
   Shape ``(3)`` with positive integer dimensions. Their product must equal
   ``nq``: the current transport solver requires a full, equal-weight grid.

``group_velocity``
   Shape ``(nq, nb, 3)``. Group-velocity vectors used for propagation and
   active-mode selection.

``temperature``
   Shape ``(nT)``. Nonempty, finite, nonnegative, strictly increasing samples
   in K. A single temperature sample is supported.

``gamma``
   Shape ``(nT, nq, nb)``. Three-phonon scattering data converted to
   lifetime through :math:`\tau=1/(4\pi\gamma)`. All values must be finite and
   nonnegative, and the first dimension must match ``temperature``.

``gamma_isotope``, ``gamma_impurity``, ``gamma_defect`` (optional)
   Independent extra resistive linewidths, in the same THz HWHM units as
   ``gamma``. Shape ``(nq, nb)`` for static tables or ``(nT, nq, nb)`` for
   temperature-dependent tables. All values must be finite and nonnegative.
   RTA automatically uses ``gamma`` plus all present extras; Callaway adds
   them to its energy-conserving resistive channel together with ``gamma_U``.
   The impurity and defect names are PhonoMC extensions. These channels must
   not overlap each other or already be included in ``gamma`` / ``gamma_U``.
   See :doc:`callaway` for the rate sums and input conventions.

If ``weight`` is present, its shape must be ``(nq)`` and every value must be
one. Reduced-grid data must be expanded before loading; a filename alone does
not establish a full grid. These checks validate grid size and weights, not
the provenance or physical completeness of the supplied data.

Total linewidths (``gamma`` plus the optional extras) are linearly interpolated in temperature before computing
lifetimes. This replaces interpolation of lifetimes and can change results
between temperature samples. A zero linewidth gives an infinite lifetime
and exactly zero finite-step relaxation weight. The solver evaluates
:math:`1-\exp(-4\pi\gamma\Delta t)` directly, avoiding infinite-lifetime
arithmetic. Negative, non-finite, and missing linewidths are not treated as
instantaneous relaxation. Zero values used to mark *uncomputed* modes must
be resolved during data preparation; the loader cannot distinguish those
from physically zero scattering rates.

Dimension mismatches, missing datasets, or an empty active-mode bank abort the
run.

Active modes
------------

Only modes whose group-velocity norm exceeds the numerical threshold enter
the active mode list. PhonoMC also builds mappings between flattened
``(q,branch)`` indices and active-mode indices, and detects compatible
degenerate partner branches for boundary scattering.

In-memory material construction
-------------------------------

The application uses ``phonomc::load_phonon_material(config, index)`` from
``material/MaterialFactory.h`` for path resolution, loading diagnostics and the
explicit test-only synthetic fallback policy. The legacy
``PhononMaterial(config, index)`` constructor delegates to the same adapter.

Other data producers can instead supply an owned ``phonomc::MaterialData``
from ``material/MaterialData.h`` directly to
``PhononMaterial(data, temperature_lookup_dt)``. This validates the same raw
values and full-grid dimensions without using files, environment flags or
fallbacks. Data retains the HDF5 units and ordering described above. Direct
construction requires finite positive lookup spacing. The standalone
``phonomc_material`` target requires no HDF5 linkage; configuration-based
construction additionally requires ``phonomc_core``.

This interface does not import a full collision matrix or expand reduced grids.

Temperature range
-----------------

The Bose energy/temperature lookup is independent of the linewidth-temperature
samples. It covers 0 K through at least 1000 K (or the largest tabulated
temperature, if higher). This allows, for example, 49/51 K reservoirs with
scattering rates held at a tabulated reference temperature of 50 K. Energy
queries beyond the thermodynamic lookup bounds are clamped to those bounds.

RTA lifetime queries outside the HDF5 temperature samples retain the endpoint
linewidths. Choose tables that cover the requested lifetime temperatures.
Callaway instead rejects a reference temperature outside its N/U table; its
small local temperature perturbations may straddle an endpoint.

Full-zone expansion tool
------------------------

``material/expand_hdf5_to_fbz.py`` converts symmetry-reduced material data to
a full Brillouin-zone representation. It requires NumPy, h5py, and phonopy.

Inspect a source file and command-line options before processing production
data:

.. code-block:: bash

   python3 material/expand_hdf5_to_fbz.py --help

Synthetic material fallback
---------------------------

Setting ``PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1`` enables a test-only synthetic
mode bank if production material loading fails. Results from this fallback are
not suitable for physical thermal-conductivity claims. Do not set this
environment variable in production workflows.
