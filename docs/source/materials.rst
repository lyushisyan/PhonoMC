Material Data
=============

Directory layout
----------------

A material directory should contain one POSCAR and one production HDF5 file:

.. code-block:: text

   material/Si/
   |-- POSCAR
   `-- kappa-fbz.hdf5

PhonoMC first looks for ``kappa-fbz.hdf5`` and then ``kappa.hdf5``. For
reproducibility, keep only the intended production HDF5 file in a material
directory. POSCAR loading is always strict.

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
   One-dimensional dataset containing at least three mesh dimensions.

``group_velocity``
   Shape ``(nq, nb, 3)``. Group-velocity vectors used for propagation and
   active-mode selection.

``temperature``
   Shape ``(nT)``. Temperature samples in K.

``gamma``
   Shape ``(nT, nq, nb)``. Temperature-dependent scattering data converted to
   lifetime through :math:`\tau=1/(4\pi\gamma)`.

Dimension mismatches, missing datasets, or an empty active-mode bank abort the
run.

Active modes
------------

Only modes whose group-velocity norm exceeds the numerical threshold enter
the active mode list. PhonoMC also builds mappings between flattened
``(q,branch)`` indices and active-mode indices, and detects compatible
degenerate partner branches for boundary scattering.

Temperature range
-----------------

Energy and lifetime queries outside the HDF5 temperature samples are clamped
to the end samples. Choose material tables that cover the expected simulation
temperature range; otherwise a run may remain numerically valid while using
endpoint material properties over a large interval.

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
