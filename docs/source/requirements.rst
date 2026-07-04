Requirements
============

Build requirements
------------------

- Linux or another Unix-like environment
- CMake 3.20 or newer
- a C++17 compiler, such as GCC 9+ or Clang 10+
- HDF5 with the C++ interface
- Qhull command-line tools, specifically ``qdelaunay``
- OpenMP for parallel execution (optional but recommended)

Ubuntu packages
---------------

.. code-block:: bash

   sudo apt update
   sudo apt install build-essential cmake libhdf5-dev qhull-bin libomp-dev

Python plotting requirements
----------------------------

The simulation executable does not require Python. The bundled plotting and
material-conversion tools use combinations of:

- Python 3.10+
- NumPy
- Matplotlib
- h5py
- phonopy (only for ``material/expand_hdf5_to_fbz.py``)

Documentation requirements
--------------------------

Install the documentation dependencies with:

.. code-block:: bash

   python3 -m pip install -r docs/requirements.txt

Material data requirements
--------------------------

Every material directory used by a simulation must contain:

- a valid ``POSCAR`` file
- ``kappa-fbz.hdf5`` or ``kappa.hdf5``

The HDF5 file must provide ``frequency``, ``qpoint``, ``mesh``,
``group_velocity``, ``temperature``, and ``gamma`` datasets with compatible
dimensions. Missing or invalid POSCAR/HDF5 data aborts the run.
