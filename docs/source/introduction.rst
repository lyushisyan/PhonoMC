Introduction
============

PhonoMC models heat transport by propagating sampled phonon particles through
a three-dimensional domain. Particle states are initialized from material
mode data, advanced to surface collisions, scattered or transferred according
to boundary rules, and relaxed using temperature-dependent phonon lifetimes.

The program is intended for research workflows where geometry, boundary
scattering, and non-equilibrium phonon transport matter. It is not a generic
continuum heat-equation solver.

Simulation pipeline
-------------------

Each run follows this sequence:

1. Parse and validate a sectioned TOML input file.
2. Build a box or load an STL/OBJ surface mesh.
3. Assign and validate thermal, periodic, and rough boundaries.
4. Create the control-volume grid.
5. Load POSCAR and HDF5 phonon material data.
6. Initialize particles, modes, velocities, temperatures, and collision data.
7. Advance particles, inject sources, apply scattering, and update observables.
8. Write convergence, geometry, and runtime results.

Core components
---------------

``SimulationConfig``
   Strictly parses supported TOML keys and converts input values to typed C++
   configuration objects.

``SurfaceMesh``
   Stores vertices, faces, merged facets, adjacency, bounding-volume data,
   ray intersections, containment queries, and volume sampling support.

``SimulationDomain``
   Builds geometry, assigns boundary properties, pairs periodic facets, and
   constructs the simulation grid.

``PhononMaterial``
   Loads lattice volume and phonon-mode properties, then builds the
   temperature-energy lookup tables used by the solver.

``MonteCarloSolver``
   Owns particle state, performs transport and scattering, injects reservoirs
   and local heat sources, and writes time-dependent observables.

Units
-----

The external input and output conventions are:

- geometry ``sizes``: nm
- STL coordinates: nm
- roughness values: nm
- temperature: K
- time step and output time: ps
- heat-source power density: W/m³
- heat flux: W/m²
- thermal conductivity: W/(m K)

Internally, geometry is represented in angstroms.

Current scope
-------------

Geometry supports ``model = "box"`` or a path to an STL/OBJ mesh. The former
cylinder generator has been removed. Material loading is strict: both POSCAR
and the required HDF5 datasets must be valid before a production run starts.
