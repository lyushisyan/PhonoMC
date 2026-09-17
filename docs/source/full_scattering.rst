Experimental full scattering matrix
===================================

Public TOML input accepts only ``rta`` and ``callaway``. The parser rejects
``full_matrix`` and ``rta_linearized``. The C++ matrix representations,
reader and collision kernels remain available for internal development;
they are not a supported route for production phono3py transient calculations.

Current limitation
------------------

The scalar reducible phono3py LBTE matrix used by the adapter follows a
heat-current odd-subspace convention. Expanding its q grid and projecting an
energy null space does not reconstruct the missing even-subspace dynamics.
Consequently the current adapter is not a validated general transient
collision operator. See the operator convention discussed in
`Chaput (2013) <https://arxiv.org/abs/1303.4062>`_.

Frequencies, group velocities and linewidths in ordinary kappa files do not
contain mode-to-mode collision couplings. Neither renaming a kappa file nor
expanding an irreducible grid supplies those couplings.

Retained implementation
-----------------------

``material/ScatteringMatrix`` owns the matrix representation and validates
its dimensions and conservation properties. ``ScatteringMatrixReader`` reads
the versioned HDF5 representation. ``FullMatrixCollisionOperator`` and
``FullMatrixSource`` contain internal update kernels; their presence in the
build does not enable a CLI model.

``material/convert_phono3py_collision.py`` and
``material/phono3py_collision.md`` retain the adapter and its historical file
format documentation. Numerical conservation checks, unit conversion and
matrix symmetry are necessary but do not resolve the operator-convention
limitation above. These files should be used for development and inspection,
not as a claim of validated full-LBTE device transport.

For production simulations, use :doc:`methodology` for RTA, or the stated
approximation and data requirements in :doc:`callaway`.
