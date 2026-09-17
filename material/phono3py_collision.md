# Importing phono3py full collision matrices

> **Theory audit, 2026-09-09:** this adapter is not valid for general transient
> distributions. The scalar reducible phono3py LBTE matrix uses a heat-current
> odd-subspace convention; full-grid expansion and `--project-energy-conservation`
> do not recover the missing even-subspace dynamics. The workflow below records
> the current implementation and is not a validated production workflow. See
> [the current scope](../docs/source/full_scattering.rst) and
> [Chaput, Eq. (2) and the discussion after Eq. (3)](https://arxiv.org/pdf/1303.4062).

`convert_phono3py_collision.py` imports the scalar collision operator on the
entire Brillouin-zone mesh. It requires Python 3, NumPy and h5py. The repository's
existing `kappa-fbz.hdf5` files contain linewidths, not full collision matrices;
the off-diagonal scattering data must be calculated with phono3py.

For a calculation with your existing force constants and phono3py configuration,
add these options, choosing the mesh and temperature appropriate to your system:

```sh
phono3py --mesh 5 5 5 --lbte --ts 300 --reducible-colmat --no-kappa-stars --write-collision
```

Retain the per-grid-point `collision-m555-g*.hdf5` files and the **same run's**
`kappa-m555.hdf5`. Calculate every grid point and every band; omit isotope,
electron-phonon and boundary mean-free-path scattering. The converter uses the
all-grid kappa file's q-points, frequencies and BZ grid-point identifiers to
identify each collision row and verifies its linewidths against the raw pieces.
Do not substitute an IBZ kappa file or a separately expanded IBZ file for this
collision metadata file. The PhonoMC material file may have a different q-point
ordering, provided its grid and frequencies match.

```sh
python3 material/convert_phono3py_collision.py \
  --collision collision-m555-g*.hdf5 \
  --collision-kappa kappa-m555.hdf5 \
  --material path/to/kappa-fbz.hdf5 \
  --temperature 300 \
  --input-format phono3py-raw \
  --out collision-300K.hdf5
```

The wildcard must select one integration method (one sigma, if applicable).
Temperature selection is exact; the converter does not interpolate matrices.
The source and material frequencies must agree to a relative tolerance of
`1e-10`. Material q-points may differ by reciprocal-lattice translations and by
up to `1e-6` in fractional coordinates, accommodating rounded FBZ coordinates.
Duplicate q-points modulo reciprocal-lattice vectors are rejected.

The converter also accepts one explicitly identified raw whole-grid tensor of
shape `(Ntemperature,Nq,Nband,Nq,Nband)` plus its separate `gamma` dataset. Use
the same `phono3py-raw` format. Per-grid-point exports are preferable because their
position in phono3py's computation is unambiguous across versions. Do not use
eigenvectors or a matrix retrieved after phono3py diagonalization as collisions.

For a matrix exported by your own phono3py API code **after adding the gamma
diagonal and completing any grid expansion, but before diagonalization**, use
`--input-format phono3py-assembled`. Supply `temperature`, `mesh`, `qpoint`,
`weight`, `frequency` and `collision_matrix` in that export, or provide matching
metadata with `--collision-kappa`. The tensor must use the same scalar shape.
The option explicitly asserts that its diagonal is already assembled; shape
alone cannot distinguish a raw matrix from an assembled one.

## Conventions and conservation

Let `Omega` be phono3py's scalar symmetric collision matrix in its native
linewidth units. Raw `--write-collision` data contain the coupling contribution,
including any self-coupling, while `gamma` supplies an additional diagonal:

```text
Omega = (Craw + Craw.T)/2 + diag(gamma)
D_i ∝ h*nu_i*sqrt(n_i*(n_i+1))
A_ij = -4*pi * (D_i/D_j) * Omega_ij
df/dt = A f
```

Here `nu` is the ordinary frequency in THz, `n` is the Bose distribution at the
selected reference temperature, and `f` is modal energy deviation. Thus the
linewidth-only diagonal decays at `4*pi*gamma` in inverse picoseconds. The
factor includes both the linewidth-to-lifetime factor of two and the
cycles-to-angular-frequency factor of `2*pi`. Signed off-diagonal entries are
retained; they must not be interpreted as probabilities or clipped to zero.

The converter verifies energy conservation column by column. A symmetric
`Omega` with `Omega*D=0` also makes `A*Cv=0`, since modal heat capacity is
proportional to `D**2`. If finite-mesh numerical integration violates the sum
rule, conversion fails by default and prints the residual. The optional
`--project-energy-conservation` explicitly replaces `Omega` with `P*Omega*P`,
where `P=I-u*u.T` and `u=D/||D||`. This restores the energy nullspace while
preserving symmetry and, for a positive-semidefinite input, positive
semidefiniteness. It changes the supplied operator: inspect the reported
relative Frobenius correction and mesh convergence before using projected
results. Projection is recorded in the output, never enabled implicitly, and
is not a substitute for correct source files.

## Supported scope and cost

- Regular, Gamma-centered meshes with all grid weights equal to one.
- Scalar full-BZ tensors; irreducible tensors with Cartesian axes are rejected.
- All branches; per-band fragments and incomplete per-grid collections are rejected.
- Pure phonon-phonon collisions. A nonzero `gamma_isotope` contains only a
  diagonal RTA loss and cannot supply the required conserving isotope kernel.
- Positive-frequency modes, including modes with zero group velocity.
  Nonpositive-frequency modes are omitted only when their rows and columns are
  numerically uncoupled; unstable coupled modes are rejected.
- One fixed positive reference temperature per output file.

This converter assembles dense work arrays before writing exact nonzero entries
to CSC. No sparsification threshold is applied. The default working-memory
estimate is limited to 2 GiB; `--max-memory-gb` changes that explicit limit.
A `31*31*31` mesh with six branches has about 179,000 modes, so even one dense
double matrix needs about 255 GB. Begin with a coarse mesh. Sparse storage does
not remove the intrinsic size of a dense full collision operator.

## Output schema

The generated HDF5 file uses numeric root datasets:

| Dataset | Shape | Meaning |
| --- | --- | --- |
| `schema_version` | scalar integer | `1` |
| `reference_temperature` | scalar double | Kelvin |
| `mode_indices` | `(M,2)` integer | Zero-based material q-point and branch |
| `frequency` | `(M,)` double | Positive mode frequencies in THz |
| `column_offsets` | `(M+1,)` int64 | CSC offsets |
| `row_indices` | `(nnz,)` int64 | CSC destination rows |
| `values` | `(nnz,)` double | Generator coefficients in `ps^-1` |

Modes are written in material q-major, branch-minor order. Each CSC column
identifies a source mode, and rows identify destination modes. Provenance,
representation and conservation diagnostics are stored as HDF5 attributes.
The converter refuses to overwrite existing output unless `--overwrite` is
specified and never overwrites an input file.

## Verification and source references

The following command requires the separately retained local `tests/` suite,
which is not part of the public source checkout.

```sh
python3 -m unittest discover -s tests/python -p 'test_phono3py_collision.py' -v
```

Tests use synthetic algebraic fixtures to verify diagonal assembly, the `4*pi`
conversion, signed couplings, Bose-to-energy transformation, q-point ordering,
per-grid collection, conservation and input rejection. They are not a physical
material benchmark.

The import conventions follow the upstream
[collision HDF5 writer](https://github.com/phonopy/phono3py/blob/develop/phono3py/file_IO.py),
[scalar collision assembly](https://github.com/phonopy/phono3py/blob/develop/phono3py/conductivity/collision_matrix_kernel.py),
[collision row calculation](https://github.com/phonopy/phono3py/blob/develop/phono3py/phonon3/collision_matrix.py)
and [CLI definitions](https://github.com/phonopy/phono3py/blob/develop/phono3py/cui/phono3py_argparse.py).
See the official [direct LBTE solution guide](https://phonopy.github.io/phono3py/direct-solution.html)
for phono3py calculation and memory requirements.
