# Material data

A simulation material folder needs `POSCAR` and a full-grid HDF5 phonon table. PhonoMC prefers `kappa-fbz.hdf5`, then `kappa.hdf5`, otherwise a single unambiguous `.hdf5` file. Required arrays, units, rate interpolation and optional disorder channels are specified in [the material documentation](../docs/source/materials.rst).

The public material set contains **Si, Ge and 3C-SiC**:

| Folder | Full q mesh | Temperature table | Scattering data |
|---|---|---|---|
| `Si/` | 31×31×31 | 0–1000 K, step 10 K | `gamma`; no separate N/U |
| `Ge/` | 31×31×31 | 0–1000 K, step 10 K | `gamma`; no separate N/U |
| `SiC/` | 31×31×31 | 100–1000 K, step 50 K | `gamma`, N/U and isotope linewidths |

Si and Ge retain their previously bundled tutorial data. [SiC](SiC/README.md) uses a SiC-specific NEP potential and matches the Si/SiC interface/device material export; its POSCAR and HDF5 hashes were verified against dongfang. Use recorded temperatures and meshes rather than assuming every dataset reproduces a particular manuscript figure. The Si/Ge tables require separately supplied consistent N/U data before Callaway can be selected. Graphene is excluded from the public source tree.

- `expand_hdf5_to_fbz.py`: full-zone preparation from symmetry-reduced data; inspect `--help` before use.
- `../tools/export_phono3py_material.py`: prepare transport material tables from phono3py inputs.
- `../tools/run_material_phono3py.py`: calculation/export utility; phono3py dependencies are optional and not needed by the C++ executable.
- `convert_phono3py_collision.py`: experimental matrix adapter, subject to the limitations in [full scattering](../docs/source/full_scattering.rst).

Publication-specific material calculations are ignored by Git. Do not combine independently tabulated channels if that double-counts scattering already present in `gamma`. For two-dimensional materials, the POSCAR cell height sets the volume convention used to report volumetric conductivity.
