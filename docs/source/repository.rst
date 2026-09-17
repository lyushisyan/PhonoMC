Repository scope
================

The public tree contains the solver, documentation, reusable tools, selected
input examples and existing tutorial material/geometry files. Publication
runs and their results are maintained separately from the source release.

Included files
--------------

* ``include/`` and ``src/``: material, geometry, transport, collision, source,
  interface, accounting and output components.
* ``CMakeLists.txt``: standalone solver build, with optional local tests.
* ``docs/source/`` and ``README.md``: current methods, usage and limitations.
* ``example/``: explicitly selected TOML examples; no generated results.
* ``model/``: reusable tutorial geometry.
* ``material/Si`` and ``material/Ge``: existing tutorial tables;
  ``material/SiC``: the 3C-SiC NEP transport table with N/U and isotope rates.
  None of these folders contains a full collision matrix.
* Material conversion and plotting utilities in ``material/`` and ``tools/``.

Local files
-----------

``tests/``, ``calculation/``, ``analysis/``, ``output/``, ``manuscript/``,
generated example results, build directories and publication-specific material
exports are ignored by Git. Ignoring or removing a path from Git's index does
not remove its local files. Previously committed results remain in Git history;
this cleanup does not rewrite that history.

Tests are a separately retained developer resource. A fresh clone builds with
``BUILD_TESTING=OFF``; see :doc:`installation` for local validation. Scientific
figures, result archives and server campaign metadata are not needed to build
or run the solver. Scripts in the public tree must not require those folders.

Supported and experimental interfaces
-------------------------------------

TOML accepts ``rta`` and ``callaway``. The experimental full-matrix classes
remain in the C++ implementation, but neither ``full_matrix`` nor
``rta_linearized`` is an accepted public model name. Do not infer a supported
CLI workflow from a retained internal class. See :doc:`full_scattering`.

The examples use bundled material folders or an explicitly generated synthetic
Callaway demonstration. A material-conversion utility does not replace a
convergence study or establish the physical validity of supplied linewidths.
