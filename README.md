# EPMC (Standalone FBZ Workflow)

This repository is a standalone implementation of the EPMC simulation flow.
It contains executable code, material data, example inputs, batch scripts, and plotting tools.

## 1. Repository Layout

- `src/`, `include/`: C++ source and headers
- `Material/`: material packs (for example `Material/Si/`)
- `Model/`: geometry files (for example STL/OBJ)
- `input*.toml`: simulation input files
- `build/`: CMake build directory
- `results/`: simulation outputs
- `tools/plot_convergence.py`: convergence and steady-state plotting script
- `Material/expand_hdf5_to_fbz.py`: IBZ to FBZ HDF5 expansion script

## 2. Requirements

Core build/runtime:
- C++17 compiler
- CMake >= 3.20
- HDF5 C++ library (`find_package(HDF5 REQUIRED COMPONENTS CXX)`)
- OpenMP runtime (optional, auto-detected when available)

Optional (but recommended):
- `qdelaunay` (Qhull CLI) for volume tetrahedralization
  - If not available, code falls back to a star-tetrahedra method.

Python tooling (for FBZ expansion and plotting):
- `numpy`
- `h5py`
- `matplotlib`
- `phonopy` (needed by `expand_hdf5_to_fbz.py`)

If you use conda, run Python tools with your `phmc` env:

```bash
conda run -n phmc python ...
```

## 3. Build

```bash
cd /Users/sxliu/Documents/New\ project/epmc
./rebuild.sh
```

Equivalent manual commands:

```bash
cd /Users/sxliu/Documents/New\ project/epmc
cmake -S . -B build
cmake --build build -j
```

OpenMP is enabled by default when found. To disable explicitly:

```bash
cmake -S . -B build -DEPMC_ENABLE_OPENMP=OFF
cmake --build build -j
```

Optional: copy executable outside `build/`:

```bash
cp /Users/sxliu/Documents/New\ project/epmc/build/epmc \
   /Users/sxliu/Documents/New\ project/epmc/epmc
```

## 4. Run

Default run (auto-find config, prefer `input.toml` then `input.txt`):

```bash
./build/epmc
```

Set OpenMP thread count (example):

```bash
OMP_NUM_THREADS=8 ./build/epmc input_cross_100nm.toml
```

Run with explicit input:

```bash
./build/epmc input_cross_100nm.toml
./build/epmc input_inplane_fbz_x1000nm_yz100nm.toml
./build/epmc input_finfet_heat_source.toml
```

Notes:
- If no CLI argument is given, executable searches `input.toml` first, then `input.txt`.
- Relative paths in input (for example `--mat_folder`, `--results_folder`, model path) are resolved from the input file directory.
- `--results_folder` is auto-versioned as `<base>_<index>` (for example `results/Inplane_FBZ_x1000nm_yz100nm_1`).

## 5. Input File Format

Preferred format: TOML (`.toml`).

Legacy format (`.txt`) is still supported for backward compatibility.

### 5.1 TOML (recommended)

Use sectioned keys like:

```toml
[geometry]
model = "box"
dimensions = [100.0, 100.0, 100.0]

[simulation]
particles = 100000
timestep = 0.1
iterations = 6000
compute_thermal_conductivity = true
temp_dist = "cold"
subvolumes_mode = "slice"
subvolumes_count = 10
subvolumes_axis = 0

[boundary]
bound_pos_mode = "relative"
bound_pos_points = [[0.0, 0.5, 0.5], [1.0, 0.5, 0.5]]
bound_cond = ["T", "T", "P"]
bound_values = [302.0, 298.0]
connect_pos_mode = "relative"
connect_pos_points = [[0.5, 0.0, 0.5], [0.5, 1.0, 0.5], [0.5, 0.5, 0.0], [0.5, 0.5, 1.0]]

[io]
mat_folder = "Material/Si/"
results_folder = "results/Cross_10nm_300K/"
```

Also available in repo:
- `input_cross_100nm.toml`
- `input_inplane_fbz_x1000nm_yz100nm.toml`
- `input_finfet_heat_source.toml`

### 5.2 Legacy TXT (still supported)

The parser is option-token based:
- Lines beginning with `#` are comments.
- Options begin with `--`.
- Values continue until the next `--option`.

Supported options:
- `--model`: `box`, `cylinder`, or mesh file (`.stl` / `.obj`)
- `--dimensions`: geometry parameters
  - `box`: `Lx Ly Lz`
  - `cylinder`: `L R N`
- `--particles`: number of particles
- `--timestep`: time step
- `--iterations`: total timesteps
- `--compute_thermal_conductivity`: `true/false` (or `1/0`)
- `--subvolumes`:
  - `slice N axis`
  - `grid Nx Ny Nz`
- `--temp_dist`: `cold`, `hot`, `mean`, `random`, `linear`
- `--bound_pos`: `relative` or `absolute`, then 3D point triplets
- `--bound_cond`: per boundary point condition code
- `--bound_values`: values associated with non-periodic boundaries
- `--connect_pos`: periodic pair points (`relative|absolute` + even number of points)
- `--mat_folder`: material folder
- `--results_folder`: output base folder

Boundary condition codes:
- `T` or `F`: reservoir facet (absorbing/reinjecting through reservoir model)
- `P`: periodic facet (must be paired by `--connect_pos`)
- `R`: rough facet

Important for `--bound_values`:
- Values are consumed in order of non-`P` boundary points.
- For rough boundaries (`R`), value is used as roughness parameter `eta`.
- To avoid ambiguity, provide explicit values for all non-periodic boundary points in the same order as `--bound_pos`.

### 5.3 Local Volumetric Heat Source

You can enable a local body heat source in a box region.

TOML example:

```toml
[heat_source]
enabled = true
mode = "relative"                # relative | absolute
min = [0.45, 0.0, 0.0]           # region lower corner
max = [0.55, 1.0, 1.0]           # region upper corner
power_density = 1.0e15           # W/m^3
```

Legacy TXT keys:
- `--heat_source_enabled`
- `--heat_source_mode`
- `--heat_source_min x y z`
- `--heat_source_max x y z`
- `--heat_source_power_density value`

Notes:
- `mode="relative"` means coordinates are normalized by geometry bounds.
- `mode="absolute"` means direct model coordinates.
- Source is applied to subvolumes whose centers fall inside the region.
- `compute_thermal_conductivity` defaults to `false`; set to `true` when you want to output non-zero `kappa_fit` and `kappa_end`.
- If `heat_source.enabled = true` and `compute_thermal_conductivity = true` are both set, the program automatically disables thermal-conductivity estimation (heat-source runs focus on temperature distribution).

## 6. Material Data and FBZ Loading

Phonon loader (`src/Phonon.cpp`) resolves HDF5 in this order:
1. `kappa-fbz.hdf5`
2. `kappa.hdf5`
3. first `.hdf5` found in material folder

So if both exist, FBZ is preferred automatically.

Each material folder should contain:
- `POSCAR` (used to compute unit-cell volume and reciprocal lattice)
- one HDF5 file with required datasets (`frequency`, `group_velocity`, `gamma`, `qpoint`, `temperature`, etc.)

## 7. Predefined Runs in This Repo

Thermal conductivity (cross-plane):

```bash
./build/epmc input_cross_100nm.toml
```

Thermal conductivity (inplane FBZ, x=1000nm, y=z=100nm):

```bash
./build/epmc input_inplane_fbz_x1000nm_yz100nm.toml
```

Temperature distribution with local heat source (FinFET):

```bash
./build/epmc input_finfet_heat_source.toml
```

## 8. Output Files

Per result directory:
- `convergence.txt`
- `geometry_summary.txt`
- optional `plots/` (if plotting script is run)

`convergence.txt` columns are:
- `timestep`
- `T_sv_0 ... T_sv_N`
- `heatflux`
- `kappa_fit` (linear fit over middle subvolumes)
- `kappa_end` (endpoint temperature-difference gradient)

If `compute_thermal_conductivity = false`, `kappa_fit` and `kappa_end` are written as `0`.

The old mean-energy trailing column is intentionally removed.

## 9. Plotting Convergence and Steady Temperature

Example for two conductivity cases:

```bash
conda run -n phmc python tools/plot_convergence.py \
  results/Cross_100nm_300K_<idx> \
  results/Inplane_FBZ_x1000nm_yz100nm_<idx> \
  --tail 50 --error sem --summary-dir results/plots_summary
```

Per-case outputs in each `results/.../plots/`:
- `heatflux_convergence.png`
- `kappa_convergence.png`
- `temperature_convergence.png`
- `temperature_steady_tail50.png`
- `temperature_steady_tail50.csv`

Cross-case summary in `results/plots_summary/`:
- `compare_heatflux.png`
- `compare_kappa_end.png`
- `compare_kappa_fit.png`
- `compare_steady_temperature.png`

## 10. IBZ to FBZ HDF5 Expansion

Use this when you only have IBZ HDF5 and need explicit FBZ:

```bash
conda run -n phmc python Material/expand_hdf5_to_fbz.py \
  --hdf Material/Si_IBZ/kappa.hdf5 \
  --poscar Material/Si_IBZ/POSCAR \
  --out Material/Si/kappa-fbz.hdf5 \
  --plot --outdir Material/Si
```

The script:
- expands q-point stars with reciprocal symmetry,
- folds to FBZ using nearest-|k+G| rule,
- rotates group velocities consistently,
- writes FBZ datasets for C++ loader,
- optionally exports CSV/figures.

## 11. 3D Temperature Plot (Voxel Style)

For heat-source cases (for example FinFET), plot subvolume temperatures as 3D blocks:

```bash
conda run -n phmc python tools/plot_temperature_3d.py \
  --input input_finfet_heat_source.toml \
  --results results/FinFET_HeatSource_A_<idx> \
  --tail 50
```

## 12. Quick Smoke Test

```bash
./build/epmc input_cross_100nm.toml
```

This writes a run into `results/Cross_100nm_300K_*` for sanity check.

## 13. Generate FinFET STL

Use the built-in generator for a simplified FinFET geometry:

```bash
python3 tools/generate_finfet_stl.py --output Model/finfet.stl
```

Common parameters:
- `--base-width`, `--base-height`
- `--stem-width`, `--stem-height`
- `--thickness-y`
- `--scale`
- `--preview --preview-output Model/finfet_preview.png`

Then set your input model to:

```toml
[geometry]
model = "Model/finfet.stl"
```

Note: this generator builds a single simple FinFET outline prism (base + stem), with an unsplit flat top surface.
