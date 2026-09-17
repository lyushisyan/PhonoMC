# PhonoMC

**First-principles-based deviational phonon Monte Carlo with temperature-dependent scattering.**

[English](#overview) · [中文](#中文说明) · [Русский](#русский)

## Overview

PhonoMC is a C++17 solver for phonon heat transport in thin films and device geometries. It imports mode-resolved frequencies, group velocities and scattering linewidths from full-Brillouin-zone HDF5 material tables. Computational carriers store occupations relative to a fixed equilibrium reference; local temperatures are reconstructed from their energy deviations.

The primary workflow uses energy-conserving RTA occupation updates with fixed or locally evaluated scattering rates. The reference temperature and the temperature used to look up rates are separate settings. Collisions and heat-source updates act on existing carriers; thermal reservoirs can change the population through absorption and injection.

## Capabilities

- **Transport:** cross-plane films, finite-length in-plane films, and infinite-length in-plane conductivity through periodic-gradient driving.
- **Collisions:** public `rta` and `callaway` selections. Callaway requires separate N/U rates; supplied isotope, impurity and defect linewidths contribute to resistive scattering. Gradient-driven transport uses linear response at a fixed temperature.
- **Geometry:** boxes and STL/OBJ surface meshes, rough reflecting surfaces, periodic boundaries and thermal reservoirs.
- **Interfaces:** layered materials with DMM, an acoustic-channel AMM approximation, or MMM mixing. The MMM parameter is the AMM mixing fraction, not an RMS roughness height.
- **Heating:** prescribed box or Gaussian spatial sources, power-density or total-power normalization, and thermal or weighted frequency/branch excitation.
- **Outputs:** cell temperatures, heat flux, conductivity estimates, source and boundary energy accounting, and OpenMP profiling.

The full-matrix classes remain experimental internal code and are **not** a public input mode. See [scope and limitations](docs/source/full_scattering.rst). Callaway and interface approximations require their own physical and discretization validation.

## Build

Requirements: CMake 3.20+, a C++17 compiler, and HDF5 with its C++ interface. OpenMP is optional. General surface-mesh workflows also use the Qhull `qdelaunay` executable; box geometry does not require it.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPHONOMC_ENABLE_OPENMP=ON
cmake --build build -j 4
OMP_NUM_THREADS=4 ./build/PhonoMC example/input_cross_100nm.toml
```

Check CMake's output to confirm OpenMP was found. Linux/macOS setup and the optional local regression suite are described in [installation](docs/source/installation.rst). Public builds default to `BUILD_TESTING=OFF` and do not require `tests/`.

## Choose an example

| Task | Input |
|---|---|
| Cross-plane Si film | [input_cross_100nm.toml](example/input_cross_100nm.toml) |
| Finite in-plane Si film | [input_inplane_x1000nm_z100nm_r1nm.toml](example/input_inplane_x1000nm_z100nm_r1nm.toml) |
| Infinite in-plane Si film, RTA | [input_inplane_gradient.toml](example/input_inplane_gradient.toml) |
| Gradient drive with stationarity screening | [input_inplane_accelerated.toml](example/input_inplane_accelerated.toml) |
| FinFET with localized heating | [input_finfet_stl_heat1e20.toml](example/input_finfet_stl_heat1e20.toml) |
| DMM Si/Ge bilayer | [input_bilayer_si_ge.toml](example/input_bilayer_si_ge.toml) |
| MMM Si/Ge bilayer | [input_bilayer_mmm.toml](example/input_bilayer_mmm.toml) |
| Prescribed total power | [input_box_total_power.toml](example/input_box_total_power.toml) |
| Gaussian frequency weighting | [input_gaussian_spectrum.toml](example/input_gaussian_spectrum.toml) |
| Synthetic Callaway demonstration | [input_callaway.toml](example/input_callaway.toml) |

Examples are starting configurations, not converged publication results. The Callaway demonstration requires `python3 tools/create_callaway_demo.py`; real Callaway calculations require validated `gamma_N` and `gamma_U` tables. Do not interpret synthetic material outputs as physical predictions.

## Material and input conventions

The bundled materials are **Si, Ge and 3C-SiC**. Si/Ge retain the existing tutorial tables; SiC uses a SiC-specific NEP potential with N/U and isotope linewidths (100–1000 K). Each material directory contains `POSCAR` and a compatible HDF5 table. File precedence is `kappa-fbz.hdf5`, then `kappa.hdf5`, then a single unambiguous `.hdf5` file. The loader requires a full, equal-weight q grid; it does not automatically expand an irreducible kappa file. See [materials](docs/source/materials.rst) and [material/README.md](material/README.md).

Use sectioned TOML; unknown keys and invalid combinations are rejected. Geometry, STL coordinates and roughness are in **nm**, time in **ps**, temperature in **K**, power density in **W/m³**, total power in **W**, and the imposed gradient in **K/m**.

| RTA formulation | `background_temperature` | `lifetime_temperature` |
|---|---:|---|
| Full distribution | `0` | `"local"` |
| Deviational, fixed scattering | `300` | `300` |
| Deviational, temperature-dependent scattering | `300` | `"local"` |

Holding rates at 300 K does not by itself linearize the Bose occupation update. Periodic-gradient driving is a separate linear-response formulation. See [methodology](docs/source/methodology.rst), [input reference](docs/source/input_files.rst), and [gradient driving](docs/source/temperature_gradient.rst).

## Results

Each run creates an indexed output directory with `summary.txt`, `grid_centers.csv` and `convergence.txt`. Gradient runs additionally write `grid_heat_flux.csv`.

- Reservoir-driven `kappa_eff` uses the reservoir temperature difference; `kappa_int` uses a fitted cell-center gradient.
- Gradient-driven `kappa_eff` is `-<q>/G`; `kappa_int` is not applicable.
- Localized-heating cases report temperatures and heat flow; a single conductivity is not automatically a device material property.

```sh
python3 tools/plot_convergence.py example/results/Cross_100nm_0
python3 tools/plot_temperature_3d.py --input example/input_finfet_stl_heat1e20.toml --results example/results/FinFET_heat1e20_0
```

Use the directory actually printed by the run. Plotting requires NumPy and Matplotlib. Check energy balance, averaging windows, time step, spatial resolution and carrier-count convergence before interpreting results.

## Documentation and source layout

```sh
python3 -m pip install -r docs/requirements.txt
python3 -m sphinx -W --keep-going -b html docs/source docs/build/html
```

Open `docs/build/html/index.html`.

| Directory | Purpose |
|---|---|
| `include/`, `src/` | Material data, transport, collisions, sources, interfaces and output |
| `example/`, `model/` | Curated input examples and tutorial geometries |
| `material/`, `tools/` | Bundled tutorial data and preparation/analysis utilities |
| `docs/source/` | Theory, configuration, tutorials and architecture |

Local `tests/`, `calculation/`, `analysis/`, `output/`, manuscript assets, generated example results and build products are excluded from Git. They remain on the local filesystem. The public source tree builds without these folders. See [repository scope](docs/source/repository.rst) and [architecture](docs/source/architecture.rst).

## 中文说明

PhonoMC 使用第一性原理声子数据和固定参考温度的偏差蒙特卡洛方法，计算薄膜及器件中的声子热输运。主要功能包括温度相关 RTA、Callaway 双弛豫、有限/无限长度面内输运、多材料 DMM/AMM/MMM 界面，以及空间和频率可加权的局域热源。

参考温度固定，散射率可按局域温度查询；碰撞与热源在现有载体上更新占据数。开放热浴的吸收和注入仍可引起总载体数波动。对外输入仅保留 `rta` 和 `callaway`，全散射矩阵仍为内部实验功能。

按上面的命令编译并运行示例。测试、计算、结果及论文文件夹仅本地保留，不随代码提交；公开材料保留 Si、Ge 和 SiC；SiC 为专用 NEP 数据，包含 N/U 和同位素散射率。Graphene 不随代码上传，其他计算材料目录继续忽略。当前整理范围与使用限制见 [仓库说明](docs/source/repository.rst)。

## Русский

PhonoMC — C++17 программа для моделирования фононного теплопереноса методом Монте-Карло с фиксированным равновесным фоном и температурно-зависимыми скоростями рассеяния. Поддерживаются RTA, модель Callaway, тонкие плёнки, локальные источники тепла и многослойные структуры. Инструкции по сборке и примеры приведены выше; подробная документация — в `docs/source/`. Локальные расчёты, результаты и тесты не входят в публикуемое дерево исходного кода.
