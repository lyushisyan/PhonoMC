# PhonoMC

[English](#english-version) | [中文](#中文版) | [Русский](#русский)

<a id="english-version"></a>
## English Version

---

### Overview

PhonoMC is a C++ phonon Monte Carlo simulator for semiconductor heat transport (e.g., Si/Ge/GaN). It supports box/mesh geometries, HDF5 material data, rough/periodic boundaries, and thermal conductivity estimation.

### Core Features

- Phonon particle transport with boundary scattering and mode updates.
- Geometry handling with surface mesh processing and grid partitioning.
- Real material loading from HDF5 (strict by default; invalid paths abort run).
- Parallel execution with OpenMP.
- Outputs for temperature profile, heat flux, and conductivity (`kappa_int`, `kappa_eff`).

### Requirements (Linux)

- Compiler: GCC 9+ / Clang 10+
- CMake 3.20+
- HDF5 1.10+ (`libhdf5-dev`)
- Qhull (`qhull-bin`, and `libqhull-dev` if linking headers/libs)
- OpenMP (`libomp-dev` on Ubuntu)

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPHONOMC_ENABLE_OPENMP=ON
cmake --build build -j
```

Executable: `build/PhonoMC`

Run the basic configuration tests:

```bash
ctest --test-dir build --output-on-failure
```

### Documentation

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html
```

Open `docs/build/html/index.html`. Read the Docs is configured through
`.readthedocs.yaml`.

### Run

```bash
./build/PhonoMC example/input_cross_100nm.toml
```

Optional thread setting:

```bash
OMP_NUM_THREADS=64 ./build/PhonoMC example/input_cross_100nm.toml
```

### Input Rules (Current)

- Only sectioned `.toml` input is supported; legacy `.txt` / `--option` input and top-level key aliases are not accepted.
- Keys must be placed under `[geometry]`, `[simulation]`, `[boundary]`, `[heat_source]`, or `[io]`; unknown keys abort parsing.
- Geometry supports `model = "box"` or a mesh file path; the cylinder generator is not supported.
- Required simulation grid key: `grid_xyz = [nx, ny, nz]`.
- `boundary_position`, `boundary_conditions`, and `boundary_values` must have identical lengths. Use value `0` for `P`, and list every periodic region exactly once in `periodic_pair`.
- Every declared boundary region must match at least one mesh facet, and every periodic facet must be paired.
- Temperature export stride: `convergence_write_interval = 10` (write `convergence.txt` every N steps).
- Reproducible runs: set `random_seed = 12345` under `[simulation]`; exact OpenMP replay also requires the same thread configuration.
- Progress print mode: `progress_temperature_summary_only = true` prints only `Tmin/Tavg/Tmax`.
- `initial_temperature` supports:
  - `300` (uniform initial temperature, unit: K)
  - `"linear"` (linear profile from cold reservoir to hot reservoir)
- Temperature reference controls use an invariant numeric `background_temperature`; the removed `"local"` background is rejected:
  - full-phonon sampling: `background_temperature = 0`, `lifetime_temperature = "local"`.
  - fixed-reference linearized sampling: `background_temperature = 300`, `lifetime_temperature = 300`.
  - fixed-reference sampling with local lifetimes: `background_temperature = 300`, `lifetime_temperature = "local"`.
- Reservoir refill uses `one_to_one` (inject by previous-step leaving counts).
- If no valid thermal reservoirs are present, fallback temperature range is `299/301 K`.
- `sizes` in TOML use **nm**.
- STL coordinates use **nm**.
- Rough boundary value (`R` in `boundary_conditions`) in `boundary_values` uses **nm**.
- Both POSCAR and HDF5 material data are required; either loading failure aborts the run.

### Outputs

- `convergence.txt`
  - Columns include: `timestep`, `time_ps`, `T_1 ... T_n`, `heatflux`, `kappa_int`, `kappa_eff`, particle/reservoir balance, source energy, lifetime-scattering residual, step energy-balance residual, and total thermal energy.
- `summary.txt`
  - Consolidated input/geometry/grid/boundary/runtime summary.
- `grid_centers.csv`
  - Grid center coordinates in **nm** (`x_nm,y_nm,z_nm`) and `volume_nm3`.

### Plotting

1D cross-plane/in-plane results:

```bash
python3 tools/plot_convergence.py example/results/Cross_100nm_0
python3 tools/plot_convergence.py example/results/Inplane_x10000nm_z100nm_r1nm_0
```

Each run writes three figures under `plots_1d/`: `temperature_vs_time.png`, `heatflux_vs_time.png`, and `kappa_vs_time.png`.

3D FinFET results:

```bash
python3 tools/plot_temperature_3d.py \
  --input example/input_finfet_stl_heat1e20.toml \
  --results example/results/FinFET_heat1e20_0
```

This writes `temperature_3d.png`, `temperature_slice_xrel0.500_yz.png`, and `temperature_slice_yrel0.500_xz.png` under `plots_3d/`. Slice locations can be changed with `--x-slice-rel` and `--y-slice-rel`.

---

<a id="中文版"></a>
## 中文版

---

### 项目简介

PhonoMC 是一个基于 C++ 的声子蒙特卡洛仿真程序，用于半导体（如 Si/Ge/GaN）热输运计算。支持 box/网格模型、HDF5 材料数据、粗糙/周期边界和热导率评估。

### 核心功能

- 声子粒子输运、边界散射与模式更新。
- 几何处理与网格划分。
- 从 HDF5 加载真实材料数据（默认严格模式，路径错误会直接终止）。
- OpenMP 并行计算。
- 输出温度、热流和热导率（`kappa_int`、`kappa_eff`）。

### 环境依赖（Linux）

- 编译器：GCC 9+ / Clang 10+
- CMake 3.20+
- HDF5 1.10+（`libhdf5-dev`）
- Qhull（`qhull-bin`；若需链接头文件/库再安装 `libqhull-dev`）
- OpenMP（Ubuntu 推荐 `libomp-dev`）

### 编译

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPHONOMC_ENABLE_OPENMP=ON
cmake --build build -j
```

可执行文件：`build/PhonoMC`

运行基础配置测试：

```bash
ctest --test-dir build --output-on-failure
```

### 文档

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html
```

打开 `docs/build/html/index.html`。Read the Docs 通过 `.readthedocs.yaml` 配置。

### 运行

```bash
./build/PhonoMC example/input_cross_100nm.toml
```

可指定线程数：

```bash
OMP_NUM_THREADS=64 ./build/PhonoMC example/input_cross_100nm.toml
```

### 当前输入规则

- 仅支持分节 `.toml` 输入；不再支持旧 `.txt` / `--option` 格式和顶层键别名。
- 配置键必须位于 `[geometry]`、`[simulation]`、`[boundary]`、`[heat_source]` 或 `[io]` 中；未知键会直接报错。
- 几何仅支持 `model = "box"` 或网格文件路径，不再支持圆柱生成器。
- 网格参数必填：`grid_xyz = [nx, ny, nz]`。
- `boundary_position`、`boundary_conditions` 和 `boundary_values` 长度必须一致；`P` 边界的值填 `0`，并在 `periodic_pair` 中恰好列出一次。
- 每个边界区域必须至少命中一个 facet，每个周期 facet 必须成对。
- 温度导出步长：`convergence_write_interval = 10`（每 N 步写一次 `convergence.txt`）。
- 可重现运行：在 `[simulation]` 中设置 `random_seed = 12345`；OpenMP 精确重现还需保持相同线程配置。
- 进度打印模式：`progress_temperature_summary_only = true` 时只打印 `Tmin/Tavg/Tmax`。
- `initial_temperature` 支持：
  - `300`（全域统一初始温度，单位 K）
  - `"linear"`（按冷热库方向线性初始化）
- 温度参考必须使用数值型固定 `background_temperature`；已移除并拒绝 `"local"` 背景：
  - 全声子：`background_temperature = 0`，`lifetime_temperature = "local"`。
  - 偏差 + 固定参考/寿命：`background_temperature = 300`，`lifetime_temperature = 300`。
  - 偏差 + 固定背景/局部寿命：`background_temperature = 300`，`lifetime_temperature = "local"`。
- 热库回填固定为 `one_to_one`（按上一步离开热库的粒子数回填）。
- 若没有可用热库，默认回退温度范围为 `299/301 K`。
- TOML 中 `sizes` 单位是 **nm**。
- STL 坐标单位是 **nm**。
- 粗糙边界（`boundary_conditions` 中的 `R`）在 `boundary_values` 里单位是 **nm**。
- POSCAR 和 HDF5 材料数据都必须存在且有效，任一加载失败都会终止运行。

### 输出文件

- `convergence.txt`
  - 包含：`timestep`、`time_ps`、`T_1 ... T_n`、`heatflux`、`kappa_int`、`kappa_eff`、粒子/热库收支、热源能量、寿命散射残差与总热能。
- `summary.txt`
  - 汇总关键输入与几何/网格/边界/运行时间信息。
- `grid_centers.csv`
  - 网格中心坐标（**nm**）：`x_nm,y_nm,z_nm`，以及 `volume_nm3`。

### 绘图

1D cross-plane / in-plane 结果：

```bash
python3 tools/plot_convergence.py example/results/Cross_100nm_0
python3 tools/plot_convergence.py example/results/Inplane_x10000nm_z100nm_r1nm_0
```

每次会在对应结果目录的 `plots_1d/` 中输出三张图：`temperature_vs_time.png`、`heatflux_vs_time.png`、`kappa_vs_time.png`。

3D FinFET 结果：

```bash
python3 tools/plot_temperature_3d.py \
  --input example/input_finfet_stl_heat1e20.toml \
  --results example/results/FinFET_heat1e20_0
```

会在 `plots_3d/` 中输出：`temperature_3d.png`、`temperature_slice_xrel0.500_yz.png`、`temperature_slice_yrel0.500_xz.png`。截面位置可用 `--x-slice-rel` 和 `--y-slice-rel` 调整。

---

<a id="русский"></a>
## Русский

---

### Обзор

PhonoMC — C++ симулятор фононного Монте-Карло для теплопереноса в полупроводниках (например, Si/Ge/GaN). Поддерживаются box/mesh геометрии, HDF5-материалы, шероховатые/периодические границы и оценка теплопроводности.

### Основные возможности

- Транспорт фононных частиц, рассеяние на границах и обновление мод.
- Обработка геометрии и разбиение на сетку.
- Загрузка реальных материалов из HDF5 (по умолчанию строгий режим, при ошибке путь/файл — остановка).
- Параллельный запуск через OpenMP.
- Вывод профиля температуры, теплового потока и теплопроводности (`kappa_int`, `kappa_eff`).

### Зависимости (Linux)

- Компилятор: GCC 9+ / Clang 10+
- CMake 3.20+
- HDF5 1.10+ (`libhdf5-dev`)
- Qhull (`qhull-bin`, а `libqhull-dev` нужен для линковки заголовков/библиотек)
- OpenMP (`libomp-dev` для Ubuntu)

### Сборка

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPHONOMC_ENABLE_OPENMP=ON
cmake --build build -j
```

Исполняемый файл: `build/PhonoMC`

Запуск базовых тестов конфигурации:

```bash
ctest --test-dir build --output-on-failure
```

### Документация

```bash
python3 -m pip install -r docs/requirements.txt
make -C docs html
```

Откройте `docs/build/html/index.html`. Read the Docs настраивается через
`.readthedocs.yaml`.

### Запуск

```bash
./build/PhonoMC example/input_cross_100nm.toml
```

Ограничение потоков:

```bash
OMP_NUM_THREADS=64 ./build/PhonoMC example/input_cross_100nm.toml
```

### Правила входных данных (текущие)

- Поддерживаются только секционные файлы `.toml`; старые `.txt` / `--option` и псевдонимы ключей верхнего уровня удалены.
- Ключи должны находиться в `[geometry]`, `[simulation]`, `[boundary]`, `[heat_source]` или `[io]`; неизвестные ключи вызывают ошибку.
- Геометрия поддерживает `model = "box"` или путь к сетке; генератор цилиндра удалён.
- Обязательный ключ сетки: `grid_xyz = [nx, ny, nz]`.
- `boundary_position`, `boundary_conditions` и `boundary_values` должны иметь одинаковую длину; для `P` используется значение `0`, а все периодические области ровно один раз указываются в `periodic_pair`.
- Каждая область границы должна совпасть хотя бы с одним facet, а каждая периодическая грань должна иметь пару.
- Шаг вывода температуры: `convergence_write_interval = 10` (запись `convergence.txt` каждые N шагов).
- Воспроизводимость: задайте `random_seed = 12345` в `[simulation]`; для точного повтора OpenMP также нужна та же конфигурация потоков.
- Режим вывода прогресса: `progress_temperature_summary_only = true` печатает только `Tmin/Tavg/Tmax`.
- `initial_temperature` поддерживает:
  - `300` (равномерная начальная температура, K)
  - `"linear"` (линейный профиль между холодным и горячим резервуарами)
- Температура фона всегда задаётся фиксированным числом; `background_temperature = "local"` больше не поддерживается:
  - полные фононы: `background_temperature = 0`, `lifetime_temperature = "local"`.
  - фиксированные фон и время жизни: `background_temperature = 300`, `lifetime_temperature = 300`.
  - фиксированный фон и локальное время жизни: `background_temperature = 300`, `lifetime_temperature = "local"`.
- Пополнение резервуаров фиксировано как `one_to_one` (инжекция по числу частиц, покинувших резервуары на предыдущем шаге).
- Если валидные терморезервуары отсутствуют, используется диапазон по умолчанию `299/301 K`.
- `sizes` в TOML задаются в **nm**.
- STL-координаты задаются в **nm**.
- Для шероховатой границы (`R` в `boundary_conditions`) значение в `boundary_values` задается в **nm**.
- POSCAR и HDF5 обязательны; ошибка загрузки любого из них прерывает расчёт.

### Выходные файлы

- `convergence.txt`
  - Колонки: `timestep`, `time_ps`, `T_1 ... T_n`, `heatflux`, `kappa_int`, `kappa_eff`, баланс частиц/резервуаров, энергия источника, остаток релаксации и полная тепловая энергия.
- `summary.txt`
  - Сводка ключевых входных параметров и статистики geometry/grid/boundary/runtime.
- `grid_centers.csv`
  - Координаты центров сетки в **nm** (`x_nm,y_nm,z_nm`) и `volume_nm3`.

### Построение графиков

1D результаты cross-plane/in-plane:

```bash
python3 tools/plot_convergence.py example/results/Cross_100nm_0
python3 tools/plot_convergence.py example/results/Inplane_x10000nm_z100nm_r1nm_0
```

Команда создаёт три файла в `plots_1d/`: `temperature_vs_time.png`, `heatflux_vs_time.png` и `kappa_vs_time.png`.

3D результаты FinFET:

```bash
python3 tools/plot_temperature_3d.py \
  --input example/input_finfet_stl_heat1e20.toml \
  --results example/results/FinFET_heat1e20_0
```

Команда создаёт `temperature_3d.png`, `temperature_slice_xrel0.500_yz.png` и `temperature_slice_yrel0.500_xz.png` в `plots_3d/`. Положение срезов задаётся параметрами `--x-slice-rel` и `--y-slice-rel`.
