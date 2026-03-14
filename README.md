# Sequence-MC

[English](#english-version) | [中文](#中文版) | [Русский](#русский)

<a id="english-version"></a>
## English Version

---

### Overview

Sequence-MC (NTMC) is a C++ phonon Monte Carlo simulator for semiconductor heat transport (e.g., Si/Ge/GaN). It supports box/mesh geometries, HDF5 material data, rough/periodic boundaries, and thermal conductivity estimation.

### Core Features

- Phonon particle transport with boundary scattering and mode updates.
- Geometry handling with surface mesh processing and grid partitioning.
- Real material loading from HDF5 (strict by default; invalid paths abort run).
- Parallel execution with OpenMP.
- Outputs for temperature profile, heat flux, and conductivity (`kappa_fit`, `kappa_end`).

### Requirements (Linux)

- Compiler: GCC 9+ / Clang 10+
- CMake 3.20+
- HDF5 1.10+ (`libhdf5-dev`)
- Qhull (`qhull-bin`, and `libqhull-dev` if linking headers/libs)
- OpenMP (`libomp-dev` on Ubuntu)

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNTMC_ENABLE_OPENMP=ON
cmake --build build -j
```

Executable: `build/ntmc`

### Run

```bash
./build/ntmc input_cross_100nm.toml
```

Optional thread setting:

```bash
OMP_NUM_THREADS=64 ./build/ntmc input_cross_100nm.toml
```

### Input Rules (Current)

- Strict key set is expected.
- Required simulation grid key: `grid_xyz = [nx, ny, nz]`.
- `initial_temperature` supports only: `cold`, `mean`.
- `sizes` in TOML are in **nm**; internally converted to **Angstrom**.
- STL coordinates are read as **nm**; internally converted to **Angstrom**.
- Rough boundary value (`R` in `boundary_conditions`) is entered in **nm**, internally multiplied by 10.
- If HDF5 material loading fails, run aborts by default.
  - Test-only fallback is allowed with: `NTMC_ALLOW_SYNTHETIC_MATERIAL=1`.

### Outputs

- `convergence.txt`
  - Columns include: `timestep`, `time_ps`, `T_1 ... T_n`, `heatflux`, `kappa_fit`, `kappa_end`.
- `summary.txt`
  - Consolidated input/geometry/grid/boundary/runtime summary.
- `grid_centers.csv`
  - Grid center coordinates in **nm** (`x_nm,y_nm,z_nm`) and `volume_nm3`.

---

<a id="中文版"></a>
## 中文版

---

### 项目简介

Sequence-MC（NTMC）是一个基于 C++ 的声子蒙特卡洛仿真程序，用于半导体（如 Si/Ge/GaN）热输运计算。支持 box/网格模型、HDF5 材料数据、粗糙/周期边界和热导率评估。

### 核心功能

- 声子粒子输运、边界散射与模式更新。
- 几何处理与网格划分。
- 从 HDF5 加载真实材料数据（默认严格模式，路径错误会直接终止）。
- OpenMP 并行计算。
- 输出温度、热流和热导率（`kappa_fit`、`kappa_end`）。

### 环境依赖（Linux）

- 编译器：GCC 9+ / Clang 10+
- CMake 3.20+
- HDF5 1.10+（`libhdf5-dev`）
- Qhull（`qhull-bin`；若需链接头文件/库再安装 `libqhull-dev`）
- OpenMP（Ubuntu 推荐 `libomp-dev`）

### 编译

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNTMC_ENABLE_OPENMP=ON
cmake --build build -j
```

可执行文件：`build/ntmc`

### 运行

```bash
./build/ntmc input_cross_100nm.toml
```

可指定线程数：

```bash
OMP_NUM_THREADS=64 ./build/ntmc input_cross_100nm.toml
```

### 当前输入规则

- 采用严格关键词集合。
- 网格参数必填：`grid_xyz = [nx, ny, nz]`。
- `initial_temperature` 仅支持：`cold`、`mean`。
- TOML 中 `sizes` 单位是 **nm**，程序内部转换为 **Å**。
- STL 坐标按 **nm** 读取，内部转换为 **Å**。
- 粗糙边界（`boundary_conditions` 中的 `R`）在 `boundary_values` 里按 **nm** 填写，内部乘 10。
- HDF5 材料加载失败时默认直接报错终止。
  - 仅测试场景可设：`NTMC_ALLOW_SYNTHETIC_MATERIAL=1` 启用回退。

### 输出文件

- `convergence.txt`
  - 包含：`timestep`、`time_ps`、`T_1 ... T_n`、`heatflux`、`kappa_fit`、`kappa_end`。
- `summary.txt`
  - 汇总关键输入与几何/网格/边界/运行时间信息。
- `grid_centers.csv`
  - 网格中心坐标（**nm**）：`x_nm,y_nm,z_nm`，以及 `volume_nm3`。

---

<a id="русский"></a>
## Русский

---

### Обзор

Sequence-MC (NTMC) — C++ симулятор фононного Монте-Карло для теплопереноса в полупроводниках (например, Si/Ge/GaN). Поддерживаются box/mesh геометрии, HDF5-материалы, шероховатые/периодические границы и оценка теплопроводности.

### Основные возможности

- Транспорт фононных частиц, рассеяние на границах и обновление мод.
- Обработка геометрии и разбиение на сетку.
- Загрузка реальных материалов из HDF5 (по умолчанию строгий режим, при ошибке путь/файл — остановка).
- Параллельный запуск через OpenMP.
- Вывод профиля температуры, теплового потока и теплопроводности (`kappa_fit`, `kappa_end`).

### Зависимости (Linux)

- Компилятор: GCC 9+ / Clang 10+
- CMake 3.20+
- HDF5 1.10+ (`libhdf5-dev`)
- Qhull (`qhull-bin`, а `libqhull-dev` нужен для линковки заголовков/библиотек)
- OpenMP (`libomp-dev` для Ubuntu)

### Сборка

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DNTMC_ENABLE_OPENMP=ON
cmake --build build -j
```

Исполняемый файл: `build/ntmc`

### Запуск

```bash
./build/ntmc input_cross_100nm.toml
```

Ограничение потоков:

```bash
OMP_NUM_THREADS=64 ./build/ntmc input_cross_100nm.toml
```

### Правила входных данных (текущие)

- Используется строгий набор ключей.
- Обязательный ключ сетки: `grid_xyz = [nx, ny, nz]`.
- `initial_temperature` поддерживает только: `cold`, `mean`.
- `sizes` в TOML задаются в **nm**, внутри переводятся в **Å**.
- STL-координаты читаются как **nm**, внутри переводятся в **Å**.
- Для шероховатой границы (`R` в `boundary_conditions`) значение в `boundary_values` задается в **nm**, внутри умножается на 10.
- При ошибке загрузки HDF5 расчёт по умолчанию прерывается.
  - Тестовый fallback включается только при `NTMC_ALLOW_SYNTHETIC_MATERIAL=1`.

### Выходные файлы

- `convergence.txt`
  - Колонки: `timestep`, `time_ps`, `T_1 ... T_n`, `heatflux`, `kappa_fit`, `kappa_end`.
- `summary.txt`
  - Сводка ключевых входных параметров и статистики geometry/grid/boundary/runtime.
- `grid_centers.csv`
  - Координаты центров сетки в **nm** (`x_nm,y_nm,z_nm`) и `volume_nm3`.
