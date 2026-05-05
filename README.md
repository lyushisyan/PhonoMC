# EPMC

[English](#english-version) | [中文](#中文版) | [Русский](#русский)

<a id="english-version"></a>
## English Version

---

### Overview

EPMC is a C++ phonon Monte Carlo simulator for semiconductor heat transport (e.g., Si/Ge/GaN). It supports box/mesh geometries, HDF5 material data, rough/periodic boundaries, and thermal conductivity estimation.

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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEPMC_ENABLE_OPENMP=ON
cmake --build build -j
```

Executable: `build/EPMC`

### Run

```bash
./build/EPMC input_cross_100nm.toml
```

Optional thread setting:

```bash
OMP_NUM_THREADS=64 ./build/EPMC input_cross_100nm.toml
```

### Input Rules (Current)

- Strict key set is expected.
- Required simulation grid key: `grid_xyz = [nx, ny, nz]`.
- Temperature export stride: `convergence_write_interval = 10` (write `convergence.txt` every N steps).
- Progress print mode: `progress_temperature_summary_only = true` prints only `Tmin/Tavg/Tmax`.
- `initial_temperature` supports:
  - `300` (uniform initial temperature, unit: K)
  - `"linear"` (linear profile from cold reservoir to hot reservoir)
- Temperature reference controls:
  - default: `background_temperature_mode = "local"` and `lifetime_temperature_mode = "local"` use `E(Tgrid)` and `tau(Tgrid)`.
  - fixed background with local lifetime: `background_temperature_mode = "fixed"`, `background_temperature = 300`, `lifetime_temperature_mode = "local"`.
  - fixed background and fixed lifetime: also set `lifetime_temperature_mode = "fixed"`, `lifetime_temperature = 300`.
- Reservoir refill uses `one_to_one` (inject by previous-step leaving counts).
- If no valid thermal reservoirs are present, fallback temperature range is `299/301 K`.
- `sizes` in TOML use **nm**.
- STL coordinates use **nm**.
- Rough boundary value (`R` in `boundary_conditions`) in `boundary_values` uses **nm**.
- If HDF5 material loading fails, run aborts by default.

### Outputs

- `convergence.txt`
  - Columns include: `timestep`, `time_ps`, `T_1 ... T_n`, `heatflux`, `kappa_fit`, `kappa_end`, `absorbed`, `injected`, `recovered`, `net`.
- `summary.txt`
  - Consolidated input/geometry/grid/boundary/runtime summary.
- `grid_centers.csv`
  - Grid center coordinates in **nm** (`x_nm,y_nm,z_nm`) and `volume_nm3`.

---

<a id="中文版"></a>
## 中文版

---

### 项目简介

EPMC 是一个基于 C++ 的声子蒙特卡洛仿真程序，用于半导体（如 Si/Ge/GaN）热输运计算。支持 box/网格模型、HDF5 材料数据、粗糙/周期边界和热导率评估。

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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEPMC_ENABLE_OPENMP=ON
cmake --build build -j
```

可执行文件：`build/EPMC`

### 运行

```bash
./build/EPMC input_cross_100nm.toml
```

可指定线程数：

```bash
OMP_NUM_THREADS=64 ./build/EPMC input_cross_100nm.toml
```

### 当前输入规则

- 采用严格关键词集合。
- 网格参数必填：`grid_xyz = [nx, ny, nz]`。
- 温度导出步长：`convergence_write_interval = 10`（每 N 步写一次 `convergence.txt`）。
- 进度打印模式：`progress_temperature_summary_only = true` 时只打印 `Tmin/Tavg/Tmax`。
- `initial_temperature` 支持：
  - `300`（全域统一初始温度，单位 K）
  - `"linear"`（按冷热库方向线性初始化）
- 温度参考控制：
  - 默认：`background_temperature_mode = "local"` 且 `lifetime_temperature_mode = "local"`，即使用 `E(Tgrid)` 和 `tau(Tgrid)`。
  - 固定背景、局部寿命：`background_temperature_mode = "fixed"`，`background_temperature = 300`，`lifetime_temperature_mode = "local"`。
  - 固定背景、固定寿命：再设置 `lifetime_temperature_mode = "fixed"`，`lifetime_temperature = 300`。
- 热库回填固定为 `one_to_one`（按上一步离开热库的粒子数回填）。
- 若没有可用热库，默认回退温度范围为 `299/301 K`。
- TOML 中 `sizes` 单位是 **nm**。
- STL 坐标单位是 **nm**。
- 粗糙边界（`boundary_conditions` 中的 `R`）在 `boundary_values` 里单位是 **nm**。
- HDF5 材料加载失败时默认直接报错终止。

### 输出文件

- `convergence.txt`
  - 包含：`timestep`、`time_ps`、`T_1 ... T_n`、`heatflux`、`kappa_fit`、`kappa_end`、`absorbed`、`injected`、`recovered`、`net`。
- `summary.txt`
  - 汇总关键输入与几何/网格/边界/运行时间信息。
- `grid_centers.csv`
  - 网格中心坐标（**nm**）：`x_nm,y_nm,z_nm`，以及 `volume_nm3`。

---

<a id="русский"></a>
## Русский

---

### Обзор

EPMC — C++ симулятор фононного Монте-Карло для теплопереноса в полупроводниках (например, Si/Ge/GaN). Поддерживаются box/mesh геометрии, HDF5-материалы, шероховатые/периодические границы и оценка теплопроводности.

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
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DEPMC_ENABLE_OPENMP=ON
cmake --build build -j
```

Исполняемый файл: `build/EPMC`

### Запуск

```bash
./build/EPMC input_cross_100nm.toml
```

Ограничение потоков:

```bash
OMP_NUM_THREADS=64 ./build/EPMC input_cross_100nm.toml
```

### Правила входных данных (текущие)

- Используется строгий набор ключей.
- Обязательный ключ сетки: `grid_xyz = [nx, ny, nz]`.
- Шаг вывода температуры: `convergence_write_interval = 10` (запись `convergence.txt` каждые N шагов).
- Режим вывода прогресса: `progress_temperature_summary_only = true` печатает только `Tmin/Tavg/Tmax`.
- `initial_temperature` поддерживает:
  - `300` (равномерная начальная температура, K)
  - `"linear"` (линейный профиль между холодным и горячим резервуарами)
- Настройки температурной привязки:
  - по умолчанию: `background_temperature_mode = "local"` и `lifetime_temperature_mode = "local"` используют `E(Tgrid)` и `tau(Tgrid)`.
  - фиксированный фон и локальное время жизни: `background_temperature_mode = "fixed"`, `background_temperature = 300`, `lifetime_temperature_mode = "local"`.
  - фиксированный фон и фиксированное время жизни: дополнительно `lifetime_temperature_mode = "fixed"`, `lifetime_temperature = 300`.
- Пополнение резервуаров фиксировано как `one_to_one` (инжекция по числу частиц, покинувших резервуары на предыдущем шаге).
- Если валидные терморезервуары отсутствуют, используется диапазон по умолчанию `299/301 K`.
- `sizes` в TOML задаются в **nm**.
- STL-координаты задаются в **nm**.
- Для шероховатой границы (`R` в `boundary_conditions`) значение в `boundary_values` задается в **nm**.
- При ошибке загрузки HDF5 расчёт по умолчанию прерывается.

### Выходные файлы

- `convergence.txt`
  - Колонки: `timestep`, `time_ps`, `T_1 ... T_n`, `heatflux`, `kappa_fit`, `kappa_end`, `absorbed`, `injected`, `recovered`, `net`.
- `summary.txt`
  - Сводка ключевых входных параметров и статистики geometry/grid/boundary/runtime.
- `grid_centers.csv`
  - Координаты центров сетки в **nm** (`x_nm,y_nm,z_nm`) и `volume_nm3`.
