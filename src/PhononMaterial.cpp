#include "PhononMaterial.h"

#include <H5Cpp.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace {
using Mat3 = std::array<std::array<double, 3>, 3>;
using Vec3 = std::array<double, 3>;

// 函数说明：读取 HDF5 数值数据集并展开为一维数组。
std::vector<double> read_dataset_nd_double(const H5::H5File& file, const std::string& name, std::vector<hsize_t>& dims_out) {
    H5::DataSet ds = file.openDataSet(name);
    H5::DataSpace sp = ds.getSpace();
    const int nd = sp.getSimpleExtentNdims();
    if (nd <= 0) {
        throw std::runtime_error("Dataset " + name + " has invalid rank");
    }
    dims_out.assign(static_cast<size_t>(nd), 0);
    sp.getSimpleExtentDims(dims_out.data());
    size_t n = 1;
    for (hsize_t d : dims_out) {
        n *= static_cast<size_t>(d);
    }
    std::vector<double> out(n);
    ds.read(out.data(), H5::PredType::NATIVE_DOUBLE);
    return out;
}

// 函数说明：清理输入文本两端空白，保证配置与数据解析的稳健性。
std::string trim(const std::string& s) {
    const auto b = std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c) != 0; });
    if (b == s.end()) {
        return "";
    }
    const auto e = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c) != 0; }).base();
    return std::string(b, e);
}

// 函数说明：由散射率 gamma 换算模态寿命 tau。
double tau_from_gamma(double gamma) {
    if (gamma <= 0.0 || std::isnan(gamma) || std::isinf(gamma)) {
        return 0.0;
    }
    return 1.0 / (4.0 * M_PI * gamma);  // ps
}

// 函数说明：计算 3x3 矩阵逆，用于晶格变换或四面体坐标变换。
bool inverse3x3(const Mat3& a, Mat3& inv_out) {
    const double det =
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if (std::abs(det) <= 1e-20) {
        return false;
    }
    const double id = 1.0 / det;
    inv_out[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * id;
    inv_out[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * id;
    inv_out[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * id;
    inv_out[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * id;
    inv_out[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * id;
    inv_out[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * id;
    inv_out[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * id;
    inv_out[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * id;
    inv_out[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * id;
    return true;
}

// 函数说明：将分数 q 点映射到倒易空间 k 向量。
Vec3 q_to_k(const Mat3& reciprocal_lattice, const Vec3& q) {
    return {
        reciprocal_lattice[0][0] * q[0] + reciprocal_lattice[0][1] * q[1] + reciprocal_lattice[0][2] * q[2],
        reciprocal_lattice[1][0] * q[0] + reciprocal_lattice[1][1] * q[1] + reciprocal_lattice[1][2] * q[2],
        reciprocal_lattice[2][0] * q[0] + reciprocal_lattice[2][1] * q[1] + reciprocal_lattice[2][2] * q[2]
    };
}

// 函数说明：计算三维向量模长。
double vec_norm(const Vec3& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
}  // namespace

// 函数说明：加载材料声子数据并建立温度-能量查找表。
PhononMaterial::PhononMaterial(const SimulationConfig& args, int mat_index) {
    if (std::isfinite(args.temperature_lookup_dt) && args.temperature_lookup_dt > 0.0) {
        temperature_lookup_dt_ = args.temperature_lookup_dt;
    }
    std::string err;
    if (!load_hdf5_data(args, mat_index, &err)) {
        const char* allow_fallback = std::getenv("EPMC_ALLOW_SYNTHETIC_MATERIAL");
        const bool use_synthetic = (allow_fallback != nullptr) && (std::string(allow_fallback) == "1");
        if (!use_synthetic) {
            throw std::runtime_error(
                "Failed to load HDF5 phonon data (" + err + "). "
                "Simulation is aborted to avoid invalid thermal conductivity results. "
                "Fix material_folder / HDF5 / POSCAR paths, or set EPMC_ALLOW_SYNTHETIC_MATERIAL=1 for test-only fallback.");
        }
        std::cerr << "Warning: failed to load HDF5 phonon data (" << err
                  << "). EPMC_ALLOW_SYNTHETIC_MATERIAL=1 is set, falling back to synthetic mode bank.\n";
        build_fallback_modes(args);
    }
    initialize_temperature_lookup();
    std::cout << "PhononMaterial initialized: active_mode_count=" << active_mode_count_ << '\n';
}

// 函数说明：在离散采样点中查找与目标值最近的索引。
int PhononMaterial::nearest_index(const std::vector<double>& arr, double x) {
    if (arr.empty()) {
        return 0;
    }
    int best = 0;
    double best_abs = std::abs(arr[0] - x);
    for (int i = 1; i < static_cast<int>(arr.size()); ++i) {
        const double a = std::abs(arr[static_cast<size_t>(i)] - x);
        if (a < best_abs) {
            best_abs = a;
            best = i;
        }
    }
    return best;
}

// 函数说明：将 (q, branch) 模态索引压平到线性索引。
int PhononMaterial::flatten_mode_index(const Mode& mode) const {
    return mode[0] * branch_count_ + mode[1];
}

// 函数说明：按活跃模态索引返回具体 (q, branch) 模式。
PhononMaterial::Mode PhononMaterial::active_mode_at(int active_index) const {
    if (active_index < 0 || active_index >= static_cast<int>(active_mode_list_.size())) {
        return {0, 0};
    }
    return active_mode_list_[static_cast<size_t>(active_index)];
}

// 函数说明：将模态反查为活跃模态索引。
int PhononMaterial::active_index_for_mode(const Mode& mode) const {
    const int fi = flatten_mode_index(mode);
    if (fi < 0 || fi >= static_cast<int>(flat_to_active_index_.size())) {
        return -1;
    }
    return flat_to_active_index_[static_cast<size_t>(fi)];
}

// 函数说明：查询简并模态对应的伴随支。
int PhononMaterial::degenerate_partner_branch(const Mode& mode) const {
    const int fi = flatten_mode_index(mode);
    if (fi < 0 || fi >= static_cast<int>(degenerate_partner_branch_data_.size())) {
        return -1;
    }
    return degenerate_partner_branch_data_[static_cast<size_t>(fi)];
}

// 函数说明：执行线性插值，用于温度与散射参数之间的连续过渡。
double PhononMaterial::lerp(double x0, double x1, double y0, double y1, double x) {
    if (std::abs(x1 - x0) <= 1e-18) {
        return y0;
    }
    const double t = (x - x0) / (x1 - x0);
    return y0 * (1.0 - t) + y1 * t;
}

// 函数说明：在边界钳位条件下做线性插值。
double PhononMaterial::interp_linear_clamped(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
    if (xs.empty() || ys.empty()) {
        return 0.0;
    }
    if (xs.size() == 1 || ys.size() == 1) {
        return ys.front();
    }
    if (x <= xs.front()) {
        return ys.front();
    }
    if (x >= xs.back()) {
        return ys.back();
    }
    auto it = std::upper_bound(xs.begin(), xs.end(), x);
    const size_t hi = static_cast<size_t>(std::distance(xs.begin(), it));
    const size_t lo = hi - 1;
    return lerp(xs[lo], xs[hi], ys[lo], ys[hi], x);
}

// 函数说明：从 POSCAR 读取晶格并计算原胞体积与倒易基矢。
bool PhononMaterial::load_poscar_lattice_volume(const std::string& folder, std::string* err) {
    try {
        namespace fs = std::filesystem;
        const fs::path poscar = fs::path(folder) / "POSCAR";
        if (!fs::exists(poscar)) {
            throw std::runtime_error("POSCAR not found: " + poscar.string());
        }
        std::ifstream in(poscar);
        if (!in) {
            throw std::runtime_error("Failed to open POSCAR: " + poscar.string());
        }

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line)) {
            const std::string t = trim(line);
            if (!t.empty()) {
                lines.push_back(t);
            }
        }
        if (lines.size() < 5) {
            throw std::runtime_error("POSCAR has insufficient lines.");
        }
        const double scale = std::stod(lines[1]);
        std::array<std::array<double, 3>, 3> a {};
        for (int i = 0; i < 3; ++i) {
            std::istringstream iss(lines[2 + i]);
            if (!(iss >> a[i][0] >> a[i][1] >> a[i][2])) {
                throw std::runtime_error("Failed to parse POSCAR lattice vector.");
            }
            a[i][0] *= scale;
            a[i][1] *= scale;
            a[i][2] *= scale;
        }

        const double det =
            a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
            a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
            a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
        unit_cell_volume_ = std::abs(det);
        if (!(unit_cell_volume_ > 0.0)) {
            throw std::runtime_error("POSCAR unit-cell volume is non-positive.");
        }
        Mat3 inv {};
        if (!inverse3x3(a, inv)) {
            throw std::runtime_error("Failed to invert POSCAR lattice.");
        }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                reciprocal_lattice_[i][j] = inv[i][j] * (2.0 * M_PI);
            }
        }
        return true;
    } catch (const std::exception& ex) {
        if (err != nullptr) {
            *err = ex.what();
        }
        unit_cell_volume_ = 1.0;
        return false;
    }
}

// 函数说明：从 HDF5 加载频率、群速度、散射率等声子数据库。
bool PhononMaterial::load_hdf5_data(const SimulationConfig& args, int mat_index, std::string* err) {
    try {
        (void) mat_index;
        namespace fs = std::filesystem;
        fs::path folder;
        if (!args.material_folder.empty()) {
            folder = fs::path(args.material_folder);
        } else {
            folder = fs::path("Material/Si");
        }
        if (!folder.is_absolute()) {
            fs::path p1 = fs::current_path() / folder;
            fs::path p2 = fs::path(args.input_directory) / folder;
            folder = fs::exists(p1) ? p1 : p2;
        }
        folder = folder.lexically_normal();
        material_folder_path_ = folder.string();
        if (!material_folder_path_.empty()) {
            std::cout << "PhononMaterial material folder: " << material_folder_path_ << '\n';
        }
        std::string poscar_err;
        if (!load_poscar_lattice_volume(material_folder_path_, &poscar_err)) {
            std::cerr << "Warning: " << poscar_err << ". Using unit_cell_volume=1.0 A^3.\n";
        }

        fs::path hdf_path = folder / "kappa-fbz.hdf5";
        if (!fs::exists(hdf_path)) {
            hdf_path = folder / "kappa.hdf5";
        }
        if (!fs::exists(hdf_path)) {
            // Legacy material packs may use different file names (e.g. kappa-m313131.hdf5).
            for (const auto& entry : fs::directory_iterator(folder)) {
                if (!entry.is_regular_file()) {
                    continue;
                }
                if (entry.path().extension() == ".hdf5") {
                    hdf_path = entry.path();
                    break;
                }
            }
        }
        if (!fs::exists(hdf_path)) {
            throw std::runtime_error("No HDF5 file found under material folder: " + folder.string());
        }
        std::cout << "PhononMaterial HDF5 file: " << hdf_path.string() << '\n';

        H5::H5File file(hdf_path.string(), H5F_ACC_RDONLY);

        std::vector<hsize_t> d_omega;
        std::vector<double> omega_thz = read_dataset_nd_double(file, "frequency", d_omega);
        if (d_omega.size() != 2) {
            throw std::runtime_error("frequency rank mismatch");
        }
        qpoint_count_ = static_cast<int>(d_omega[0]);
        branch_count_ = static_cast<int>(d_omega[1]);
        const int nm = qpoint_count_ * branch_count_;
        mode_angular_frequency_data_.resize(static_cast<size_t>(nm));
        mode_wavevector_norm_data_.assign(static_cast<size_t>(nm), 0.0);
        mode_frequency_window_data_.assign(static_cast<size_t>(nm), 0.0);
        for (int i = 0; i < nm; ++i) {
            mode_angular_frequency_data_[static_cast<size_t>(i)] = std::max(0.0, omega_thz[static_cast<size_t>(i)]) * 2.0 * M_PI;
        }

        std::vector<hsize_t> d_qp;
        std::vector<double> qp = read_dataset_nd_double(file, "qpoint", d_qp);
        if (d_qp.size() != 2 || static_cast<int>(d_qp[0]) != qpoint_count_ || d_qp[1] != 3) {
            throw std::runtime_error("qpoint shape mismatch");
        }
        qpoint_fractions_.assign(static_cast<size_t>(qpoint_count_), {0.0, 0.0, 0.0});
        for (int q = 0; q < qpoint_count_; ++q) {
            const size_t b = static_cast<size_t>(q) * 3u;
            qpoint_fractions_[static_cast<size_t>(q)] = {qp[b + 0], qp[b + 1], qp[b + 2]};
        }

        std::vector<hsize_t> d_mesh;
        std::vector<double> mesh_data = read_dataset_nd_double(file, "mesh", d_mesh);
        if (d_mesh.size() == 1 && d_mesh[0] >= 3) {
            mesh_grid_size_[0] = std::max(1, static_cast<int>(std::llround(mesh_data[0])));
            mesh_grid_size_[1] = std::max(1, static_cast<int>(std::llround(mesh_data[1])));
            mesh_grid_size_[2] = std::max(1, static_cast<int>(std::llround(mesh_data[2])));
        } else {
            mesh_grid_size_ = {1, 1, 1};
        }

        std::vector<hsize_t> d_gv;
        std::vector<double> gv = read_dataset_nd_double(file, "group_velocity", d_gv);
        if (d_gv.size() != 3 || static_cast<int>(d_gv[0]) != qpoint_count_ ||
            static_cast<int>(d_gv[1]) != branch_count_ || d_gv[2] != 3) {
            throw std::runtime_error("group_velocity shape mismatch");
        }
        mode_group_velocity_data_.resize(static_cast<size_t>(nm));
        for (int q = 0; q < qpoint_count_; ++q) {
            for (int b = 0; b < branch_count_; ++b) {
                const size_t base = (static_cast<size_t>(q) * static_cast<size_t>(branch_count_) + static_cast<size_t>(b)) * 3u;
                mode_group_velocity_data_[static_cast<size_t>(q * branch_count_ + b)] = {gv[base + 0], gv[base + 1], gv[base + 2]};
            }
        }
        const Vec3 mesh_q {
            1.0 / (2.0 * static_cast<double>(mesh_grid_size_[0])),
            1.0 / (2.0 * static_cast<double>(mesh_grid_size_[1])),
            1.0 / (2.0 * static_cast<double>(mesh_grid_size_[2]))
        };
        const Vec3 k_grid = q_to_k(reciprocal_lattice_, mesh_q);
        for (int q = 0; q < qpoint_count_; ++q) {
            const Vec3 kvec = q_to_k(reciprocal_lattice_, qpoint_fractions_[static_cast<size_t>(q)]);
            const double kn = vec_norm(kvec);
            for (int b = 0; b < branch_count_; ++b) {
                const int fi = q * branch_count_ + b;
                const Vec3 v = mode_group_velocity_data_[static_cast<size_t>(fi)];
                mode_wavevector_norm_data_[static_cast<size_t>(fi)] = kn;
                mode_frequency_window_data_[static_cast<size_t>(fi)] = std::sqrt(
                    (v[0] * k_grid[0]) * (v[0] * k_grid[0]) +
                    (v[1] * k_grid[1]) * (v[1] * k_grid[1]) +
                    (v[2] * k_grid[2]) * (v[2] * k_grid[2]));
            }
        }

        std::vector<hsize_t> d_temp;
        temperature_samples_ = read_dataset_nd_double(file, "temperature", d_temp);
        if (d_temp.size() != 1) {
            throw std::runtime_error("temperature rank mismatch");
        }

        std::vector<hsize_t> d_gamma;
        gamma_table_ = read_dataset_nd_double(file, "gamma", d_gamma);
        if (d_gamma.size() != 3 || static_cast<int>(d_gamma[1]) != qpoint_count_ ||
            static_cast<int>(d_gamma[2]) != branch_count_) {
            throw std::runtime_error("gamma shape mismatch");
        }

        active_mode_list_.clear();
        for (int q = 0; q < qpoint_count_; ++q) {
            for (int b = 0; b < branch_count_; ++b) {
                const auto& v = mode_group_velocity_data_[static_cast<size_t>(q * branch_count_ + b)];
                const double vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                if (vn > 1e-12) {
                    active_mode_list_.push_back({q, b});
                }
            }
        }
        active_mode_count_ = static_cast<int>(active_mode_list_.size());
        if (active_mode_count_ == 0) {
            throw std::runtime_error("No active modes found in HDF5 data.");
        }
        flat_to_active_index_.assign(static_cast<size_t>(nm), -1);
        active_to_flat_index_.assign(active_mode_list_.size(), -1);
        for (size_t ai = 0; ai < active_mode_list_.size(); ++ai) {
            const int fi = flatten_mode_index(active_mode_list_[ai]);
            if (fi >= 0 && fi < nm) {
                flat_to_active_index_[static_cast<size_t>(fi)] = static_cast<int>(ai);
                active_to_flat_index_[ai] = fi;
            }
        }
        build_degenerate_mode_map();
        return true;
    } catch (const std::exception& ex) {
        if (err != nullptr) {
            *err = ex.what();
        }
        return false;
    }
}

// 函数说明：构建测试用合成模态库（非物理生产用途）。
void PhononMaterial::build_fallback_modes(const SimulationConfig& args) {
    if (!args.material_folder.empty()) {
        std::cout << "PhononMaterial material folder: " << args.material_folder << '\n';
    }
    const int n_q = 64;
    qpoint_count_ = n_q;
    branch_count_ = 3;
    active_mode_count_ = n_q * branch_count_;
    unit_cell_volume_ = 1.0;
    mode_angular_frequency_data_.resize(static_cast<size_t>(active_mode_count_));
    mode_wavevector_norm_data_.assign(static_cast<size_t>(active_mode_count_), 1.0);
    mode_frequency_window_data_.assign(static_cast<size_t>(active_mode_count_), 1e-3);
    mode_group_velocity_data_.resize(static_cast<size_t>(active_mode_count_));
    temperature_samples_ = {100.0, 200.0, 300.0, 400.0, 500.0};
    gamma_table_.assign(temperature_samples_.size() * static_cast<size_t>(active_mode_count_), 0.05);
    active_mode_list_.clear();
    active_mode_list_.reserve(static_cast<size_t>(active_mode_count_));

    std::mt19937_64 rng(20250209);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int q = 0; q < n_q; ++q) {
        for (int b = 0; b < branch_count_; ++b) {
            const int m = q * branch_count_ + b;
            const double frac = (static_cast<double>(q) + 0.5) / static_cast<double>(n_q);
            mode_angular_frequency_data_[static_cast<size_t>(m)] = (2.0 + 18.0 * frac) * (1.0 + 0.15 * b) * 2.0 * M_PI;
            const double vmag = 20.0 + 90.0 * U(rng);
            const auto dir = random_unit_vector(rng);
            mode_group_velocity_data_[static_cast<size_t>(m)] = {dir[0] * vmag, dir[1] * vmag, dir[2] * vmag};
            active_mode_list_.push_back({q, b});
        }
    }
    flat_to_active_index_.assign(static_cast<size_t>(active_mode_count_), -1);
    active_to_flat_index_.assign(active_mode_list_.size(), -1);
    for (size_t ai = 0; ai < active_mode_list_.size(); ++ai) {
        const int fi = flatten_mode_index(active_mode_list_[ai]);
        if (fi >= 0 && fi < active_mode_count_) {
            flat_to_active_index_[static_cast<size_t>(fi)] = static_cast<int>(ai);
            active_to_flat_index_[ai] = fi;
        }
    }
    build_degenerate_mode_map();
}

// 函数说明：采样随机单位方向用于合成模式或方向扰动。
PhononMaterial::Vec3 PhononMaterial::random_unit_vector(std::mt19937_64& rng) {
    std::normal_distribution<double> N(0.0, 1.0);
    Vec3 v {N(rng), N(rng), N(rng)};
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n <= 1e-14) {
        return {1.0, 0.0, 0.0};
    }
    return {v[0] / n, v[1] / n, v[2] / n};
}

// 函数说明：在活跃模态集合中采样一个声子模式。
PhononMaterial::Mode PhononMaterial::sample_active_mode(std::mt19937_64& rng) const {
    if (active_mode_list_.empty()) {
        throw std::runtime_error("No active modes in PhononMaterial.");
    }
    std::uniform_int_distribution<int> pick(0, static_cast<int>(active_mode_list_.size()) - 1);
    return active_mode_list_[static_cast<size_t>(pick(rng))];
}

// 函数说明：查询给定模态的群速度矢量。
PhononMaterial::Vec3 PhononMaterial::mode_group_velocity(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_group_velocity_data_.size())) {
        return {0.0, 0.0, 0.0};
    }
    return mode_group_velocity_data_[static_cast<size_t>(m)];
}

// 函数说明：查询给定模态的角频率。
double PhononMaterial::mode_angular_frequency(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_angular_frequency_data_.size())) {
        return 1.0;
    }
    return mode_angular_frequency_data_[static_cast<size_t>(m)];
}

// 函数说明：查询给定模态波矢模长。
double PhononMaterial::mode_wavevector_norm(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_wavevector_norm_data_.size())) {
        return 0.0;
    }
    return mode_wavevector_norm_data_[static_cast<size_t>(m)];
}

// 函数说明：查询给定模态对应频率窗宽。
double PhononMaterial::mode_frequency_window(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_frequency_window_data_.size())) {
        return 0.0;
    }
    return mode_frequency_window_data_[static_cast<size_t>(m)];
}

// 函数说明：预计算简并模态映射以支持粗糙边界散射匹配。
void PhononMaterial::build_degenerate_mode_map() {
    const int nm = qpoint_count_ * branch_count_;
    degenerate_partner_branch_data_.assign(static_cast<size_t>(std::max(0, nm)), -1);
    if (qpoint_count_ <= 0 || branch_count_ <= 1) {
        return;
    }
    constexpr double tol = 1e-10;
    for (int q = 0; q < qpoint_count_; ++q) {
        for (int j1 = 0; j1 < branch_count_; ++j1) {
            const int f1 = q * branch_count_ + j1;
            for (int j2 = j1 + 1; j2 < branch_count_; ++j2) {
                const int f2 = q * branch_count_ + j2;
                if (std::abs(mode_angular_frequency_data_[static_cast<size_t>(f1)] - mode_angular_frequency_data_[static_cast<size_t>(f2)]) < tol) {
                    degenerate_partner_branch_data_[static_cast<size_t>(f1)] = j2;
                    degenerate_partner_branch_data_[static_cast<size_t>(f2)] = j1;
                    break;
                }
            }
        }
    }
}

// 函数说明：计算给定温度下模态的玻色占据数。
double PhononMaterial::bose_occupation(double temperature, const Mode& mode) const {
    const double w = mode_angular_frequency(mode);
    if (temperature <= 0.0 || w <= 0.0) {
        return 0.0;
    }
    const double x = hbar_ * w / (kb_ * temperature);
    if (x > 700.0) {
        return 0.0;
    }
    const double ex = std::exp(x);
    return 1.0 / std::max(ex - 1.0, 1e-12);
}

// 函数说明：计算模态能量（含占据项与零点项）。
double PhononMaterial::mode_energy(double temperature, const Mode& mode) const {
    return hbar_ * mode_angular_frequency(mode) * (bose_occupation(temperature, mode) + 0.5);
}

// 函数说明：按温度和散射率表计算模态寿命。
double PhononMaterial::mode_lifetime(double temperature, const Mode& mode) const {
    if (gamma_table_.empty() || temperature_samples_.empty()) {
        const double T = std::max(1.0, temperature);
        const double w = std::max(1e-9, mode_angular_frequency(mode));
        const double tau = 12.0 / (1.0 + 0.05 * w + 0.002 * (T - 300.0) * (T - 300.0) / 100.0);
        return std::max(0.1, tau);
    }
    const int m = flatten_mode_index(mode);
    if (m < 0) {
        return 0.0;
    }
    const size_t nmode = static_cast<size_t>(qpoint_count_ * branch_count_);
    if (nmode == 0 || gamma_table_.empty()) {
        return 0.0;
    }
    if (temperature <= temperature_samples_.front()) {
        const size_t idx = static_cast<size_t>(m);
        if (idx >= gamma_table_.size()) {
            return 0.0;
        }
        return tau_from_gamma(gamma_table_[idx]);
    }
    if (temperature >= temperature_samples_.back()) {
        const size_t idx = (temperature_samples_.size() - 1) * nmode + static_cast<size_t>(m);
        if (idx >= gamma_table_.size()) {
            return 0.0;
        }
        return tau_from_gamma(gamma_table_[idx]);
    }
    auto it = std::upper_bound(temperature_samples_.begin(), temperature_samples_.end(), temperature);
    const size_t hi = static_cast<size_t>(std::distance(temperature_samples_.begin(), it));
    const size_t lo = hi - 1;
    const size_t idx0 = lo * nmode + static_cast<size_t>(m);
    const size_t idx1 = hi * nmode + static_cast<size_t>(m);
    if (idx1 >= gamma_table_.size()) {
        return 0.0;
    }
    const double tau0 = tau_from_gamma(gamma_table_[idx0]);
    const double tau1 = tau_from_gamma(gamma_table_[idx1]);
    return lerp(temperature_samples_[lo], temperature_samples_[hi], tau0, tau1, temperature);
}

// 函数说明：将模态总量归一化为晶体能量密度量纲。
double PhononMaterial::normalize_to_energy_density(double x) const {
    const double den = std::max(1.0, static_cast<double>(qpoint_count_) * unit_cell_volume_);
    return x / den;
}

// 函数说明：将模态总量归一化为晶体能量密度量纲。
PhononMaterial::Vec3 PhononMaterial::normalize_to_energy_density(const Vec3& x) const {
    const double den = std::max(1.0, static_cast<double>(qpoint_count_) * unit_cell_volume_);
    return {x[0] / den, x[1] / den, x[2] / den};
}

// 函数说明：累计所有活跃模态得到晶体能量密度。
double PhononMaterial::crystal_energy_density(double temperature) const {
    double e = 0.0;
    for (const Mode& mode : active_mode_list_) {
        const double w = mode_angular_frequency(mode);
        if (w <= 0.0) {
            continue;
        }
        const double n = bose_occupation(temperature, mode);
        e += hbar_ * w * n;
    }
    return normalize_to_energy_density(e) + zero_point_energy_density_;
}

// 函数说明：由温度正向查表得到能量密度（插值），缺表时回退全模态积分。
double PhononMaterial::energy_density_from_temperature(double temperature) const {
    if (temperature_lookup_table_.empty() || energy_lookup_table_.empty()) {
        return crystal_energy_density(temperature);
    }
    return interp_linear_clamped(temperature_lookup_table_, energy_lookup_table_, temperature);
}

// 函数说明：由能量密度反查温度（查找表插值）。
double PhononMaterial::temperature_from_energy_density(double energy_density) const {
    return interp_linear_clamped(energy_lookup_table_, temperature_lookup_table_, energy_density);
}

// 函数说明：构建温度-能量单调查找表用于快速反演。
void PhononMaterial::initialize_temperature_lookup() {
    zero_point_energy_density_ = 0.0;
    for (double w : mode_angular_frequency_data_) {
        zero_point_energy_density_ += 0.5 * hbar_ * std::max(0.0, w);
    }
    zero_point_energy_density_ = normalize_to_energy_density(zero_point_energy_density_);

    const double tmin = temperature_samples_.empty() ? 0.0 : *std::min_element(temperature_samples_.begin(), temperature_samples_.end());
    const double tmax = temperature_samples_.empty() ? 1000.0 : *std::max_element(temperature_samples_.begin(), temperature_samples_.end());
    const double dT = std::max(1e-6, temperature_lookup_dt_);
    const int n = std::max(2, static_cast<int>(std::floor((tmax - tmin) / dT + 0.5)) + 1);

    temperature_lookup_table_.assign(static_cast<size_t>(n), 0.0);
    energy_lookup_table_.assign(static_cast<size_t>(n), 0.0);
    std::vector<double> raw_energy(static_cast<size_t>(n), 0.0);

#ifdef EPMC_USE_OPENMP
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i) {
        const double T = tmin + static_cast<double>(i) * dT;
        raw_energy[static_cast<size_t>(i)] = crystal_energy_density(T);
    }
#else
    const unsigned int workers = std::max(1U, std::thread::hardware_concurrency());
    if (workers <= 1 || n < 256) {
        for (int i = 0; i < n; ++i) {
            const double T = tmin + static_cast<double>(i) * dT;
            raw_energy[static_cast<size_t>(i)] = crystal_energy_density(T);
        }
    } else {
        std::vector<std::thread> pool;
        pool.reserve(workers);
        for (unsigned int w = 0; w < workers; ++w) {
            const int begin = static_cast<int>((static_cast<long long>(w) * n) / workers);
            const int end = static_cast<int>((static_cast<long long>(w + 1U) * n) / workers);
            pool.emplace_back([&, begin, end]() {
                for (int i = begin; i < end; ++i) {
                    const double T = tmin + static_cast<double>(i) * dT;
                    raw_energy[static_cast<size_t>(i)] = crystal_energy_density(T);
                }
            });
        }
        for (auto& th : pool) {
            if (th.joinable()) {
                th.join();
            }
        }
    }
#endif

    for (int i = 0; i < n; ++i) {
        const double T = tmin + static_cast<double>(i) * dT;
        temperature_lookup_table_[static_cast<size_t>(i)] = T;
        double E = raw_energy[static_cast<size_t>(i)];
        if (i > 0) {
            E = std::max(E, std::nextafter(energy_lookup_table_[static_cast<size_t>(i - 1)], std::numeric_limits<double>::infinity()));
        }
        energy_lookup_table_[static_cast<size_t>(i)] = E;
    }
}
