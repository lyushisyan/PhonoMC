#include "material/MaterialDataReader.h"
#include <H5Cpp.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace phonomc {
namespace {
using Mat3 = std::array<std::array<double, 3>, 3>;
using Vec3 = std::array<double, 3>;

// 函数说明：读取 HDF5 数值数据集并展开为一维数组。
std::vector<double> read_dataset_nd_double(const H5::H5File& file, const std::string& name, std::vector<hsize_t>& dims_out) {
    if (H5Lexists(file.getId(), name.c_str(), H5P_DEFAULT) <= 0) {
        throw std::runtime_error("Missing HDF5 dataset: " + name);
    }
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
        if (d == 0 || d > std::numeric_limits<size_t>::max() / n) {
            throw std::runtime_error("Dataset " + name + " is empty or too large");
        }
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

} // namespace

LatticeData read_poscar_lattice(const std::filesystem::path& path) {
    LatticeData lattice;
    namespace fs = std::filesystem;
    const auto& poscar = path;
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
    if (!std::isfinite(scale)) throw std::runtime_error("Non-finite POSCAR scale.");
    std::array<std::array<double, 3>, 3> a {};
    for (int i = 0; i < 3; ++i) {
        std::istringstream iss(lines[2 + i]);
        if (!(iss >> a[i][0] >> a[i][1] >> a[i][2])) {
            throw std::runtime_error("Failed to parse POSCAR lattice vector.");
        }
        a[i][0] *= scale;
        a[i][1] *= scale;
        a[i][2] *= scale;
        for (double x : a[i]) if (!std::isfinite(x)) throw std::runtime_error("Non-finite POSCAR lattice.");
    }

    const double det =
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
        a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
        a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    lattice.volume_a3 = std::abs(det);
    if (!std::isfinite(lattice.volume_a3) || !(lattice.volume_a3 > 0.0)) {
        throw std::runtime_error("POSCAR unit-cell volume is non-positive.");
    }
    Mat3 inv {};
    if (!inverse3x3(a, inv)) {
        throw std::runtime_error("Failed to invert POSCAR lattice.");
    }
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            lattice.reciprocal[i][j] = inv[i][j] * (2.0 * M_PI);
        }
    }
    for (const auto& row : lattice.reciprocal)
        for (double x : row) if (!std::isfinite(x)) throw std::runtime_error("Non-finite reciprocal lattice.");
    return lattice;
}

MaterialData read_material_data(const std::filesystem::path& folder, std::ostream* progress) {
    try {
        namespace fs = std::filesystem;
        MaterialData data;
        try { data.lattice = read_poscar_lattice(folder / "POSCAR"); }
        catch (const std::exception& e) {
            throw std::runtime_error(std::string("Strict POSCAR loading failed: ") + e.what());
        }
        fs::path hdf_path = folder / "kappa-fbz.hdf5";
        if (!fs::exists(hdf_path)) {
            hdf_path = folder / "kappa.hdf5";
        }
        if (!fs::exists(hdf_path)) {
            // Legacy material packs may use different file names (e.g. kappa-m313131.hdf5).
            std::vector<fs::path> candidates;
            for (const auto& entry : fs::directory_iterator(folder)) {
                if (!entry.is_regular_file()) {
                    continue;
                }
                if (entry.path().extension() == ".hdf5") {
                    candidates.push_back(entry.path());
                }
            }
            if (candidates.size() > 1) {
                throw std::runtime_error("Ambiguous HDF5 files: provide kappa-fbz.hdf5 or kappa.hdf5.");
            }
            if (candidates.size() == 1) hdf_path = candidates.front();
        }
        if (!fs::exists(hdf_path)) {
            throw std::runtime_error("No HDF5 file found under material folder: " + folder.string());
        }
        if (progress) *progress << "PhononMaterial HDF5 file: " << hdf_path.string() << '\n';

        H5::H5File file(hdf_path.string(), H5F_ACC_RDONLY);

        std::vector<hsize_t> d_omega;
        std::vector<double> omega_thz = read_dataset_nd_double(file, "frequency", d_omega);
        if (d_omega.size() != 2) {
            throw std::runtime_error("frequency rank mismatch");
        }
        if (d_omega[0] > static_cast<hsize_t>(std::numeric_limits<int>::max()) / d_omega[1]) {
            throw std::runtime_error("frequency has too many modes for integer mode indices");
        }
        data.qpoints = static_cast<int>(d_omega[0]);
        data.branches = static_cast<int>(d_omega[1]);
        const int nm = data.qpoints * data.branches;
        data.frequency_thz = std::move(omega_thz);

        std::vector<hsize_t> d_qp;
        std::vector<double> qp = read_dataset_nd_double(file, "qpoint", d_qp);
        if (d_qp.size() != 2 || d_qp[0] != static_cast<hsize_t>(data.qpoints) || d_qp[1] != 3) {
            throw std::runtime_error("qpoint shape mismatch");
        }
        data.qpoint_fractions.assign(static_cast<size_t>(data.qpoints), {0.0, 0.0, 0.0});
        for (int q = 0; q < data.qpoints; ++q) {
            const size_t b = static_cast<size_t>(q) * 3u;
            data.qpoint_fractions[static_cast<size_t>(q)] = {qp[b + 0], qp[b + 1], qp[b + 2]};
        }

        std::vector<hsize_t> d_mesh;
        std::vector<double> mesh_data = read_dataset_nd_double(file, "mesh", d_mesh);
        if (d_mesh.size() != 1 || d_mesh[0] != 3) {
            throw std::runtime_error("mesh must have shape (3)");
        }
        for (int axis = 0; axis < 3; ++axis) {
            const double count = mesh_data[static_cast<size_t>(axis)];
            if (!std::isfinite(count) || count < 1 || count != std::floor(count) ||
                count > std::numeric_limits<int>::max()) {
                throw std::runtime_error("mesh entries must be positive integers");
            }
            data.mesh[axis] = static_cast<int>(count);
        }
        if (H5Lexists(file.getId(), "weight", H5P_DEFAULT) > 0) {
            std::vector<hsize_t> d_weight;
            const auto weights = read_dataset_nd_double(file, "weight", d_weight);
            if (d_weight.size() != 1 || d_weight[0] != d_omega[0] ||
                std::any_of(weights.begin(), weights.end(), [](double w) { return w != 1.0; })) {
                throw std::runtime_error("Full Brillouin-zone data requires one unit weight per qpoint");
            }
        }

        std::vector<hsize_t> d_gv;
        std::vector<double> gv = read_dataset_nd_double(file, "group_velocity", d_gv);
        if (d_gv.size() != 3 || d_gv[0] != static_cast<hsize_t>(data.qpoints) ||
            d_gv[1] != static_cast<hsize_t>(data.branches) || d_gv[2] != 3) {
            throw std::runtime_error("group_velocity shape mismatch");
        }
        data.group_velocity.resize(static_cast<size_t>(nm));
        for (int q = 0; q < data.qpoints; ++q) {
            for (int b = 0; b < data.branches; ++b) {
                const size_t base = (static_cast<size_t>(q) * static_cast<size_t>(data.branches) + static_cast<size_t>(b)) * 3u;
                data.group_velocity[static_cast<size_t>(q * data.branches + b)] = {gv[base + 0], gv[base + 1], gv[base + 2]};
            }
        }
        std::vector<hsize_t> d_temp;
        data.temperatures = read_dataset_nd_double(file, "temperature", d_temp);
        if (d_temp.size() != 1) {
            throw std::runtime_error("temperature rank mismatch");
        }

        std::vector<hsize_t> d_gamma;
        data.gamma = read_dataset_nd_double(file, "gamma", d_gamma);
        if (d_gamma.size() != 3 || d_gamma[0] != data.temperatures.size() || d_gamma[1] != static_cast<hsize_t>(data.qpoints) ||
            d_gamma[2] != static_cast<hsize_t>(data.branches)) {
            throw std::runtime_error("gamma shape mismatch");
        }
        for (const auto& entry : {std::make_pair("gamma_N", &data.gamma_normal),
                                  std::make_pair("gamma_U", &data.gamma_umklapp)}) {
            if (H5Lexists(file.getId(), entry.first, H5P_DEFAULT)<=0) continue;
            std::vector<hsize_t> dims;
            *entry.second=read_dataset_nd_double(file,entry.first,dims);
            if (dims!=d_gamma) throw std::runtime_error(std::string(entry.first)+" shape mismatch");
        }
        for (const char* name : {"gamma_isotope", "gamma_impurity", "gamma_defect"}) {
            if (H5Lexists(file.getId(), name, H5P_DEFAULT) <= 0) continue;
            std::vector<hsize_t> dims;
            auto values = read_dataset_nd_double(file, name, dims);
            if (dims != d_omega && dims != d_gamma)
                throw std::runtime_error(std::string(name) + " shape mismatch: expected (q, branch) or (temperature, q, branch)");
            data.gamma_resistive.emplace(name, std::move(values));
            if (progress) *progress << "Including resistive linewidth " << name << "\n";
        }
        validate_material_data(data);
        return data;
    } catch (const H5::Exception& ex) {
        throw std::runtime_error("HDF5: " + ex.getDetailMsg());
    }
}
} // namespace phonomc
