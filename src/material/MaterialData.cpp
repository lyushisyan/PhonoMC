#include "material/MaterialData.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace phonomc {
namespace {
void require_finite(const std::vector<double>& values, const char* name) {
    for (double x : values)
        if (!std::isfinite(x)) throw std::runtime_error(std::string(name) + " contains non-finite values");
}
void require_finite(const std::vector<MaterialVec3>& values, const char* name) {
    for (const auto& v : values)
        for (double x : v)
            if (!std::isfinite(x)) throw std::runtime_error(std::string(name) + " contains non-finite values");
}
} // namespace

void validate_material_data(const MaterialData& data) {
    if (!std::isfinite(data.lattice.volume_a3) || data.lattice.volume_a3 <= 0)
        throw std::runtime_error("Unit-cell volume must be finite and positive");
    for (const auto& row : data.lattice.reciprocal)
        for (double x : row)
            if (!std::isfinite(x)) throw std::runtime_error("Non-finite reciprocal lattice");
    if (data.qpoints <= 0 || data.branches <= 0 ||
        data.qpoints > std::numeric_limits<int>::max() / data.branches)
        throw std::runtime_error("Material mode dimensions must be positive and fit integer indices");
    const size_t nm = static_cast<size_t>(data.qpoints) * static_cast<size_t>(data.branches);
    if (data.frequency_thz.size() != nm) throw std::runtime_error("frequency shape mismatch");
    if (data.qpoint_fractions.size() != static_cast<size_t>(data.qpoints))
        throw std::runtime_error("qpoint shape mismatch");
    if (data.group_velocity.size() != nm) throw std::runtime_error("group_velocity shape mismatch");
    size_t full_grid_size = 1;
    for (int count : data.mesh) {
        if (count <= 0) throw std::runtime_error("mesh entries must be positive integers");
        if (static_cast<size_t>(count) > std::numeric_limits<size_t>::max() / full_grid_size)
            throw std::runtime_error("mesh size overflow");
        full_grid_size *= static_cast<size_t>(count);
    }
    if (full_grid_size != static_cast<size_t>(data.qpoints))
        throw std::runtime_error("Full Brillouin-zone grid required: qpoint count differs from mesh product. Expand reduced data first.");
    require_finite(data.frequency_thz, "frequency");
    require_finite(data.qpoint_fractions, "qpoint");
    require_finite(data.group_velocity, "group_velocity");
    require_finite(data.temperatures, "temperature");
    if (data.temperatures.empty() || data.temperatures.front() < 0 ||
        !std::is_sorted(data.temperatures.begin(), data.temperatures.end()) ||
        std::adjacent_find(data.temperatures.begin(), data.temperatures.end()) != data.temperatures.end())
        throw std::runtime_error("temperature must be nonempty, nonnegative and strictly increasing");
    if (data.temperatures.size() > std::numeric_limits<size_t>::max() / nm ||
        data.gamma.size() != data.temperatures.size() * nm)
        throw std::runtime_error("gamma shape mismatch");
    require_finite(data.gamma, "gamma");
    if (std::any_of(data.gamma.begin(), data.gamma.end(), [](double g) { return g < 0; }))
        throw std::runtime_error("gamma must be nonnegative");
    if (data.gamma_normal.empty() != data.gamma_umklapp.empty())
        throw std::runtime_error("gamma_N and gamma_U must both be supplied or both omitted");
    for (const auto* values : {&data.gamma_normal, &data.gamma_umklapp}) {
        if (values->empty()) continue;
        if (values->size()!=data.gamma.size()) throw std::runtime_error("N/U linewidth shape mismatch");
        require_finite(*values, "N/U linewidth");
        if (std::any_of(values->begin(),values->end(),[](double g){ return g<0; }))
            throw std::runtime_error("N/U linewidths must be nonnegative");
    }
    for (const auto& [name, values] : data.gamma_resistive) {
        if (name != "gamma_isotope" && name != "gamma_impurity" && name != "gamma_defect")
            throw std::runtime_error("Unsupported resistive linewidth: " + name);
        if (values.size() != nm && values.size() != data.gamma.size())
            throw std::runtime_error(name + " shape mismatch");
        require_finite(values, name.c_str());
        if (std::any_of(values.begin(), values.end(), [](double g) { return g < 0; }))
            throw std::runtime_error(name + " must be nonnegative");
    }
    for (size_t i = 0; i < data.gamma.size(); ++i) {
        double total = data.gamma[i];
        for (const auto& [name, values] : data.gamma_resistive)
            total += values[values.size() == nm ? i % nm : i];
        if (!std::isfinite(total))
            throw std::runtime_error("Total linewidth overflow including resistive scattering");
    }
}
} // namespace phonomc
