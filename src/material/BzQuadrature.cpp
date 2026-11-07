#include "material/BzQuadrature.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace phonomc {
unsigned complete_callaway_bz(MaterialData& d) {
    validate_material_data(d);
    const auto& b = d.lattice.reciprocal;
    const double det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1]) -
                       b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0]) +
                       b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
    if (!std::isfinite(det) || det == 0) {
        throw std::invalid_argument("Callaway reciprocal lattice is singular.");
    }
    MaterialMat3 inv{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            inv[i][j] = (b[(j + 1) % 3][(i + 1) % 3] * b[(j + 2) % 3][(i + 2) % 3] -
                         b[(j + 1) % 3][(i + 2) % 3] * b[(j + 2) % 3][(i + 1) % 3]) /
                        det;
        }
    }
    const auto length2 = [&](const MaterialVec3& q) {
        double r = 0;
        for (int i = 0; i < 3; ++i) {
            double x = 0;
            for (int j = 0; j < 3; ++j) {
                x += b[i][j] * q[j];
            }
            r += x * x;
        }
        return r;
    };
    std::set<std::array<int, 3>> unique;
    std::vector<std::vector<MaterialVec3>> representatives;
    unsigned refinement = 1;
    for (const auto& q : d.qpoint_fractions) {
        std::array<int, 3> key{};
        for (int a = 0; a < 3; ++a) {
            const double x = (q[a] - std::floor(q[a])) * d.mesh[a];
            const long index = std::lround(x);
            if (std::abs(x - index) > 2e-6 * d.mesh[a]) {
                throw std::invalid_argument(
                    "Callaway requires a complete Gamma-centred full q mesh.");
            }
            key[a] = static_cast<int>(index % d.mesh[a]);
        }
        if (!unique.insert(key).second) {
            throw std::invalid_argument(
                "Callaway input q mesh contains duplicate points modulo reciprocal translations.");
        }
        const double original = length2(q);
        int lo[3], hi[3];
        double count = 1;
        for (int a = 0; a < 3; ++a) {
            double norm = 0;
            for (double x : inv[a]) {
                norm += x * x;
            }
            const double extent = std::sqrt(original * norm);
            if (!std::isfinite(extent) || extent > 100 || std::abs(q[a]) > 100) {
                throw std::invalid_argument(
                    "Callaway first-BZ validation exceeds supported lattice conditioning.");
            }
            lo[a] = static_cast<int>(std::ceil(q[a] - extent - 1e-10));
            hi[a] = static_cast<int>(std::floor(q[a] + extent + 1e-10));
            count *= hi[a] - lo[a] + 1;
        }
        if (count > 100000) {
            throw std::invalid_argument("Callaway first-BZ validation work limit exceeded.");
        }
        std::vector<std::pair<double, MaterialVec3>> candidates;
        double shortest = original;
        for (int i = lo[0]; i <= hi[0]; ++i) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int k = lo[2]; k <= hi[2]; ++k) {
                    MaterialVec3 p{q[0] - i, q[1] - j, q[2] - k};
                    const double r = length2(p);
                    shortest = std::min(shortest, r);
                    candidates.push_back({r, p});
                }
            }
        }
        const double tolerance = 1e-10 * std::max(shortest, 1e-20);
        if (original > shortest + tolerance) {
            throw std::invalid_argument("Callaway qpoint is outside the first Brillouin zone; "
                                        "supply first-BZ representatives.");
        }
        std::vector<MaterialVec3> reps;
        for (const auto& c : candidates) {
            if (std::abs(c.first - shortest) <= tolerance) {
                reps.push_back(c.second);
            }
        }
        if (reps.empty()) {
            throw std::runtime_error("Empty Callaway BZ representative set.");
        }
        const unsigned m = static_cast<unsigned>(reps.size());
        const auto next = static_cast<unsigned long long>(refinement / std::gcd(refinement, m)) * m;
        if (next > 96 || next * static_cast<unsigned long long>(d.qpoints) * d.branches >
                             std::numeric_limits<int>::max()) {
            throw std::length_error("Callaway equal-weight BZ refinement exceeds supported "
                                    "storage; use a smaller mesh.");
        }
        refinement = static_cast<unsigned>(next);
        representatives.push_back(std::move(reps));
    }
    if (refinement == 1) {
        return 1;
    }
    const size_t old_nm = d.frequency_thz.size();
    std::vector<size_t> source;
    source.reserve(old_nm * refinement);
    std::vector<MaterialVec3> qpoints;
    qpoints.reserve(d.qpoints * refinement);
    for (int q = 0; q < d.qpoints; ++q) {
        for (const auto& rep : representatives[q]) {
            for (unsigned copy = 0; copy < refinement / representatives[q].size(); ++copy) {
                qpoints.push_back(rep);
                for (int branch = 0; branch < d.branches; ++branch) {
                    source.push_back(static_cast<size_t>(q) * d.branches + branch);
                }
            }
        }
    }
    const auto expand = [&](const std::vector<double>& values) {
        if (values.empty()) {
            return values;
        }
        std::vector<double> out;
        out.reserve(values.size() * refinement);
        for (size_t t = 0; t < values.size() / old_nm; ++t) {
            for (size_t i : source) {
                out.push_back(values[t * old_nm + i]);
            }
        }
        return out;
    };
    std::vector<MaterialVec3> velocity;
    velocity.reserve(source.size());
    for (size_t i : source) {
        velocity.push_back(d.group_velocity[i]);
    }
    d.frequency_thz = expand(d.frequency_thz);
    d.gamma = expand(d.gamma);
    d.gamma_normal = expand(d.gamma_normal);
    d.gamma_umklapp = expand(d.gamma_umklapp);
    for (auto& entry : d.gamma_resistive) {
        entry.second = expand(entry.second);
    }
    d.qpoints *= refinement;
    d.qpoint_fractions = std::move(qpoints);
    d.group_velocity = std::move(velocity);
    return refinement;
}
} // namespace phonomc
