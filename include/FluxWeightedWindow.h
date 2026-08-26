#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace phonomc_detail {

// Remove the equilibrium outgoing flux already supplied by the discrete
// specular map. The remaining nonnegative rates define the diffuse channel.
// Mode indices in specular_match refer to the common active-mode indexing of
// the incoming and outgoing arrays.
inline std::vector<double> residual_diffuse_creation_rates(
    const std::vector<double>& outgoing_flux,
    const std::vector<double>& incoming_flux,
    const std::vector<double>& specularity,
    const std::vector<int>& specular_match) {
    std::vector<double> residual = outgoing_flux;
    const size_t count = std::min(
        incoming_flux.size(),
        std::min(specularity.size(), specular_match.size()));
    for (size_t incoming = 0; incoming < count; ++incoming) {
        const int outgoing = specular_match[incoming];
        if (outgoing < 0 || outgoing >= static_cast<int>(residual.size())) {
            continue;
        }
        const double flux = incoming_flux[incoming];
        const double probability = specularity[incoming];
        if (!std::isfinite(flux) || !std::isfinite(probability) ||
            flux <= 0.0 || probability <= 0.0) {
            continue;
        }
        residual[static_cast<size_t>(outgoing)] -=
            flux * std::clamp(probability, 0.0, 1.0);
    }
    for (double& rate : residual) {
        if (!std::isfinite(rate) || rate < 0.0) {
            rate = 0.0;
        }
    }
    return residual;
}

// Build a prefix sum in an explicitly supplied mode order. Negative and
// non-finite rates are treated as zero because they cannot define a sampling
// probability. Keeping the ordering explicit lets the caller use contiguous
// frequency windows without duplicating the rate array.
inline std::vector<double> ordered_nonnegative_prefix(
    const std::vector<double>& weights,
    const std::vector<int>& order) {
    std::vector<double> prefix(order.size() + 1u, 0.0);
    for (size_t i = 0; i < order.size(); ++i) {
        const int index = order[i];
        double weight = (index >= 0 && index < static_cast<int>(weights.size()))
            ? weights[static_cast<size_t>(index)] : 0.0;
        if (!std::isfinite(weight) || weight < 0.0) {
            weight = 0.0;
        }
        prefix[i + 1u] = prefix[i] + weight;
    }
    return prefix;
}

// Select one interval [prefix[i], prefix[i + 1]) inside [begin, end),
// using a unit-interval sample. Returns -1 for an invalid or zero-weight
// window. Keeping this deterministic core separate makes the rough-boundary
// probability law directly testable without exposing solver internals.
inline int flux_weighted_window_index(
    const std::vector<double>& prefix,
    int begin,
    int end,
    double unit_sample) {
    if (begin < 0 || end <= begin ||
        end >= static_cast<int>(prefix.size())) {
        return -1;
    }
    const double flux_begin = prefix[static_cast<size_t>(begin)];
    const double flux_end = prefix[static_cast<size_t>(end)];
    const double window_flux = flux_end - flux_begin;
    if (!(window_flux > 0.0) || !std::isfinite(window_flux)) {
        return -1;
    }
    const double u = std::clamp(
        unit_sample,
        0.0,
        std::nextafter(1.0, 0.0));
    const double target = flux_begin + u * window_flux;
    const auto first = prefix.begin() + begin + 1;
    const auto last = prefix.begin() + end + 1;
    const auto it = std::upper_bound(first, last, target);
    return std::clamp(
        static_cast<int>(std::distance(prefix.begin(), it)) - 1,
        begin,
        end - 1);
}

}  // namespace phonomc_detail
