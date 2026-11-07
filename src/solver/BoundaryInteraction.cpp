#include "solver/BoundaryInteraction.h"

#include "PhononMaterial.h"

#include <algorithm>
#include <cmath>

namespace phonomc {
namespace {
using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 mul(const Vec3& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

Vec3 random_unit_vector(std::mt19937_64& rng) {
    // Keep the distribution local: caching a fourth Gaussian across events
    // would change the original RNG consumption.
    std::normal_distribution<double> N(0.0, 1.0);
    Vec3 v{N(rng), N(rng), N(rng)};
    const double n = norm(v);
    if (n <= 1e-12) {
        return {1.0, 0.0, 0.0};
    }
    return mul(v, 1.0 / n);
}

double
specularity(const BoundaryState& in, const BoundaryHit& hit, const PhononMaterial& material) {
    const double eta = std::max(0.0, hit.roughness);
    const double speed = std::max(1e-12, norm(in.velocity));
    const double incidence_cos =
        std::min(1.0, std::max(0.0, std::abs(dot(in.velocity, hit.outward_normal)) / speed));
    const double k_norm = std::max(1e-12, material.mode_angular_frequency(in.mode) / speed);
    const double x = 2.0 * eta * incidence_cos * k_norm;
    const double p = std::exp(-(x * x));
    return std::clamp(p, 0.0, 1.0);
}
} // namespace

BoundaryOutcome interact_boundary(const BoundaryState& incoming,
                                  const BoundaryHit& hit,
                                  const BoundaryPolicy& policy,
                                  const PhononMaterial& material,
                                  double energy_weight,
                                  const RoughBoundarySampler& sampler,
                                  std::mt19937_64& rng) {
    BoundaryOutcome result;
    result.state = incoming;
    auto fail = [&]() {
        result.action = BoundaryAction::Fail;
        result.error = BoundaryError::InvalidOutgoingMode;
    };
    if (!hit.valid) {
        if (policy.require_discrete_mode) {
            fail();
        }
        return result;
    }
    auto& out = result.state;
    const auto& n = hit.outward_normal;
    const double vn = dot(incoming.velocity, n);
    if (hit.condition == 'T') {
        result.action = BoundaryAction::Absorb;
        constexpr double kHbarEvPs = 6.582119569e-4;
        const double neq = material.bose_occupation(policy.reference_temperature, incoming.mode);
        const double raw_deviation_ev = kHbarEvPs * material.mode_angular_frequency(incoming.mode) *
                                        (incoming.occupation - neq);
        result.absorbed_energy_ev = raw_deviation_ev * energy_weight;
        return result;
    }
    if (hit.condition == 'P' && hit.has_periodic_pair) {
        result.action = BoundaryAction::Periodic;
        result.position_shift = hit.periodic_shift;
        return result;
    }
    result.action = BoundaryAction::Reflect;
    if (policy.require_discrete_mode) {
        if (hit.rough_table < 0) {
            fail();
            return result;
        }
        const auto reflected = sampler.sample(material,
                                              hit.rough_table,
                                              incoming.mode,
                                              incoming.occupation,
                                              policy.reference_temperature,
                                              rng,
                                              true);
        result.rough = reflected.diagnostics;
        const Vec3 velocity = material.mode_group_velocity(reflected.mode);
        const double omega = material.mode_angular_frequency(reflected.mode);
        if (!std::isfinite(reflected.occupation) || !std::isfinite(omega) || !(omega > 0.0) ||
            !std::isfinite(norm(velocity)) || !(dot(velocity, n) < 0.0)) {
            fail();
            return result;
        }
        out.mode = reflected.mode;
        out.velocity = velocity;
        out.omega = omega;
        out.occupation = reflected.occupation;
        result.clear_cached_energy = true;
    } else if (hit.condition == 'R') {
        if (hit.rough_table >= 0) {
            const auto reflected = sampler.sample(material,
                                                  hit.rough_table,
                                                  incoming.mode,
                                                  incoming.occupation,
                                                  policy.reference_temperature,
                                                  rng,
                                                  false);
            result.rough = reflected.diagnostics;
            out.mode = reflected.mode;
            out.velocity = material.mode_group_velocity(out.mode);
            if (dot(out.velocity, n) > 0.0) {
                out.velocity = sub(out.velocity, mul(n, 2.0 * dot(out.velocity, n)));
            }
            out.omega = material.mode_angular_frequency(out.mode);
            out.occupation = reflected.occupation;
            result.clear_cached_energy = true;
        } else {
            const double speed = norm(incoming.velocity);
            if (!std::isfinite(speed) || speed <= 0.0 ||
                !std::isfinite(hit.roughness) || hit.roughness < 0.0) {
                fail();
                return result;
            }
            result.rough.events = 1;
            result.rough.missing_data = 1;
            const double p_spec = specularity(incoming, hit, material);
            std::uniform_real_distribution<double> U(0.0, 1.0);
            if (U(rng) <= p_spec) {
                result.rough.specular = 1;
                out.velocity = sub(out.velocity, mul(n, 2.0 * vn));
            } else {
                result.rough.diffuse = 1;
                result.rough.global_random = 1;
                Vec3 dir = random_unit_vector(rng);
                if (dot(dir, n) > 0.0) {
                    dir = mul(dir, -1.0);
                }
                out.velocity = mul(dir, speed);
                out.omega = material.mode_angular_frequency(out.mode);
            }
        }
    } else {
        out.velocity = sub(out.velocity, mul(n, 2.0 * vn));
    }
    return result;
}
} // namespace phonomc
