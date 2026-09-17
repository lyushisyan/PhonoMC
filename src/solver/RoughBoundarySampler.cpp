#include "solver/RoughBoundarySampler.h"
#include "PhononMaterial.h"
#include "FluxWeightedWindow.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

namespace phonomc {
namespace {
using Vec3 = std::array<double, 3>;
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
Vec3 mul(const Vec3& a, double s) { return {a[0]*s, a[1]*s, a[2]*s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
double norm(const Vec3& a) { return std::sqrt(dot(a,a)); }
}
void RoughBoundarySampler::configure(const std::vector<RoughFacetSpec>& facets,
    const std::vector<const PhononMaterial*>& materials, std::ostream* progress,
    bool require_discrete_reflection, double equilibrium_temperature) {
    if (!std::isfinite(equilibrium_temperature) || equilibrium_temperature<0)
        throw std::invalid_argument("Invalid boundary equilibrium temperature.");
    for (const auto* mat : materials)
        if (mat == nullptr) throw std::invalid_argument("Null rough-boundary material.");
    for (const auto& facet : facets) {
        if (!facet.enabled) continue;
        if (!std::isfinite(facet.roughness) || facet.roughness < 0)
            throw std::invalid_argument("Invalid roughness.");
        for (double x : facet.outward_normal)
            if (!std::isfinite(x)) throw std::invalid_argument("Non-finite rough normal.");
        if (std::abs(norm(facet.outward_normal) - 1.0) > 1e-8)
            throw std::invalid_argument("Rough normal must have unit length.");
    }
    RoughBoundarySampler next;
    next.build(facets, materials, progress, require_discrete_reflection, equilibrium_temperature);
    *this = std::move(next);
}
int RoughBoundarySampler::table_index(int material, int facet) const {
    if (material < 0 || material >= static_cast<int>(active_counts_.size()) ||
        facet < 0 || facet >= facet_count_) return -1;
    return facet_to_rough_data_[static_cast<size_t>(material) * facet_count_ + facet];
}

// 函数说明：构建粗糙边界散射查找表，包含镜面率与漫反射候选映射。
void RoughBoundarySampler::build(const std::vector<RoughFacetSpec>& facets,
    const std::vector<const PhononMaterial*>& materials, std::ostream* progress,
    bool require_discrete_reflection, double equilibrium_temperature) {
    const int nfacets = static_cast<int>(facets.size());
    facet_count_ = nfacets;
    for (const auto* mat : materials) active_counts_.push_back(mat->active_mode_count());
    facet_to_rough_data_.assign(static_cast<size_t>(nfacets) * materials.size(), -1);
    rough_boundary_data_.clear();
    int rough_total = 0;
    for (int facet = 0; facet < nfacets; ++facet) {
        if (facets[static_cast<size_t>(facet)].enabled) ++rough_total;
    }
    if (rough_total == 0) {
        if (progress) *progress << "[init] Rough-boundary preprocessing skipped (0 rough facets).\n";
        return;
    }

    for (int material_id = 0; material_id < static_cast<int>(materials.size()); ++material_id) {
    const PhononMaterial& phonon = *materials[static_cast<size_t>(material_id)];
    const int na = phonon.active_mode_count();
    if (na <= 0) {
        continue;
    }
    const auto& active = phonon.active_mode_list();
    std::vector<Vec3> v_active(static_cast<size_t>(na), {0.0, 0.0, 0.0});
    std::vector<double> vnorm_active(static_cast<size_t>(na), 0.0);
    std::vector<double> omega_active(static_cast<size_t>(na), 0.0);
    std::vector<double> knorm_active(static_cast<size_t>(na), 0.0);
    std::vector<double> domega_active(static_cast<size_t>(na), 0.0);
    for (int ai = 0; ai < na; ++ai) {
        v_active[static_cast<size_t>(ai)] = phonon.mode_group_velocity(active[static_cast<size_t>(ai)]);
        vnorm_active[static_cast<size_t>(ai)] = norm(v_active[static_cast<size_t>(ai)]);
        omega_active[static_cast<size_t>(ai)] = phonon.mode_angular_frequency(active[static_cast<size_t>(ai)]);
        knorm_active[static_cast<size_t>(ai)] = phonon.mode_wavevector_norm(active[static_cast<size_t>(ai)]);
        domega_active[static_cast<size_t>(ai)] = phonon.mode_frequency_window(active[static_cast<size_t>(ai)]);
    }

    // Disjoint equal-frequency shells, so every diffuse row in a shell sees
    // the same residual distribution. Overlapping per-mode windows do not
    // in general preserve equilibrium.
    std::vector<int> shell(na), frequency_order(na);
    for (int i=0; i<na; ++i) frequency_order[i]=i;
    std::sort(frequency_order.begin(),frequency_order.end(),[&](int a,int b){return omega_active[a]<omega_active[b];});
    int shells=0;
    double shell_frequency=-1;
    for (int ai : frequency_order) {
        if (shell_frequency<0 || omega_active[ai]-shell_frequency>1e-8*std::max(shell_frequency,1e-12)) {
            ++shells; shell_frequency=omega_active[ai];
        }
        shell[ai]=shells-1;
    }
    int rough_done = 0;
    const int report_stride = std::max(1, rough_total / 10);

    for (int facet = 0; facet < nfacets; ++facet) {
        if (!facets[static_cast<size_t>(facet)].enabled) {
            continue;
        }
        ++rough_done;
        if (rough_done == 1 || rough_done == rough_total || (rough_done % report_stride) == 0) {
            if (progress) *progress << "[init]   rough facet " << rough_done << "/" << rough_total
                      << " (mesh facet id=" << facet << ")\n";
        }
        RoughFacetData rd;
        rd.facet = facet;
        rd.material_id = material_id;
        rd.strict_elastic = equilibrium_temperature>0;
        rd.specularity.assign(static_cast<size_t>(na), 0.0);
        rd.spec_match_active.assign(static_cast<size_t>(na), -1);
        rd.diffuse_begin.assign(static_cast<size_t>(na), -1);
        rd.diffuse_end.assign(static_cast<size_t>(na), -1);

        const double eta = std::max(0.0, facets[static_cast<size_t>(facet)].roughness);
        const Vec3 n_out = facets[static_cast<size_t>(facet)].outward_normal;
        const Vec3 n_in {-n_out[0], -n_out[1], -n_out[2]};
        std::vector<int> incoming;
        std::vector<int> outgoing;
        std::vector<double> incoming_destruction(static_cast<size_t>(na), 0.0);
        std::vector<double> outgoing_creation(static_cast<size_t>(na), 0.0);
        incoming.reserve(static_cast<size_t>(na));
        outgoing.reserve(static_cast<size_t>(na));

        for (int ai = 0; ai < na; ++ai) {
            const Vec3 v = v_active[static_cast<size_t>(ai)];
            const double vn = dot(v, n_in);
            if (vn < 0.0) {
                incoming.push_back(ai);
                incoming_destruction[static_cast<size_t>(ai)] = -vn;
            } else if (vn > 0.0) {
                outgoing.push_back(ai);
                outgoing_creation[static_cast<size_t>(ai)] = vn;
            }
            const double vnorm = std::max(1e-12, vnorm_active[static_cast<size_t>(ai)]);
            const double incidence_cos = std::abs(vn) / vnorm;
            const double x = rd.strict_elastic ? 2*eta*std::abs(dot(phonon.mode_wavevector(active[ai]),n_in))
                : 2.0 * eta * incidence_cos * knorm_active[static_cast<size_t>(ai)];
            double p = std::exp(-(x * x));
            if (!std::isfinite(p)) {
                p = 0.0;
            }
            // Specularity is only meaningful for incoming modes on this facet.
            rd.specularity[static_cast<size_t>(ai)] = (vn < 0.0) ? std::clamp(p, 0.0, 1.0) : 0.0;
            if (rd.strict_elastic) {
                const double c=phonon.mode_heat_capacity(equilibrium_temperature,active[ai]);
                incoming_destruction[ai]*=c;
                outgoing_creation[ai]*=c;
            }
        }
        if (require_discrete_reflection && !incoming.empty() && outgoing.empty())
            throw std::invalid_argument("Full-matrix reflecting facet " + std::to_string(facet) +
                " has incoming modes but no inward outgoing material mode for material " +
                std::to_string(material_id) + ". Supply a complete full-Brillouin-zone mode bank.");
        rd.outgoing_active = outgoing;
        rd.outgoing_sorted_active = outgoing;
        std::sort(
            rd.outgoing_sorted_active.begin(),
            rd.outgoing_sorted_active.end(),
            [&omega_active](int a, int b) { return omega_active[static_cast<size_t>(a)] < omega_active[static_cast<size_t>(b)]; });
        rd.outgoing_sorted_omega.resize(rd.outgoing_sorted_active.size());
        for (size_t k = 0; k < rd.outgoing_sorted_active.size(); ++k) {
            const int ao = rd.outgoing_sorted_active[k];
            rd.outgoing_sorted_omega[k] = omega_active[static_cast<size_t>(ao)];
        }

        std::vector<int> shell_begin(shells,-1),shell_end(shells,-1);
        if (rd.strict_elastic) {
            std::vector<double> in_flux(shells,0),out_flux(shells,0);
            for (int ai : incoming) in_flux[shell[ai]]+=incoming_destruction[ai];
            for (int ao : outgoing) out_flux[shell[ao]]+=outgoing_creation[ao];
            // A nominally stationary shell can contain ~1e-30 velocity noise.
            // Relative error alone assigns it a 100% imbalance. Allow only a
            // roundoff-sized absolute floor relative to resolved shell fluxes;
            // substantial asymmetry still fails and no inelastic fallback is added.
            const double resolved_flux=std::max(*std::max_element(in_flux.begin(),in_flux.end()),
                                                *std::max_element(out_flux.begin(),out_flux.end()));
            const double roundoff_flux=64*std::numeric_limits<double>::epsilon()*resolved_flux;
            for (int s=0; s<shells; ++s)
                if (std::abs(in_flux[s]-out_flux[s])>1e-7*std::max(in_flux[s],out_flux[s])+roundoff_flux)
                    throw std::invalid_argument("Callaway reflecting boundary lacks balanced equal-frequency incoming/outgoing flux. Check full-BZ material symmetry and velocities.");
            for (int pos=0; pos<static_cast<int>(rd.outgoing_sorted_active.size()); ++pos) {
                const int s=shell[rd.outgoing_sorted_active[pos]];
                if (shell_begin[s]<0) shell_begin[s]=pos;
                shell_end[s]=pos+1;
            }
        }
        for (int ai : incoming) {
            if (rd.outgoing_sorted_active.empty()) {
                continue;
            }
            const double w_in = omega_active[static_cast<size_t>(ai)];
            const double tol = std::max(1e-12, domega_active[static_cast<size_t>(ai)]);
            const double w_lo = w_in - tol;
            const double w_hi = w_in + tol;
            auto it_lo = std::lower_bound(rd.outgoing_sorted_omega.begin(), rd.outgoing_sorted_omega.end(), w_lo);
            auto it_hi = std::upper_bound(rd.outgoing_sorted_omega.begin(), rd.outgoing_sorted_omega.end(), w_hi);
            const int ib = rd.strict_elastic ? shell_begin[shell[ai]] : static_cast<int>(std::distance(rd.outgoing_sorted_omega.begin(), it_lo));
            const int ie = rd.strict_elastic ? shell_end[shell[ai]] : static_cast<int>(std::distance(rd.outgoing_sorted_omega.begin(), it_hi));
            if (ie <= ib) {
                continue;
            }
            rd.diffuse_begin[static_cast<size_t>(ai)] = ib;
            rd.diffuse_end[static_cast<size_t>(ai)] = ie;

            const Vec3 vin = v_active[static_cast<size_t>(ai)];
            const double vin_dot_n = dot(vin, n_in);
            const Vec3 vtry = sub(vin, mul(n_in, 2.0 * vin_dot_n));

            int best_ao = -1;
            double best_metric = std::numeric_limits<double>::infinity();
            for (int pos = ib; pos < ie; ++pos) {
                const int ao = rd.outgoing_sorted_active[static_cast<size_t>(pos)];
                const Vec3 vout = v_active[static_cast<size_t>(ao)];
                const Vec3 dv = sub(vtry, vout);
                const double refn = std::max({1e-12, norm(vtry), norm(vout)});
                const double rel = norm(dv) / refn;
                const double metric = rel;
                if (metric < best_metric) {
                    best_metric = metric;
                    best_ao = ao;
                }
            }
            if (best_ao >= 0) {
                rd.spec_match_active[static_cast<size_t>(ai)] = best_ao;
            }
        }

        // Population.py style: only truly matched incoming modes can be specular.
        for (int ai : incoming) {
            if (rd.spec_match_active[static_cast<size_t>(ai)] < 0) {
                rd.specularity[static_cast<size_t>(ai)] = 0.0;
            }
        }

        if (rd.strict_elastic) {
            // Limit each many-to-one specular column to the outgoing capacity.
            // Overflow becomes diffuse flux, instead of clipping a negative
            // residual and silently creating an angular equilibrium defect.
            std::vector<double> assigned(na,0);
            for (int ai : incoming) if (rd.spec_match_active[ai]>=0)
                assigned[rd.spec_match_active[ai]]+=incoming_destruction[ai]*rd.specularity[ai];
            for (int ai : incoming) {
                const int ao=rd.spec_match_active[ai];
                if (ao>=0 && assigned[ao]>outgoing_creation[ao])
                    rd.specularity[ai]*=outgoing_creation[ao]/assigned[ao];
            }
        }
        // Detailed-balance residual:
        // creation_rate(out) = C_total(out) - sum(specular_D(in -> out)).
        const auto diffuse_creation_rate =
            phonomc_detail::residual_diffuse_creation_rates(
                outgoing_creation,
                incoming_destruction,
                rd.specularity,
                rd.spec_match_active);

        double rate_sum = 0.0;
        for (int ao : outgoing) {
            const double r = diffuse_creation_rate[static_cast<size_t>(ao)];
            rate_sum += r;
        }
        rd.outgoing_sorted_residual_flux_prefix =
            phonomc_detail::ordered_nonnegative_prefix(
                diffuse_creation_rate,
                rd.outgoing_sorted_active);
        if (rate_sum > 0.0) {
            double cdf = 0.0;
            for (int ao : outgoing) {
                const double r = diffuse_creation_rate[static_cast<size_t>(ao)];
                if (r <= 0.0) {
                    continue;
                }
                cdf += r;
                rd.diffuse_roulette_active.push_back(ao);
                rd.diffuse_roulette_cdf.push_back(cdf);
            }
        }

        const int rid = static_cast<int>(rough_boundary_data_.size());
        rough_boundary_data_.push_back(std::move(rd));
        facet_to_rough_data_[static_cast<size_t>(material_id) * static_cast<size_t>(nfacets) +
            static_cast<size_t>(facet)] = rid;
    }
    }
}

// 函数说明：在粗糙边界条件下采样漫反射后的活跃模态索引。
// source: 1=precomputed window/roulette, 2=outgoing-pool fallback,
// 3=global-random fallback.
int RoughBoundarySampler::sample_diffuse(int rough_idx, int in_ai, int* source,
    std::mt19937_64& rng, RoughSampleEvents& events) const {
    const int rough_material = (rough_idx >= 0 && rough_idx < static_cast<int>(rough_boundary_data_.size()))
        ? rough_boundary_data_[static_cast<size_t>(rough_idx)].material_id : 0;
    const int na = active_counts_.at(static_cast<size_t>(rough_material));
    if (na <= 0) {
        if (source != nullptr) {
            *source = 3;
        }
        return 0;
    }
    if (rough_idx < 0 || rough_idx >= static_cast<int>(rough_boundary_data_.size())) {
        std::uniform_int_distribution<int> U(0, na - 1);
        if (source != nullptr) {
            *source = 3;
        }
        return U(rng);
    }
    const auto& rd = rough_boundary_data_[static_cast<size_t>(rough_idx)];

    // Adiabatic rough scattering is elastic. Draw within the incoming mode's
    // frequency window from the residual diffuse creation flux. The full
    // outgoing |v_g.n| distribution is not valid for partial specularity:
    // the specular channel has already supplied part of that equilibrium
    // outgoing flux and must be subtracted to preserve detailed balance.
    if (in_ai >= 0 && in_ai < static_cast<int>(rd.diffuse_begin.size()) &&
        in_ai < static_cast<int>(rd.diffuse_end.size())) {
        const int begin = rd.diffuse_begin[static_cast<size_t>(in_ai)];
        const int end = rd.diffuse_end[static_cast<size_t>(in_ai)];
        if (begin >= 0 && end > begin && end <= static_cast<int>(rd.outgoing_sorted_active.size())) {
            std::uniform_real_distribution<double> U01(0.0, 1.0);
            const int pos = phonomc_detail::flux_weighted_window_index(
                rd.outgoing_sorted_residual_flux_prefix,
                begin,
                end,
                U01(rng));
            if (pos >= begin && pos < end) {
                ++events.residual_window;
                if (source != nullptr) {
                    *source = 1;
                }
                return rd.outgoing_sorted_active[static_cast<size_t>(pos)];
            }
            // A zero-residual elastic window can occur because the discrete
            // mode list does not provide a complete frequency-shell match.
            // Use the global residual creation distribution before falling
            // back to an unbalanced uniform outgoing draw.
            ++events.residual_fallback;
        }
    }

    if (rd.strict_elastic) return -1; // no cross-frequency or uniform fallback
    if (!rd.diffuse_roulette_active.empty() &&
        !rd.diffuse_roulette_cdf.empty() &&
        rd.diffuse_roulette_active.size() == rd.diffuse_roulette_cdf.size() &&
        rd.diffuse_roulette_cdf.back() > 0.0) {
        std::uniform_real_distribution<double> U(0.0, rd.diffuse_roulette_cdf.back());
        const double r = U(rng);
        auto it = std::lower_bound(rd.diffuse_roulette_cdf.begin(), rd.diffuse_roulette_cdf.end(), r);
        size_t pos = static_cast<size_t>(std::distance(rd.diffuse_roulette_cdf.begin(), it));
        if (pos >= rd.diffuse_roulette_active.size()) {
            pos = rd.diffuse_roulette_active.size() - 1;
        }
        if (source != nullptr) {
            *source = 1;
        }
        return rd.diffuse_roulette_active[pos];
    }

    if (!rd.outgoing_active.empty()) {
        std::uniform_int_distribution<int> U(0, static_cast<int>(rd.outgoing_active.size()) - 1);
        if (source != nullptr) {
            *source = 2;
        }
        return rd.outgoing_active[static_cast<size_t>(U(rng))];
    }
    std::uniform_int_distribution<int> U(0, na - 1);
    if (source != nullptr) {
        *source = 3;
    }
    return U(rng);
}

// 函数说明：按镜面/漫反射概率选择碰撞后的模态与占据数。
RoughReflection RoughBoundarySampler::sample(
    const PhononMaterial& phonon,
    int rough_idx,
    const std::array<int, 2>& in_mode,
    double in_occupation, double reference_temperature, std::mt19937_64& rng, bool signed_weights) const {
    double out_occupation = in_occupation;
    RoughSampleEvents events;
    const int na = phonon.active_mode_count();
    ++events.events;
    if (rough_idx < 0 || rough_idx >= static_cast<int>(rough_boundary_data_.size()) || na <= 0) {
        ++events.missing_data;
        out_occupation = in_occupation;
        return {in_mode, out_occupation, events};
    }
    const auto& rd = rough_boundary_data_[static_cast<size_t>(rough_idx)];
    int in_ai = phonon.active_index_for_mode(in_mode);
    if (in_ai < 0 || in_ai >= na) {
        in_ai = -1;
    }

    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double p_spec = (in_ai >= 0 && in_ai < static_cast<int>(rd.specularity.size()))
        ? std::clamp(rd.specularity[static_cast<size_t>(in_ai)], 0.0, 1.0)
        : 0.0;
    if (in_ai >= 0 && U01(rng) <= p_spec) {
        const int out_ai = rd.spec_match_active[static_cast<size_t>(in_ai)];
        if (out_ai >= 0 && out_ai < na) {
            ++events.specular;
            const std::array<int, 2> out_mode = phonon.active_mode_at(out_ai);
            const double Tref = reference_temperature;
            const double win = phonon.mode_angular_frequency(in_mode);
            const double wout = phonon.mode_angular_frequency(out_mode);
            const double neq_in = phonon.bose_occupation(Tref, in_mode);
            const double neq_out = phonon.bose_occupation(Tref, out_mode);
            const double mapped = neq_out + win * (in_occupation - neq_in) / std::max(1e-30, wout);
            if (std::isfinite(mapped) && (signed_weights || mapped >= 0.0)) {
                out_occupation = mapped;
                return {out_mode, out_occupation, events};
            }
            out_occupation = in_occupation;
            return {in_mode, out_occupation, events};
        }
        ++events.missing_spec_match;
    }

    int source = 3;
    const int out_ai = sample_diffuse(rough_idx, in_ai, &source, rng, events);
    if (out_ai<0) return {in_mode,in_occupation,events}; // transport reports failed outgoing direction
    ++events.diffuse;
    if (source == 2) {
        ++events.outgoing_pool;
    } else if (source == 3) {
        ++events.global_random;
    }
    const std::array<int, 2> out_mode = phonon.active_mode_at(out_ai);
    const double Tref = reference_temperature;
    const double win = phonon.mode_angular_frequency(in_mode);
    const double wout = phonon.mode_angular_frequency(out_mode);
    const double neq_in = phonon.bose_occupation(Tref, in_mode);
    const double neq_out = phonon.bose_occupation(Tref, out_mode);
    const double mapped = neq_out + win * (in_occupation - neq_in) / std::max(1e-30, wout);
    if (std::isfinite(mapped) && (signed_weights || mapped >= 0.0)) {
        out_occupation = mapped;
        return {out_mode, out_occupation, events};
    }
    // A signed deviational mapping can be incompatible with a non-negative
    // physical occupation for a poorly matched mode.  Fall back to an elastic
    // reflection in the incoming mode rather than destroying energy.
    out_occupation = in_occupation;
    return {in_mode, out_occupation, events};
}

} // namespace phonomc
