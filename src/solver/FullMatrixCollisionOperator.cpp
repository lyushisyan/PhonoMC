#include "solver/FullMatrixCollisionOperator.h"

#include "PhononMaterial.h"
#include "solver/ParticleStorage.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace phonomc {
namespace {
constexpr double kHbarEvPs = 6.582119569e-4;
constexpr double kTwoPi = 6.283185307179586476925286766559;

bool same_value(double a, double b) {
    return std::isfinite(a) && std::isfinite(b) &&
        std::abs(a - b) <= 1e-9 * std::max(std::abs(a), std::abs(b));
}
} // namespace

void FullMatrixCollisionOperator::configure(std::vector<ScatteringMatrix> matrices,
    const std::vector<const PhononMaterial*>& materials, double reference_temperature) {
    if (!std::isfinite(reference_temperature) || reference_temperature <= 0.0 ||
        matrices.empty() || matrices.size() != materials.size())
        throw std::invalid_argument("Full scattering requires one matrix per material and a positive fixed reference temperature.");
    std::vector<ModalCache> caches(materials.size());
    for (size_t m = 0; m < materials.size(); ++m) {
        if (materials[m] == nullptr) throw std::invalid_argument("Null full-scattering material.");
        const auto& material = *materials[m];
        const auto& data = matrices[m].data();
        if (!same_value(data.reference_temperature, reference_temperature))
            throw std::invalid_argument("Scattering matrix temperature differs from the fixed background temperature.");
        const auto n = matrices[m].size();
        if (n > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("Full-scattering mode count exceeds int indexing range.");
        auto& cache = caches[m];
        cache.flat_to_matrix.assign(static_cast<size_t>(material.qpoint_count()) * material.branch_count(), -1);
        cache.hw.resize(n);
        cache.background.resize(n);
        std::vector<long double> capacity(n), stationary(n, 0.0L), scale(n, 0.0L);
        for (size_t i = 0; i < n; ++i) {
            const auto& mode = data.modes[i];
            if (mode[0] < 0 || mode[0] >= material.qpoint_count() ||
                mode[1] < 0 || mode[1] >= material.branch_count())
                throw std::invalid_argument("Scattering matrix mode is outside the material q/branch grid.");
            const auto flat = static_cast<size_t>(mode[0]) * material.branch_count() + mode[1];
            const double omega = material.mode_angular_frequency(mode);
            if (!(omega > 0.0) || cache.flat_to_matrix[flat] >= 0 ||
                !same_value(data.frequency_thz[i], omega / kTwoPi))
                throw std::invalid_argument("Scattering matrix modes/frequencies do not match the positive-frequency material modes.");
            if (material.active_index_for_mode(mode) < 0)
                throw std::invalid_argument("Full scattering requires stationary positive-frequency modes in the material active bank.");
            cache.flat_to_matrix[flat] = static_cast<int>(i);
            cache.hw[i] = kHbarEvPs * omega;
            cache.background[i] = material.bose_occupation(reference_temperature, mode);
            capacity[i] = material.mode_heat_capacity(reference_temperature, mode);
            if (!std::isfinite(capacity[i]) || capacity[i] < 0.0L)
                throw std::invalid_argument("Non-finite full-scattering equilibrium heat capacity.");
        }
        for (int q = 0; q < material.qpoint_count(); ++q)
            for (int b = 0; b < material.branch_count(); ++b)
                if (material.mode_angular_frequency({q, b}) > 0.0 &&
                    cache.flat_to_matrix[static_cast<size_t>(q) * material.branch_count() + b] < 0)
                    throw std::invalid_argument("Full scattering matrix omits a positive-frequency material mode.");
        if (static_cast<size_t>(material.active_mode_count()) != n)
            throw std::invalid_argument("Full-scattering active bank must contain exactly the positive-frequency matrix modes.");
        for (size_t col = 0; col < n; ++col)
            for (size_t p = data.column_offsets[col]; p < data.column_offsets[col + 1]; ++p) {
                const auto row = data.row_indices[p];
                const long double term = data.values[p] * capacity[col];
                stationary[row] += term;
                scale[row] += std::abs(term);
            }
        const long double largest_scale = *std::max_element(scale.begin(), scale.end());
        for (size_t i = 0; i < n; ++i)
            if (std::abs(stationary[i]) > 1e-8L * scale[i] +
                64.0L * std::numeric_limits<double>::epsilon() * largest_scale)
                throw std::invalid_argument("Full scattering matrix does not preserve the equilibrium temperature perturbation (A*c != 0).");
    }
    matrices_ = std::move(matrices);
    callaway_.clear();
    materials_ = materials;
    modal_cache_ = std::move(caches);
    reference_temperature_ = reference_temperature;
    fallback_cells_.invalidate();
}

void FullMatrixCollisionOperator::configure_callaway(
    const std::vector<const PhononMaterial*>& materials, double reference_temperature) {
    configure_relaxation(materials, reference_temperature, false);
}

void FullMatrixCollisionOperator::configure_linearized_rta(
    const std::vector<const PhononMaterial*>& materials, double reference_temperature) {
    configure_relaxation(materials, reference_temperature, true);
}

void FullMatrixCollisionOperator::configure_relaxation(
    const std::vector<const PhononMaterial*>& materials, double reference_temperature, bool rta) {
    std::vector<CallawayOperator> kernels;
    std::vector<std::vector<double>> rates;
    std::vector<ScatteringMatrix> metadata;
    for (const auto* material : materials) {
        if (!material) throw std::invalid_argument("Null Callaway material.");
        const auto relaxation = rta ? material->linearized_rta_data(reference_temperature)
                                    : material->callaway_data(reference_temperature);
        kernels.emplace_back(relaxation);
        auto rate=relaxation.normal_rate;
        for (size_t i=0;i<rate.size();++i) rate[i]+=relaxation.resistive_rate[i];
        rates.push_back(std::move(rate));
        // An empty CSC holds the validated mode metadata only; no M*M storage.
        ScatteringMatrixData data;
        data.reference_temperature=reference_temperature;
        data.modes=material->active_mode_list();
        for (const auto& mode : data.modes)
            data.frequency_thz.push_back(material->mode_angular_frequency(mode)/kTwoPi);
        data.column_offsets.assign(data.modes.size()+1,0);
        metadata.emplace_back(std::move(data));
    }
    configure(std::move(metadata),materials,reference_temperature);
    callaway_=std::move(kernels);
    for (size_t i=0;i<rates.size();++i) modal_cache_[i].relaxation_rate=std::move(rates[i]);
}

FullMatrixCollisionResult FullMatrixCollisionOperator::apply(ParticleStorage& particles,
    const std::vector<const PhononMaterial*>& materials, const CollisionGridView& grid,
    double time_step_ps, double particle_volume_a3, double background_temperature,
    std::mt19937_64& rng, bool resolve_carrier_decay) {
    if (matrices_.empty() || materials != materials_)
        throw std::invalid_argument("Full scattering operator is not configured for these materials.");
    if (resolve_carrier_decay && callaway_.empty())
        throw std::invalid_argument("Carrier decay resolution requires a relaxation-model kernel.");
    if (!std::isfinite(time_step_ps) || time_step_ps < 0.0 ||
        !std::isfinite(particle_volume_a3) || particle_volume_a3 <= 0.0 ||
        !same_value(background_temperature, reference_temperature_))
        throw std::invalid_argument("Invalid full-scattering timestep, particle volume, or background temperature.");
    if (!particles.aligned() || grid.material_ids.size() != grid.particle_counts.size() ||
        grid.temperatures.size() != grid.particle_counts.size() || grid.particle_counts.empty() ||
        grid.particle_counts.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Invalid full-scattering particle/grid dimensions.");
    const int nsv = static_cast<int>(grid.particle_counts.size());
    const CellParticleIndex* cells = grid.cells;
    if (cells == nullptr) {
        fallback_cells_.invalidate();
        fallback_cells_.ensure(particles, nsv);
        cells = &fallback_cells_;
    }
    if (!cells->matches(particles, nsv) || cells->counts() != grid.particle_counts)
        throw std::invalid_argument("Full scattering requires a current cell index and matching counts.");

    // Stage the complete update so invalid input or an exponential failure
    // cannot leave occupations partly advanced. Index buckets remain unchanged
    // until every original cell has been evaluated.
    auto next_occupations = particles.occupation;
    std::vector<ParticleRecord> appended;
    long double residual = 0.0L;
    for (int cell = 0; cell < nsv; ++cell) {
        const int material_id = grid.material_ids[cell];
        if (material_id < 0 || material_id >= static_cast<int>(materials.size()) ||
            !std::isfinite(grid.temperatures[cell]))
            throw std::invalid_argument("Invalid full-scattering cell material or temperature.");
        const auto& material = *materials[material_id];
        const auto& cache = modal_cache_[material_id];
        const auto& matrix = matrices_[material_id];
        const int begin = cells->offsets()[cell], end = cells->offsets()[cell + 1];
        if (begin == end) continue;
        std::vector<double> energy(matrix.size(), 0.0);
        std::vector<int> counts(matrix.size(), 0);
        std::vector<size_t> particle_modes(static_cast<size_t>(end - begin));
        for (int pos = begin; pos < end; ++pos) {
            const int i = cells->indices()[pos];
            const auto& mode = particles.modes[i];
            if (particles.grid_ids[i] != cell || particles.material_ids[i] != material_id ||
                particles.alive[i] == 0 || !std::isfinite(particles.occupation[i]) ||
                mode[0] < 0 || mode[0] >= material.qpoint_count() ||
                mode[1] < 0 || mode[1] >= material.branch_count())
                throw std::invalid_argument("Invalid particle in full-scattering cell.");
            const int index = cache.flat_to_matrix[static_cast<size_t>(mode[0]) * material.branch_count() + mode[1]];
            if (index < 0) throw std::invalid_argument("Particle mode is missing from the full scattering matrix.");
            const auto mode_index = static_cast<size_t>(index);
            particle_modes[static_cast<size_t>(pos - begin)] = mode_index;
            energy[mode_index] += cache.hw[mode_index] * (particles.occupation[i] - cache.background[mode_index]);
            ++counts[mode_index];
        }
        for (double value : energy)
            if (!std::isfinite(value)) throw std::overflow_error("Full-scattering modal energy overflow.");
        if (time_step_ps == 0.0 || (!resolve_carrier_decay &&
            std::all_of(energy.begin(), energy.end(), [](double e) { return e == 0.0; }))) continue;
        const auto evolved = callaway_.empty() ? matrix.evolve(energy, time_step_ps)
            : callaway_.at(static_cast<size_t>(material_id)).evolve(energy,time_step_ps);
        std::vector<double> increment(matrix.size(), 0.0), attenuation(matrix.size(),1.0);
        if (resolve_carrier_decay)
            for (size_t mode=0;mode<matrix.size();++mode)
                attenuation[mode]=std::exp(-cache.relaxation_rate[mode]*time_step_ps);
        for (size_t mode = 0; mode < matrix.size(); ++mode)
            if (counts[mode] > 0)
                increment[mode] = (evolved[mode] - attenuation[mode]*energy[mode]) / cache.hw[mode] / counts[mode];
        long double cell_residual = 0.0L;
        for (int pos = begin; pos < end; ++pos) {
            const int i = cells->indices()[pos];
            const auto mode = particle_modes[static_cast<size_t>(pos - begin)];
            // Preserve the exact cell-modal exponential while also attenuating
            // differences between carriers of the same mode at different positions.
            // Only the collision gain is shared uniformly within that modal cell.
            const double next = resolve_carrier_decay
                ? cache.background[mode]+(particles.occupation[i]-cache.background[mode])*attenuation[mode]+increment[mode]
                : particles.occupation[i]+increment[mode];
            if (!std::isfinite(next)) throw std::overflow_error("Full scattering produced a non-finite signed occupation.");
            next_occupations[i] = next;
            cell_residual += cache.hw[mode] * (static_cast<long double>(next) - particles.occupation[i]);
        }
        for (size_t mode = 0; mode < matrix.size(); ++mode) {
            if (counts[mode] > 0 || evolved[mode] == 0.0) continue;
            if (appended.size() >= static_cast<size_t>(std::numeric_limits<int>::max() - particles.size()))
                throw std::length_error("Full-scattering carrier population exceeds int range.");
            std::uniform_int_distribution<int> choose(begin, end - 1);
            const int donor = cells->indices()[choose(rng)];
            ParticleRecord p;
            p.mode = matrix.data().modes[mode];
            p.material_id = material_id;
            p.grid_id = cell;
            p.position = particles.positions[donor];
            for (double x : p.position)
                if (!std::isfinite(x)) throw std::invalid_argument("Non-finite full-scattering carrier position.");
            p.collision_position = p.position;
            p.velocity = material.mode_group_velocity(p.mode);
            p.temperature = grid.temperatures[cell];
            p.omega = material.mode_angular_frequency(p.mode);
            p.occupation = cache.background[mode] + evolved[mode] / cache.hw[mode];
            if (!std::isfinite(p.occupation))
                throw std::overflow_error("Full scattering produced a non-finite signed carrier.");
            p.energy = cache.hw[mode] * (p.occupation - cache.background[mode]);
            cell_residual += p.energy;
            appended.push_back(p);
        }
        const long double energy_weight = static_cast<double>(material.active_mode_count()) *
            particles.spatial_volume(material_id,particle_volume_a3) / material.energy_density_normalization();
        if (!std::isfinite(energy_weight) || energy_weight <= 0.0L)
            throw std::overflow_error("Invalid full-scattering carrier energy weight.");
        residual += cell_residual * energy_weight;
    }
    FullMatrixCollisionResult result;
    result.residual_ev = static_cast<double>(residual);
    if (!std::isfinite(result.residual_ev)) throw std::overflow_error("Full-scattering residual overflow.");
    result.appended_indices.reserve(appended.size());
    // append() maintains aligned storage even if memory allocation fails.
    for (const auto& p : appended) result.appended_indices.push_back(particles.append(p));
    std::copy(next_occupations.begin(), next_occupations.end(), particles.occupation.begin());
    fallback_cells_.invalidate();
    return result;
}
} // namespace phonomc
