#include "solver/GridStatistics.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#ifdef PHONOMC_USE_OPENMP
#include <omp.h>
#endif

namespace phonomc {
namespace {
Vec3 add(const Vec3& a, const Vec3& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
Vec3 mul(const Vec3& a, double s) { return {a[0]*s, a[1]*s, a[2]*s}; }
bool finite_vector(const Vec3& v) {
    return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}
void validate_mesh(const GridState& grid, const std::vector<double>& volumes, size_t materials) {
    if (grid.material_ids.empty() || grid.material_ids.size() > static_cast<size_t>(std::numeric_limits<int>::max()) ||
        grid.material_ids.size() != volumes.size()) throw std::invalid_argument("Invalid statistics grid dimensions.");
    for (size_t i = 0; i < volumes.size(); ++i) {
        if (!std::isfinite(volumes[i]) || volumes[i] <= 0 ||
            grid.material_ids[i] < 0 || static_cast<size_t>(grid.material_ids[i]) >= materials)
            throw std::invalid_argument("Invalid statistics cell volume/material.");
    }
}
void validate_materials(const std::vector<const PhononMaterial*>& materials) {
    for (const auto* material : materials)
        if (!material) throw std::invalid_argument("Null statistics material.");
}
void validate_background(const std::vector<BackgroundCache>& background) {
    for (const auto& b : background)
        if (!std::isfinite(b.energy_weight()) || b.energy_weight() < 0 || !std::isfinite(b.energy_density()))
            throw std::invalid_argument("Invalid statistics background cache.");
}
bool valid_particle(const ParticleStorage& p, int i, const GridState& grid,
    const std::vector<BackgroundCache>& background) {
    const int cell = p.grid_ids[i], material = p.material_ids[i];
    return cell >= 0 && static_cast<size_t>(cell) < grid.material_ids.size() &&
        material >= 0 && static_cast<size_t>(material) < background.size() &&
        material == grid.material_ids[cell];
}
bool valid_mode(const ParticleStorage& p, int i, const std::vector<const PhononMaterial*>& materials,
    const std::vector<BackgroundCache>& background) {
    const auto& mode = p.modes[i];
    const auto& material = *materials[p.material_ids[i]];
    return mode[0] >= 0 && mode[0] < material.qpoint_count() &&
        mode[1] >= 0 && mode[1] < material.branch_count() && background[p.material_ids[i]].contains_mode(mode);
}
} // namespace

// 函数说明：按线程数与网格数复用线程局部缓冲，降低 OpenMP 临时分配开销。
void GridStatistics::ensure_tls_buffers(int thread_count, int nsv) {
    thread_count = std::max(1, thread_count);
    nsv = std::max(1, nsv);
    if (thread_count == tls_thread_count_ && nsv == tls_nsv_) {
        return;
    }
    if (static_cast<size_t>(thread_count) > energy_tls_buffer_.max_size() / static_cast<size_t>(nsv) ||
        static_cast<size_t>(thread_count) > flux_tls_buffer_.max_size() / static_cast<size_t>(nsv))
        throw std::length_error("Statistics scratch dimensions overflow.");
    const size_t total = static_cast<size_t>(thread_count) * static_cast<size_t>(nsv);
    energy_tls_buffer_.assign(total, 0.0);
    flux_tls_buffer_.assign(total, Vec3 {0.0, 0.0, 0.0});
    tls_thread_count_ = thread_count;
    tls_nsv_ = nsv;
}

// 函数说明：由固定统计体积的粒子载体统计网格能量密度。
void GridStatistics::update_energy_density(ParticleStorage& particles, GridState& grid,
    const std::vector<double>& volumes, const std::vector<const PhononMaterial*>& materials,
    const std::vector<BackgroundCache>& background) {
    const int nsv = static_cast<int>(grid.material_ids.size());
    int invalid = 0;
    grid.energy_density.assign(static_cast<size_t>(nsv), 0.0);

#ifdef PHONOMC_USE_OPENMP
    const int thread_count = std::max(1, omp_get_max_threads());
    ensure_tls_buffers(thread_count, nsv);
    std::fill(energy_tls_buffer_.begin(), energy_tls_buffer_.end(), 0.0);
#pragma omp parallel reduction(|:invalid)
    {
        const int tid = omp_get_thread_num();
        double* local_energy = energy_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
#pragma omp for
        for (int i = 0; i < particles.size(); ++i) {
            if (!valid_particle(particles, i, grid, background) ||
                !valid_mode(particles, i, materials, background) || !std::isfinite(particles.occupation[i])) {
                invalid = 1; continue;
            }
            const PhononMaterial& phonon = *materials[particles.material_ids[i]];
            const int sv = particles.grid_ids[i];
            const double n_eq = background[particles.material_ids[i]].occupation(particles.modes[i]);
            const double dn = particles.occupation[i] - n_eq;
            particles.omega[i] = phonon.mode_angular_frequency(particles.modes[i]);
            particles.energies[i] = 6.582119569e-4 * particles.omega[i] * dn;  // hbar[eV*ps] * mode_angular_frequency[rad/ps] => eV
            local_energy[sv] += particles.energies[i] * background[particles.material_ids[i]].energy_weight();
        }
    }
    for (int tid = 0; tid < thread_count; ++tid) {
        const double* local = energy_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
        for (int sv = 0; sv < nsv; ++sv) {
            grid.energy_density[static_cast<size_t>(sv)] += local[sv];
        }
    }
#else
    for (int i = 0; i < particles.size(); ++i) {
        if (!valid_particle(particles, i, grid, background) ||
            !valid_mode(particles, i, materials, background) || !std::isfinite(particles.occupation[i])) {
            invalid = 1; continue;
        }
        const PhononMaterial& phonon = *materials[particles.material_ids[i]];
        const int sv = particles.grid_ids[i];
        const double n_eq = background[particles.material_ids[i]].occupation(particles.modes[i]);
        const double dn = particles.occupation[i] - n_eq;
        particles.omega[i] = phonon.mode_angular_frequency(particles.modes[i]);
        particles.energies[i] = 6.582119569e-4 * particles.omega[i] * dn;  // hbar[eV*ps] * mode_angular_frequency[rad/ps] => eV
        grid.energy_density[static_cast<size_t>(sv)] += particles.energies[i] * background[particles.material_ids[i]].energy_weight();
    }
#endif

    if (invalid) throw std::runtime_error("Invalid particle in grid energy statistics.");
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for reduction(|:invalid)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const double volume = std::max(1e-30, volumes[static_cast<size_t>(sv)]);
        // Each carrier keeps the spatial volume assigned at initialization.
        // Unlike M/N_g, this weight does not change when carriers cross a cell
        // boundary, so free flight cannot create or destroy represented energy.
        double e = grid.energy_density[static_cast<size_t>(sv)] / volume;
        e += background[grid.material_ids[static_cast<size_t>(sv)]].energy_density();
        grid.energy_density[static_cast<size_t>(sv)] = e;
        if (!std::isfinite(e)) invalid = 1;
    }
    if (invalid) throw std::overflow_error("Non-finite grid energy reduction.");
}

// 函数说明：完成粒子计数、网格温度反演与粒子温度回写闭环。
void GridStatistics::refresh_temperatures(ParticleStorage& particles, GridState& grid,
    const std::vector<double>& volumes, const std::vector<const PhononMaterial*>& materials,
    const std::vector<BackgroundCache>& background, bool linearized_temperature) {
    validate_mesh(grid, volumes, materials.size());
    validate_materials(materials);
    validate_background(background);
    if (linearized_temperature)
        for (const auto& b : background)
            if (!(b.heat_capacity_density() > 0.0) || !std::isfinite(b.heat_capacity_density()))
                throw std::invalid_argument("Linearized temperature requires positive heat capacity.");
    if (background.size() != materials.size()) throw std::invalid_argument("Background/material count mismatch.");
    if (!particles.aligned()) throw std::invalid_argument("Unaligned particle storage.");
    grid.temperatures.resize(grid.material_ids.size());
    const int nsv = static_cast<int>(grid.material_ids.size());
    grid.cells.ensure(particles, nsv);

    update_energy_density(particles, grid, volumes, materials, background);

#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const PhononMaterial& phonon = *materials[grid.material_ids[static_cast<size_t>(sv)]];
        const auto& reference = background[grid.material_ids[static_cast<size_t>(sv)]];
        grid.temperatures[static_cast<size_t>(sv)] = linearized_temperature
            ? reference.temperature() + (grid.energy_density[static_cast<size_t>(sv)] - reference.energy_density()) /
                reference.heat_capacity_density()
            : phonon.temperature_from_energy_density(grid.energy_density[static_cast<size_t>(sv)]);
    }
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles.size(); ++i) {
        const int sv = particles.grid_ids[i];
        particles.temperatures[i] = grid.temperatures[static_cast<size_t>(sv)];
    }
}

double GridStatistics::total_thermal_energy(const GridState& grid,
    const std::vector<double>& volumes, const std::vector<const PhononMaterial*>& materials) const {
    validate_mesh(grid, volumes, materials.size());
    validate_materials(materials);
    if (grid.energy_density.size() != volumes.size()) throw std::invalid_argument("Missing grid energy density.");
    for (double e : grid.energy_density)
        if (!std::isfinite(e)) throw std::invalid_argument("Non-finite grid energy density.");

    const size_t n = volumes.size();
    double total = 0.0;
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for reduction(+:total)
#endif
    for (int sv = 0; sv < static_cast<int>(n); ++sv) {
        const double zero_point_density = materials[grid.material_ids[static_cast<size_t>(sv)]]->zero_point_energy_density();
        total += (grid.energy_density[static_cast<size_t>(sv)] - zero_point_density) *
            volumes[static_cast<size_t>(sv)];
    }
    if (!std::isfinite(total)) throw std::overflow_error("Non-finite total thermal energy.");
    return total;
}

// Grid heat flux in W/m^2; conductivity estimates remain a caller responsibility.
void GridStatistics::update_heat_flux(const ParticleStorage& particles, GridState& grid,
    const std::vector<double>& volumes, const std::vector<BackgroundCache>& background) {
    validate_mesh(grid, volumes, background.size());
    validate_background(background);
    if (!particles.aligned()) throw std::invalid_argument("Unaligned particle storage.");
    const int nsv = static_cast<int>(grid.material_ids.size());
    int invalid = 0;
    grid.heat_flux.assign(static_cast<size_t>(nsv), {0.0, 0.0, 0.0});

#ifdef PHONOMC_USE_OPENMP
    const int thread_count = std::max(1, omp_get_max_threads());
    ensure_tls_buffers(thread_count, nsv);
    std::fill(flux_tls_buffer_.begin(), flux_tls_buffer_.end(), Vec3 {0.0, 0.0, 0.0});
#pragma omp parallel reduction(|:invalid)
    {
        const int tid = omp_get_thread_num();
        Vec3* local_flux = flux_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
#pragma omp for
        for (int i = 0; i < particles.size(); ++i) {
            if (!valid_particle(particles, i, grid, background) || !std::isfinite(particles.energies[i]) ||
                !finite_vector(particles.velocities[i])) { invalid = 1; continue; }
            const int sv = particles.grid_ids[i];
            local_flux[sv] = add(local_flux[sv], mul(
                particles.velocities[i],
                particles.energies[i] * background[particles.material_ids[i]].energy_weight()));
        }
    }
    for (int tid = 0; tid < thread_count; ++tid) {
        const Vec3* local = flux_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
        for (int sv = 0; sv < nsv; ++sv) {
            grid.heat_flux[static_cast<size_t>(sv)] =
                add(grid.heat_flux[static_cast<size_t>(sv)], local[sv]);
        }
    }
#else
    for (int i = 0; i < particles.size(); ++i) {
        if (!valid_particle(particles, i, grid, background) || !std::isfinite(particles.energies[i]) ||
            !finite_vector(particles.velocities[i])) { invalid = 1; continue; }
        const int sv = particles.grid_ids[i];
        grid.heat_flux[static_cast<size_t>(sv)] =
            add(grid.heat_flux[static_cast<size_t>(sv)], mul(
                particles.velocities[i],
                particles.energies[i] * background[particles.material_ids[i]].energy_weight()));
    }
#endif

    if (invalid) throw std::runtime_error("Invalid particle in grid heat-flux statistics.");
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for reduction(|:invalid)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const double volume = std::max(1e-30, volumes[static_cast<size_t>(sv)]);
        grid.heat_flux[static_cast<size_t>(sv)] =
            mul(grid.heat_flux[static_cast<size_t>(sv)], 1.0 / volume);
        grid.heat_flux[static_cast<size_t>(sv)] = mul(grid.heat_flux[static_cast<size_t>(sv)], 1.602176634e13);
        if (!finite_vector(grid.heat_flux[static_cast<size_t>(sv)])) invalid = 1;
    }
    if (invalid) throw std::overflow_error("Non-finite grid heat-flux reduction.");
}
} // namespace phonomc
