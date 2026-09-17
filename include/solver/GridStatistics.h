#pragma once
#include "solver/BackgroundCache.h"
#include "solver/GridState.h"

namespace phonomc {
// Owns reusable reduction scratch, not a solver, geometry, RNG or output stream.
// Views must use the original material/cache mode ordering. Not concurrently
// reentrant: one statistics instance per solver. Internal loops may use OpenMP.
class GridStatistics {
public:
    // Refreshes modal energy/omega, cell density/temperature and particle
    // temperatures. Reuses GridState::cells; never changes occupations or modes.
    void refresh_temperatures(ParticleStorage& particles, GridState& grid,
        const std::vector<double>& volumes,
        const std::vector<const PhononMaterial*>& materials,
        const std::vector<BackgroundCache>& background, bool linearized_temperature = false);
    // Requires modal energies refreshed after the latest occupation/mode change.
    void update_heat_flux(const ParticleStorage& particles, GridState& grid,
        const std::vector<double>& volumes, const std::vector<BackgroundCache>& background);
    double total_thermal_energy(const GridState& grid, const std::vector<double>& volumes,
        const std::vector<const PhononMaterial*>& materials) const;
private:
    void ensure_tls_buffers(int thread_count, int cells);
    void update_energy_density(ParticleStorage& particles, GridState& grid,
        const std::vector<double>& volumes,
        const std::vector<const PhononMaterial*>& materials,
        const std::vector<BackgroundCache>& background);
    std::vector<double> energy_tls_buffer_;
    std::vector<Vec3> flux_tls_buffer_;
    int tls_thread_count_ = 0, tls_nsv_ = 0;
};
} // namespace phonomc
