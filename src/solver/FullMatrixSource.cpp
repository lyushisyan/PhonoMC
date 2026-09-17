#include "solver/PhononSourceOperator.h"
#include "solver/ParticleStorage.h"
#include "solver/PrescribedHeatSource.h"
#include "PhononMaterial.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace phonomc {
SourceDeposit PhononSourceOperator::apply_full_matrix(
    ParticleStorage& particles, const std::vector<const PhononMaterial*>& materials,
    const SourceGridView& grid, PrescribedHeatSource& source,
    double particle_volume_a3, double background_temperature, std::mt19937_64& rng) {
    SourceDeposit result;
    if (!source.enabled() || source.pending_energy_ev() <= 0.0) return result;
    const auto& pending = source.cell_pending_ev();
    if (!particles.aligned() || grid.material_ids.size() != pending.size() ||
        grid.temperatures.size() != pending.size() || !(background_temperature > 0.0))
        throw std::invalid_argument("Invalid full-matrix source grid/reference.");
    const CellParticleIndex* cells = grid.cells;
    if (!cells) {
        fallback_cells_.invalidate();
        fallback_cells_.ensure(particles, static_cast<int>(pending.size()));
        cells = &fallback_cells_;
    }
    if (!cells->matches(particles, static_cast<int>(pending.size())))
        throw std::invalid_argument("Full-matrix source requires a current cell index.");
    // Buckets refer to the original population; append never changes their contents.
    for (size_t cell = 0; cell < pending.size(); ++cell) {
        const int begin = cells->offsets()[cell], end = cells->offsets()[cell+1];
        if (!(pending[cell] > 0.0)) continue;
        const int material_id = grid.material_ids[cell];
        const auto* material_ptr = materials.at(static_cast<size_t>(material_id));
        if (!material_ptr) throw std::invalid_argument("Null full-matrix source material.");
        const auto& mat = *material_ptr;
        const auto& modes = mat.active_mode_list();
        const double energy_weight = modes.size() * particles.spatial_volume(material_id,particle_volume_a3) / mat.energy_density_normalization();
        if (!(energy_weight > 0.0) || !std::isfinite(energy_weight))
            throw std::invalid_argument("Invalid full-matrix source energy weight.");
        std::vector<std::vector<int>> carriers(modes.size());
        for (int pos = begin; pos < end; ++pos) {
            const int i = cells->indices()[pos];
            const auto mode = particles.modes[i];
            if (particles.material_ids[i] != material_id || !particles.alive[i] ||
                mode[0] < 0 || mode[0] >= mat.qpoint_count() || mode[1] < 0 || mode[1] >= mat.branch_count())
                throw std::invalid_argument("Invalid full-matrix source carrier.");
            const int index = mat.active_index_for_mode(mode);
            if (index < 0 || !std::isfinite(particles.occupation[i]))
                throw std::invalid_argument("Invalid full-matrix source mode/weight.");
            carriers[index].push_back(i);
        }
        std::vector<double> allocation(modes.size(), 0.0);
        long double sum = 0.0L;
        for (size_t m = 0; m < modes.size(); ++m) {
            if (spectrum_ == HeatSourceSpectrum::Gaussian) {
                allocation[m] = gaussian_fractions_.at(static_cast<size_t>(material_id)).at(m);
                sum += allocation[m];
                continue;
            }
            const double weight = mode_weight(mat, modes[m]);
            if (weight > 0.0 && mat.mode_angular_frequency(modes[m]) > 0.0)
                allocation[m] = spectrum_ == HeatSourceSpectrum::Thermal
                    ? mat.mode_heat_capacity(background_temperature, modes[m]) : weight;
            sum += allocation[m];
        }
        if (!(sum > 0.0L)) continue;
        const double target = pending[cell] / energy_weight;
        if (!std::isfinite(target)) throw std::overflow_error("Full-matrix source energy overflow.");
        double assigned = 0.0;
        const size_t largest = static_cast<size_t>(std::max_element(allocation.begin(), allocation.end())-allocation.begin());
        for (double& e : allocation) { e = target * static_cast<double>(e / sum); assigned += e; }
        allocation[largest] += target-assigned;
        std::vector<std::pair<int,double>> updates;
        std::vector<ParticleRecord> additions;
        for (size_t m = 0; m < modes.size(); ++m) {
            if (allocation[m] == 0.0) continue;
            const double hw = 6.582119569e-4 * mat.mode_angular_frequency(modes[m]);
            if (!carriers[m].empty()) {
                const double dn = allocation[m] / hw / carriers[m].size();
                for (int i : carriers[m]) {
                    const double next = particles.occupation[i]+dn;
                    if (!std::isfinite(next)) throw std::overflow_error("Full-matrix source weight overflow.");
                    updates.emplace_back(i,next);
                }
            } else {
                ParticleRecord p;
                p.mode = modes[m]; p.material_id = material_id; p.grid_id = static_cast<int>(cell);
                // An empty deviational cell still contains material. Seed it
                // from the configured source geometry, without requiring a donor.
                p.position = begin==end ? source.cell_source_position(cell) :
                    particles.positions[cells->indices()[std::uniform_int_distribution<int>(begin,end-1)(rng)]];
                p.collision_position = p.position;
                p.velocity = mat.mode_group_velocity(p.mode); p.omega = mat.mode_angular_frequency(p.mode);
                p.temperature = grid.temperatures[cell];
                p.occupation = mat.bose_occupation(background_temperature,p.mode)+allocation[m]/hw;
                p.energy = allocation[m];
                if (!std::isfinite(p.occupation)) throw std::overflow_error("Full-matrix source weight overflow.");
                additions.push_back(p);
            }
        }
        if (additions.size() > static_cast<size_t>(std::numeric_limits<int>::max()-particles.size()))
            throw std::overflow_error("Full-matrix source exceeds carrier capacity.");
        for (const auto& update : updates) particles.occupation[update.first] = update.second;
        for (const auto& p : additions) result.appended_indices.push_back(particles.append(p));
        result.occupation_updates += updates.size()+additions.size();
        result.energy_ev += source.deposit_cell(static_cast<int>(cell));
    }
    return result;
}
} // namespace phonomc
