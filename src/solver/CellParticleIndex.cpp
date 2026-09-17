#include "solver/CellParticleIndex.h"
#include "solver/ParticleStorage.h"

#include <cassert>
#include <stdexcept>

namespace phonomc {
bool CellParticleIndex::matches(const ParticleStorage& particles, int cell_count) const {
    return valid_ && owner_ == &particles && cell_count == static_cast<int>(counts_.size()) &&
        indices_.size() == static_cast<std::size_t>(particles.size());
}

void CellParticleIndex::ensure(const ParticleStorage& particles, int cell_count) {
    if (matches(particles, cell_count)) return;
    valid_ = false;
    if (cell_count <= 0) throw std::invalid_argument("Cell index requires a positive grid size.");
    assert(particles.aligned());
    counts_.assign(static_cast<std::size_t>(cell_count), 0);
    for (int cell : particles.grid_ids) {
        if (cell < 0 || cell >= cell_count)
            throw std::out_of_range("Particle references an invalid grid cell; refusing to clamp it.");
        ++counts_[static_cast<std::size_t>(cell)];
    }
    offsets_.resize(static_cast<std::size_t>(cell_count) + 1);
    offsets_[0] = 0;
    for (int cell = 0; cell < cell_count; ++cell)
        offsets_[cell + 1] = offsets_[cell] + counts_[cell];
    cursor_ = offsets_;
    indices_.resize(static_cast<std::size_t>(particles.size()));
    for (int i = 0; i < particles.size(); ++i)
        indices_[cursor_[particles.grid_ids[i]]++] = i;
    owner_ = &particles;
    valid_ = true;
    ++rebuild_count_;
}
}  // namespace phonomc
