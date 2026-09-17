#pragma once

#include <cstddef>
#include <vector>

namespace phonomc {
class ParticleStorage;

// Stable cell buckets of stored particles (including any not yet compacted).
// Invalidate after transport, cell reassignment, insertion or compaction.
// Occupation/mode changes do not invalidate this index. Never mutate in a kernel.
class CellParticleIndex {
public:
    void invalidate() { valid_ = false; }
    void ensure(const ParticleStorage& particles, int cell_count);
    const std::vector<int>& counts() const { return counts_; }
    const std::vector<int>& offsets() const { return offsets_; }
    const std::vector<int>& indices() const { return indices_; }
    std::size_t rebuild_count() const { return rebuild_count_; }
    bool matches(const ParticleStorage& particles, int cell_count) const;
private:
    bool valid_ = false;
    const ParticleStorage* owner_ = nullptr;
    std::size_t rebuild_count_ = 0;
    std::vector<int> counts_, offsets_, cursor_, indices_;
};
}  // namespace phonomc
