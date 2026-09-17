#pragma once
#include "solver/ResultSnapshots.h"
#include <ostream>
#include <string>
#include <utility>

namespace phonomc {
// Synchronous output only. No simulation, material, geometry, RNG or HDF5 dependency.
// Empty folder disables file output; stream formatting remains available.
class ResultWriter {
public:
    explicit ResultWriter(std::string folder = {}) : folder_(std::move(folder)) {}
    bool enabled() const { return !folder_.empty(); }
    void write_convergence_header(std::size_t cells, bool gradient_driven = false, bool write_cell_heat_flux = false) const;
    void append_convergence(const ConvergenceSnapshot& snapshot) const;
    void write_source_ledger(const SourceLedgerSnapshot& snapshot) const;
    static void format_convergence_header(std::ostream& out, std::size_t cells, bool gradient_driven = false);
    static void format_convergence(std::ostream& out, const ConvergenceSnapshot& snapshot);
    static void format_source_ledger(std::ostream& out, const SourceLedgerSnapshot& snapshot);
    static void format_source_summary(std::ostream& out, const SourceSummarySnapshot& snapshot);
    static void format_run_summary(std::ostream& out, const RunDiagnosticsSnapshot& snapshot);
private:
    std::string folder_;
};
}  // namespace phonomc
