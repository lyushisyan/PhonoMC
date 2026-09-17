#pragma once

#include "SimulationConfig.h"
#include <array>
#include <vector>
#include <random>
#include "solver/CellParticleIndex.h"

class PhononMaterial;

namespace phonomc {
class ParticleStorage;
class PrescribedHeatSource;

struct SourceGridView {
    const std::vector<int>& material_ids;
    const std::vector<double>& temperatures;
    const CellParticleIndex* cells = nullptr;
};

struct SourceDeposit {
    long long occupation_updates = 0;
    double energy_ev = 0.0;
    std::vector<int> appended_indices;
};

// Modal allocation only: does not calculate power, move particles, or write files.
// Weighted spectra are phenomenological energy weights per sampled carrier,
// not electron-phonon matrix elements or prescribed total branch fractions.
class PhononSourceOperator {
public:
    void configure(const SimulationConfig& config,
                   const std::vector<const PhononMaterial*>& materials);
    SourceDeposit apply(ParticleStorage& particles,
                        const std::vector<const PhononMaterial*>& materials,
                        const SourceGridView& grid, PrescribedHeatSource& source,
                        double particle_volume_a3, double background_temperature);
    SourceDeposit apply_full_matrix(ParticleStorage& particles,
                        const std::vector<const PhononMaterial*>& materials,
                        const SourceGridView& grid, PrescribedHeatSource& source,
                        double particle_volume_a3, double background_temperature,
                        std::mt19937_64& rng);

private:
    double mode_weight(const PhononMaterial& material, const std::array<int, 2>& mode) const;
    HeatSourceSpectrum spectrum_ = HeatSourceSpectrum::Thermal;
    double frequency_min_ = -1.0, frequency_max_ = -1.0;
    std::vector<int> branches_;
    std::vector<double> weights_;
    bool linearized_thermal_ = false;
    double gaussian_center_ = 0, gaussian_sigma_ = 1;
    // Complete-bank Gaussian energy fractions, indexed by material then active mode.
    // Independent of the number or distribution of Monte Carlo carriers.
    std::vector<std::vector<double>> gaussian_fractions_;
    CellParticleIndex fallback_cells_;
    std::vector<int> indices_;
    std::vector<double> increments_;
};

}  // namespace phonomc
