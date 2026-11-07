#include "MonteCarloSolver.h"

void MonteCarloSolver::initialize_interface_mode_banks() {
    phonomc::InterfaceInteractionConfig config;
    config.axis = material_layout_.axis();
    config.enabled = material_layout_.has_interfaces();
    config.linearized = uses_linearized_transport(args_);
    config.model = args_.material_interface_model;
    config.frequency_bin_thz = args_.material_interface_frequency_bin_thz;
    config.amm_fraction = args_.material_interface_amm_fraction;
    config.parallel_bin_inv_a = args_.material_interface_parallel_bin_inv_a;
    config.reference_temperature = background_temperature_;
    config.mass_densities_kg_m3 = args_.material_mass_densities_kg_m3;
    interface_interaction_.configure(config, materials_, material_energy_weights());
}
