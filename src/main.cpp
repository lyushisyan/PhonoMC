#include "SimulationConfig.h"
#include "Geometry.h"
#include "Phonon.h"
#include "Simulation.h"

#include <chrono>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        namespace fs = std::filesystem;

        std::string input_file;
        if (argc > 1) {
            input_file = argv[1];
        } else {
            const std::vector<fs::path> cwd_candidates {
                fs::path("input.toml"),
                fs::path("input.txt")
            };
            bool found = false;
            for (const auto& c : cwd_candidates) {
                if (fs::exists(c)) {
                    input_file = c.string();
                    found = true;
                    break;
                }
            }
            if (!found) {
                fs::path exe = fs::weakly_canonical(fs::path(argv[0]));
                fs::path exe_dir = exe.has_parent_path() ? exe.parent_path() : fs::current_path();
                const std::vector<fs::path> candidates {
                    exe_dir / "input.toml",
                    exe_dir / "input.txt",
                    exe_dir / "../input.toml",
                    exe_dir / "../input.txt",
                    exe_dir / "../../input.toml",
                    exe_dir / "../../input.txt"
                };
                for (const auto& c : candidates) {
                    if (fs::exists(c)) {
                        input_file = fs::weakly_canonical(c).string();
                        found = true;
                        break;
                    }
                }
            }
            if (!found) {
                input_file = "input.toml";
            }
        }

        SimulationConfig args = load_simulation_config(input_file);
        if (args.heat_source_enabled && args.compute_thermal_conductivity) {
            std::cout << "Config note: heat_source.enabled=true, so thermal conductivity estimation is disabled for this run.\n";
            args.compute_thermal_conductivity = false;
        }
        if (!args.results_base_folder.empty()) {
            fs::path r(args.results_base_folder);
            if (!r.is_absolute()) {
                args.results_base_folder = (fs::path(args.input_directory) / r).lexically_normal().string();
            }
        }
        args.results_base_folder = create_indexed_results_folder(args.results_base_folder);
        std::cout << "Results folder: " << args.results_base_folder << '\n';

        const auto start = std::chrono::steady_clock::now();

        Geometry geo(args);
        Phonon phonons(args, 0);
        Simulation pop(args, geo, phonons);

        while (pop.current_timestep() < args.iterations) {
            pop.run_timestep();
        }

        const auto end = std::chrono::steady_clock::now();
        const auto sec = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
        const auto days = sec / (24 * 3600);
        const auto hours = (sec % (24 * 3600)) / 3600;
        const auto minutes = (sec % 3600) / 60;
        const auto seconds = sec % 60;

        std::cout << "Total time: " << days << " days "
                  << hours << " h "
                  << minutes << " min "
                  << seconds << " s\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Fatal error: " << ex.what() << '\n';
        return 1;
    }
}
