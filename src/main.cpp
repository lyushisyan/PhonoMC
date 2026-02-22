#include "SimulationConfig.h"
#include "SimulationDomain.h"
#include "PhononMaterial.h"
#include "MonteCarloSolver.h"

#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
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
        if (args.heat_source_enabled && args.compute_kappa) {
            std::cout << "Config note: heat_source.enabled=true, so thermal conductivity estimation is disabled for this run.\n";
            args.compute_kappa = false;
        }
        if (!args.output_folder.empty()) {
            fs::path r(args.output_folder);
            if (!r.is_absolute()) {
                args.output_folder = (fs::path(args.input_directory) / r).lexically_normal().string();
            }
        }
        args.output_folder = create_indexed_output_folder(args.output_folder);
        std::cout << "Results folder: " << args.output_folder << '\n';

        const auto start = std::chrono::steady_clock::now();

        SimulationDomain geo(args);
        PhononMaterial phonons(args, 0);
        MonteCarloSolver pop(args, geo, phonons);

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

        if (!args.output_folder.empty()) {
            std::ofstream summary(fs::path(args.output_folder) / "summary.txt", std::ios::app);
            if (summary) {
                summary << "\n[runtime]\n";
                summary << "total_seconds = " << sec << '\n';
                summary << "total_days = " << days << '\n';
                summary << "total_hours = " << hours << '\n';
                summary << "total_minutes = " << minutes << '\n';
                summary << "total_seconds_remainder = " << seconds << '\n';
                summary << "total_human_readable = "
                        << days << " days "
                        << hours << " h "
                        << minutes << " min "
                        << seconds << " s\n";
            }
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Fatal error: " << ex.what() << '\n';
        return 1;
    }
}
