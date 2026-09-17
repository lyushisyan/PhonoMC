#include "SimulationConfig.h"
#include "SimulationDomain.h"
#include "PhononMaterial.h"
#include "material/MaterialFactory.h"
#include "MonteCarloSolver.h"

#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <functional>
#include <memory>
#include <string>
#include <stdexcept>
#include <vector>

static void print_startup_banner() {
    static constexpr std::size_t kInnerWidth = 68;
    const auto print_border = []() {
        std::cout << "+" << std::string(kInnerWidth + 2, '-') << "+\n";
    };
    const auto print_line = [](const std::string& content) {
        std::cout << "| " << content;
        if (content.size() < kInnerWidth) {
            std::cout << std::string(kInnerWidth - content.size(), ' ');
        }
        std::cout << " |\n";
    };

    print_border();
    print_line("  _____  _                       __  __  _____");
    print_line(" |  __ \\| |                     |  \\/  |/ ____|");
    print_line(" | |__) | |__   ___  _ __   ___ | \\  / | |");
    print_line(" |  ___/| '_ \\ / _ \\| '_ \\ / _ \\| |\\/| | |");
    print_line(" | |    | | | | (_) | | | | (_) | |  | | |____");
    print_line(" |_|    |_| |_|\\___/|_| |_|\\___/|_|  |_|\\_____|");
    print_line("");
    print_line("                              PhonoMC");
    print_border();
    std::cout << std::flush;
}

// 函数说明：程序入口：读取输入、构建几何与材料、执行时间步并输出汇总结果。
int main(int argc, char** argv) {
    try {
        namespace fs = std::filesystem;
        print_startup_banner();

        std::string input_file;
        if (argc > 1) {
            input_file = argv[1];
        } else {
            const std::vector<fs::path> cwd_candidates {
                fs::path("input.toml")
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
                    exe_dir / "../input.toml",
                    exe_dir / "../../input.toml"
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
        const int material_count = args.material_folders.empty()
            ? 1 : static_cast<int>(args.material_folders.size());
        std::vector<std::unique_ptr<PhononMaterial>> material_storage;
        std::vector<std::reference_wrapper<const PhononMaterial>> materials;
        material_storage.reserve(static_cast<size_t>(material_count));
        materials.reserve(static_cast<size_t>(material_count));
        for (int material_index = 0; material_index < material_count; ++material_index) {
            material_storage.push_back(std::make_unique<PhononMaterial>(phonomc::load_phonon_material(args, material_index)));
            materials.push_back(std::cref(*material_storage.back()));
        }
        MonteCarloSolver pop(args, geo, materials);

        while (pop.current_timestep() < args.iterations && !pop.converged()) {
            pop.run_timestep();
        }

        const auto end = std::chrono::steady_clock::now();
        if (pop.converged()) std::cout << "Stopped after global and layer heat-flux stationarity windows passed.\n";
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
            if (!summary) throw std::runtime_error("Cannot append runtime summary.txt.");
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
            pop.append_profile_summary(summary);
            summary.close();
            if (!summary) throw std::runtime_error("Failed flushing runtime summary.txt.");
        }
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Fatal error: " << ex.what() << '\n';
        return 1;
    }
}
