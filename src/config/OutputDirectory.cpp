#include "SimulationConfig.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>

// 函数说明：按序号创建结果目录，避免覆盖历史仿真输出。
std::string create_indexed_output_folder(const std::string& output_folder) {
    if (output_folder.empty()) {
        throw std::runtime_error("output_folder cannot be empty.");
    }

    namespace fs = std::filesystem;
    fs::path base_abs = fs::absolute(fs::path(output_folder)).lexically_normal();
    if (base_abs.filename().empty()) {
        base_abs = base_abs.parent_path();
    }
    const fs::path parent = base_abs.parent_path();
    const auto stem = base_abs.filename().string();

    int max_idx = -1;
    if (fs::exists(parent) && fs::is_directory(parent)) {
        for (const auto& ent : fs::directory_iterator(parent)) {
            if (!ent.is_directory()) {
                continue;
            }
            const auto name = ent.path().filename().string();
            const auto prefix = stem + "_";
            if (name.rfind(prefix, 0) != 0) {
                continue;
            }
            const auto suffix = name.substr(prefix.size());
            if (suffix.empty()) {
                continue;
            }
            bool all_digits = std::all_of(suffix.begin(), suffix.end(), [](unsigned char ch) {
                return std::isdigit(ch);
            });
            if (!all_digits) {
                continue;
            }
            max_idx = std::max(max_idx, std::stoi(suffix));
        }
    }

    const int next_idx = max_idx + 1;
    const fs::path out = parent / (stem + "_" + std::to_string(next_idx));
    fs::create_directories(out);
    return out.string();
}
