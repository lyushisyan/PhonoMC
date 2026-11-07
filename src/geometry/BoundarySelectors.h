#pragma once
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace phonomc::geometry_detail {
using Box6 = std::array<double, 6>;

// Preserve legacy selector aliases and traversal order for assignment and output.
inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

struct BoundarySelectors {
    std::vector<Box6> boxes;
};

// Relative region sextuplets; retain the explicitly supported selector aliases.
inline BoundarySelectors parse_boundary_selectors(const std::vector<std::string>& tokens,
                                                  const std::string& name) {
    BoundarySelectors out;
    if (tokens.empty()) return out;
    const std::string mode = to_lower(tokens.front());
    const bool named = mode == "box" || mode == "boxes" || mode == "region" ||
                       mode == "regions" || mode == "relative_box" ||
                       mode == "relative_boxes" || mode == "relative_region" ||
                       mode == "relative_regions";
    const size_t start = named ? 1 : 0;
    if ((tokens.size() - start) % 6 != 0) {
        throw std::runtime_error(name + " regions require sextuplets [xmin,ymin,zmin,xmax,ymax,zmax].");
    }
    out.boxes.reserve((tokens.size() - start) / 6);
    for (size_t i = start; i < tokens.size(); i += 6) {
        Box6 box;
        for (size_t j = 0; j < box.size(); ++j) {
            const auto& token = tokens[i + j];
            size_t consumed = 0;
            try {
                box[j] = std::stod(token, &consumed);
            } catch (const std::invalid_argument&) {
                throw std::runtime_error(name + " contains a nonnumeric coordinate: " + token);
            } catch (const std::out_of_range&) {
                throw std::runtime_error(name + " coordinate is out of range: " + token);
            }
            if (consumed != token.size() || !std::isfinite(box[j])) {
                throw std::runtime_error(name + " requires finite numeric coordinates: " + token);
            }
        }
        out.boxes.push_back(box);
    }
    return out;
}

} // namespace phonomc::geometry_detail
