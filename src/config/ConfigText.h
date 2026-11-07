#pragma once
#include <algorithm>
#include <cctype>
#include <string>

namespace phonomc::config_detail {
// 函数说明：清理输入文本两端空白，保证配置与数据解析的稳健性。
inline std::string trim(const std::string& s) {
    const auto first = std::find_if_not(s.begin(), s.end(), [](unsigned char ch) {
        return std::isspace(ch);
    });
    if (first == s.end()) {
        return "";
    }
    const auto last = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char ch) {
                          return std::isspace(ch);
                      }).base();
    return std::string(first, last);
}

// 函数说明：将关键字统一为小写，避免输入大小写差异影响流程。
inline std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return s;
}

} // namespace phonomc::config_detail
