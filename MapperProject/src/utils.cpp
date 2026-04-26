#include "utils.hpp"
#include <vector>
#include <string>

std::vector<std::string> splitString(const std::string &input, char delimiter) {
    std::vector<std::string> result;
    size_t start = 0;
    size_t end = input.find(delimiter);
    while (end != std::string::npos) {
        result.push_back(input.substr(start, end - start));
        start = end + 1;
        end = input.find(delimiter, start);
    }
    result.push_back(input.substr(start)); // Add the last segment
    return result;
}