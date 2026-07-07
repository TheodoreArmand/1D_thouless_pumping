#include "csv_utils.hpp"

#include <sstream>

namespace ecg1d {

std::vector<std::string> split_simple_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }
    return fields;
}

}  // namespace ecg1d
