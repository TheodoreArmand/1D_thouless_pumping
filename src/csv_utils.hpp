#pragma once

#include <string>
#include <vector>

namespace ecg1d {

// Simple CSV convention used by this project:
// - fields are separated by literal commas
// - quoted fields and escaped commas are not supported
// - whitespace is preserved, not trimmed
//
// Example phase schedule rows:
//   s,phi
//   0.0,0.0
//   0.5,3.141592653589793
std::vector<std::string> split_simple_csv_line(const std::string& line);

}  // namespace ecg1d
