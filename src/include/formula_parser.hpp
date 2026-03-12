#pragma once

#include <string>

namespace miint {

// Parse a chemical formula string and return the monoisotopic mass.
// Formula format: repeating [A-Z][a-z]*[0-9]* groups (e.g., "H2O", "C10H22", "NaCl")
// Throws std::runtime_error on invalid input.
double ParseFormula(const std::string &formula);

} // namespace miint
