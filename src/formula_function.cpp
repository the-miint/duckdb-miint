// Monoisotopic element masses derived from pyteomics (Apache License 2.0).
// Source: pyteomics/auxiliary/constants.py, _nist_mass dictionary, key 0 (most abundant isotope).
// Repository: https://github.com/levitsky/pyteomics
// Commit: 3f1fd4afb51a5033222851666bef585c9253cd68
// See THIRD_PARTY_LICENSES.md for full license and citation information.

#include "formula_function.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <string>
#include <unordered_map>

namespace duckdb {

// Monoisotopic masses from pyteomics _nist_mass[element][0][0].
// These are the mass of the most abundant isotope for each element.
// Values must match pyteomics exactly for reproducibility with MassQL.
static const std::unordered_map<std::string, double> ELEMENT_MASSES = {
    {"H", 1.00782503207},  {"He", 4.00260325415}, {"Li", 7.01600455},    {"Be", 9.0121822},     {"B", 11.0093054},
    {"C", 12.0},           {"N", 14.0030740048},  {"O", 15.99491461956}, {"F", 18.99840322},    {"Ne", 19.9924401754},
    {"Na", 22.9897692809}, {"Mg", 23.9850417},    {"Al", 26.98153863},   {"Si", 27.9769265325}, {"P", 30.97376163},
    {"S", 31.972071},      {"Cl", 34.96885268},   {"Ar", 39.9623831225}, {"K", 38.96370668},    {"Ca", 39.96259098},
    {"Mn", 54.9380451},    {"Fe", 55.9349375},    {"Co", 58.933195},     {"Ni", 57.9353429},    {"Cu", 62.9295975},
    {"Zn", 63.9291422},    {"Se", 79.9165213},    {"Br", 78.9183371},    {"Mo", 97.9054082},    {"I", 126.904473},
};

// Parse a chemical formula string and return the monoisotopic mass.
// Formula format: repeating [A-Z][a-z]*[0-9]* groups (e.g., "H2O", "C10H22", "NaCl")
static double ParseFormula(const std::string &formula) {
	if (formula.empty()) {
		throw InvalidInputException("formula(): empty formula string");
	}

	double total_mass = 0.0;
	size_t i = 0;
	size_t len = formula.size();

	while (i < len) {
		// Must start with uppercase letter
		if (!std::isupper(static_cast<unsigned char>(formula[i]))) {
			throw InvalidInputException("formula(): unexpected character '%c' at position %zu in '%s'", formula[i], i,
			                            formula.c_str());
		}

		// Parse element symbol: [A-Z][a-z]*
		size_t elem_start = i;
		i++; // consume uppercase
		while (i < len && std::islower(static_cast<unsigned char>(formula[i]))) {
			i++;
		}
		std::string element = formula.substr(elem_start, i - elem_start);

		// Look up element mass
		auto it = ELEMENT_MASSES.find(element);
		if (it == ELEMENT_MASSES.end()) {
			throw InvalidInputException("formula(): unknown element '%s' in '%s'", element.c_str(), formula.c_str());
		}

		// Parse optional count: [0-9]*
		size_t count_start = i;
		while (i < len && std::isdigit(static_cast<unsigned char>(formula[i]))) {
			i++;
		}
		int count = 1;
		if (i > count_start) {
			std::string count_str = formula.substr(count_start, i - count_start);
			try {
				count = std::stoi(count_str);
			} catch (const std::out_of_range &) {
				throw InvalidInputException("formula(): element count '%s' is too large in '%s'", count_str.c_str(),
				                            formula.c_str());
			}
			if (count <= 0) {
				throw InvalidInputException("formula(): element count must be positive in '%s'", formula.c_str());
			}
		}

		total_mass += it->second * count;
	}

	return total_mass;
}

static void FormulaScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(),
	                                         [](string_t input) { return ParseFormula(input.GetString()); });
}

void FormulaFunction::Register(ExtensionLoader &loader) {
	ScalarFunction formula_func("formula", {LogicalType::VARCHAR}, LogicalType::DOUBLE, FormulaScalarFunction);
	loader.RegisterFunction(formula_func);
}

} // namespace duckdb
