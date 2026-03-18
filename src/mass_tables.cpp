#include "mass_tables.hpp"

#include <cctype>
#include <stdexcept>
#include <unordered_map>

namespace miint {

// Monoisotopic residue masses (same as pyteomics/Python MassQL)
static const std::unordered_map<char, double> AMINO_ACID_MASS = {
    {'G', 57.02146},  {'A', 71.03711},  {'V', 99.06841},  {'L', 113.08406}, {'I', 113.08406},
    {'P', 97.05276},  {'F', 147.06841}, {'W', 186.07931}, {'M', 131.04049}, {'S', 87.03203},
    {'T', 101.04768}, {'C', 103.00919}, {'Y', 163.06333}, {'H', 137.05891}, {'D', 115.02694},
    {'E', 129.04259}, {'N', 114.04293}, {'Q', 128.05858}, {'K', 128.09496}, {'R', 156.10111},
};

double AminoacidDeltaMass(const std::string &sequence) {
	if (sequence.empty()) {
		throw std::runtime_error("aminoaciddelta(): empty sequence");
	}
	double total = 0.0;
	for (char c : sequence) {
		char upper_c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
		auto it = AMINO_ACID_MASS.find(upper_c);
		if (it == AMINO_ACID_MASS.end()) {
			throw std::runtime_error("aminoaciddelta(): unknown amino acid '" + std::string(1, c) + "'");
		}
		total += it->second;
	}
	return total;
}

static constexpr double PROTON_MASS = 1.00728;
static constexpr double WATER_MASS = 18.01056; // H2O
static constexpr double CO_MASS = 27.99491;    // CO
static constexpr double NH3_MASS = 17.02655;   // NH3

double PeptideFragmentMass(const std::string &sequence, int charge, char ion_type) {
	if (charge < 1) {
		throw std::runtime_error("peptide(): charge must be >= 1");
	}
	double sum;
	try {
		sum = AminoacidDeltaMass(sequence);
	} catch (const std::runtime_error &e) {
		std::string msg = e.what();
		// Replace "aminoaciddelta():" prefix with "peptide():" so the error
		// references the function the user actually called.
		const std::string old_prefix = "aminoaciddelta():";
		if (msg.substr(0, old_prefix.size()) == old_prefix) {
			msg = "peptide():" + msg.substr(old_prefix.size());
		}
		throw std::runtime_error(msg);
	}

	// Ion type adjustments (neutral fragment mass)
	double neutral;
	switch (ion_type) {
	case 'b':
		neutral = sum;
		break;
	case 'y':
		neutral = sum + WATER_MASS;
		break;
	case 'a':
		neutral = sum - CO_MASS;
		break;
	case 'c':
		neutral = sum + NH3_MASS;
		break;
	case 'x':
		neutral = sum + CO_MASS + WATER_MASS - 2 * PROTON_MASS;
		break;
	case 'z':
		neutral = sum + WATER_MASS - NH3_MASS;
		break;
	default:
		throw std::runtime_error("peptide(): unknown ion type '" + std::string(1, ion_type) + "'");
	}

	// m/z = (neutral + charge * proton) / charge
	return (neutral + charge * PROTON_MASS) / charge;
}

} // namespace miint
