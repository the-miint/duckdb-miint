#pragma once

#include <string>

namespace miint {

// Sum of monoisotopic residue masses for a single-letter amino acid sequence
double AminoacidDeltaMass(const std::string &sequence);

// Theoretical m/z of a peptide fragment ion
// ion_type: 'a', 'b', 'c' (N-terminal) or 'x', 'y', 'z' (C-terminal)
double PeptideFragmentMass(const std::string &sequence, int charge, char ion_type);

} // namespace miint
