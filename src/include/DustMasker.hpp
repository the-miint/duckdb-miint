#pragma once

#include <string>

namespace miint {

// DUST low-complexity masking algorithm.
//
// Identifies low-complexity regions in DNA/RNA sequences by counting
// overlapping trinucleotide (3-mer) repeats in sliding windows. Regions
// with high trinucleotide repetition are "masked" by lowercasing their
// bases. Unmasked bases are uppercased.
//
// This is critical for k-mer indexing: k-mers containing any lowercase
// (masked) base should be skipped, preventing low-complexity regions
// (e.g., AAAAAAA, ACACAC) from flooding the index with non-informative hits.
//
// Algorithm: sliding window of 64bp, advancing 32bp per step. Within each
// window, finds the highest-scoring sub-interval using trinucleotide counts.
// Score = 10 * sum(c_t * (c_t-1) / 2) / (L-2) where c_t is the count of
// each trinucleotide type. Regions scoring above the threshold are masked.
//
// Based on the DUST algorithm as implemented in vsearch (BSD-2-Clause).
// Default threshold matches vsearch and NCBI BLAST: level=20 (scaled).

static constexpr int DUST_WINDOW = 64;
static constexpr int DUST_LEVEL = 20;
static constexpr int DUST_WORD = 3;

// Apply DUST masking to a sequence in-place.
// Masked (low-complexity) bases are lowercased.
// Unmasked bases are uppercased.
// Handles DNA (ACGT) and RNA (U treated as T).
void dust_mask(std::string &sequence);

} // namespace miint
