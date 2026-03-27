// DUST low-complexity masking algorithm.
//
// This implementation follows the algorithm as implemented in vsearch v2.30.5
// (src/mask.cc, BSD-2-Clause license). The scoring formula, window parameters,
// and masking threshold are matched to produce equivalent masking results.
//
// Original DUST algorithm:
//   Tatusov RL, Lipman DJ. unpublished (1999).
//   Morgulis A, Gertz EM, Schaffer AA, Agarwala R.
//   "A fast and symmetric DUST implementation to mask low-complexity DNA sequences."
//   J Comput Biol. 2006;13(5):1028-40.
//
// vsearch reference:
//   Rognes T, Flouri T, Nichols B, Quince C, Mahe F.
//   "VSEARCH: a versatile open source tool for metagenomics."
//   PeerJ. 2016;4:e2584.

#include "DustMasker.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <vector>

namespace miint {

// Encode a single base to 2 bits: A/a=0, C/c=1, G/g=2, T/t/U/u=3.
// Non-nucleotide characters map to 0 (matching vsearch's behavior where
// ambiguous bases get a deterministic but arbitrary encoding).
static unsigned int base_to_2bit(char c) {
	switch (c) {
	case 'A':
	case 'a':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
	case 'U':
	case 'u':
		return 3;
	default:
		return 0; // Ambiguous bases map to 0 (like vsearch's map_2bit for unknown chars)
	}
}

// Find the highest-scoring sub-interval within a window.
// Follows vsearch's wo() function exactly:
//   - Scores: 10 * sum / j where j is the iteration counter (starts at dust_word-1)
//   - Returns score and sets besti (start offset) and bestj (end offset from start)
//   - Minimum region: j starts at dust_word-1 = 2, smallest evaluated interval is 3 words
static int dust_wo(const int *words, int len, int &besti, int &bestj) {
	besti = 0;
	bestj = 0;
	int bestv = 0;

	if (len < DUST_WORD) {
		return 0;
	}

	int num_words = len - DUST_WORD + 1;
	// Smallest possible region is 8 bases (5 words beyond the minimum 3-word unit)
	int l1 = num_words - 5;

	for (int i = 0; i < l1; i++) {
		int counts[64];
		std::memset(counts, 0, sizeof(counts));
		int sum = 0;

		// Start with the first trinucleotide at position i
		counts[words[i]]++;

		// Scan forward. j is the offset from the start position i.
		// j starts at dust_word-1 (=2) matching vsearch's loop bounds.
		for (int j = DUST_WORD - 1; j < num_words - i; j++) {
			int word = words[i + j];
			int c = counts[word];
			if (c != 0) {
				sum += c;
				// Score: 10 * sum / j (vsearch's formula, j is absolute iteration counter)
				int v = (10 * sum) / j;
				if (v > bestv) {
					bestv = v;
					besti = i;
					bestj = j;
				}
			}
			counts[word]++;
		}
	}

	return bestv;
}

void dust_mask(std::string &sequence) {
	int seq_len = static_cast<int>(sequence.size());

	if (seq_len < DUST_WORD) {
		for (auto &c : sequence) {
			c = std::toupper(static_cast<unsigned char>(c));
		}
		return;
	}

	// Pre-encode all trinucleotides for the entire sequence using rolling 2-bit encoding
	// (matching vsearch's approach)
	int num_words = seq_len - DUST_WORD + 1;
	std::vector<int> words(num_words);
	{
		unsigned int word = 0;
		unsigned int bitmask = (1U << (2 * DUST_WORD)) - 1; // 0x3F for 3-mers
		for (int j = 0; j < seq_len; j++) {
			word <<= 2U;
			word |= base_to_2bit(sequence[j]);
			if (j >= DUST_WORD - 1) {
				words[j - DUST_WORD + 1] = static_cast<int>(word & bitmask);
			}
		}
	}

	// Track which positions are masked
	std::vector<bool> masked(seq_len, false);

	int half_window = DUST_WINDOW / 2;

	// Slide window across the sequence (matching vsearch's dust() function)
	for (int window_start = 0; window_start < seq_len; window_start += half_window) {
		int window_len = std::min(DUST_WINDOW, seq_len - window_start);
		if (window_len < DUST_WORD) {
			break;
		}

		int besti, bestj;
		int score = dust_wo(&words[window_start], window_len, besti, bestj);

		if (score > DUST_LEVEL) {
			// Mask from besti to besti+bestj (inclusive, matching vsearch)
			int mask_from = window_start + besti;
			int mask_to = window_start + besti + bestj + DUST_WORD - 1;
			mask_to = std::min(mask_to, seq_len - 1);
			for (int j = mask_from; j <= mask_to; j++) {
				masked[j] = true;
			}
		}
	}

	// Apply masking: lowercase masked bases, uppercase unmasked bases
	for (int i = 0; i < seq_len; i++) {
		if (masked[i]) {
			sequence[i] = std::tolower(static_cast<unsigned char>(sequence[i]));
		} else {
			sequence[i] = std::toupper(static_cast<unsigned char>(sequence[i]));
		}
	}
}

} // namespace miint
