#pragma once

#include <array>
#include <stdexcept>
#include <string>

namespace miint {

// DNA complement lookup table for all 256 ASCII characters
// 0 means invalid character
static constexpr std::array<char, 256> CreateDnaComplementTable() {
	std::array<char, 256> table = {};
	for (size_t i = 0; i < 256; i++) {
		table[i] = 0;
	}
	table['A'] = 'T';
	table['T'] = 'A';
	table['G'] = 'C';
	table['C'] = 'G';
	table['R'] = 'Y';
	table['Y'] = 'R';
	table['S'] = 'S';
	table['W'] = 'W';
	table['K'] = 'M';
	table['M'] = 'K';
	table['B'] = 'V';
	table['D'] = 'H';
	table['H'] = 'D';
	table['V'] = 'B';
	table['N'] = 'N';
	table['a'] = 't';
	table['t'] = 'a';
	table['g'] = 'c';
	table['c'] = 'g';
	table['r'] = 'y';
	table['y'] = 'r';
	table['s'] = 's';
	table['w'] = 'w';
	table['k'] = 'm';
	table['m'] = 'k';
	table['b'] = 'v';
	table['d'] = 'h';
	table['h'] = 'd';
	table['v'] = 'b';
	table['n'] = 'n';
	table['-'] = '-';
	table['.'] = '.';
	return table;
}

// RNA complement lookup table for all 256 ASCII characters
// 0 means invalid character
static constexpr std::array<char, 256> CreateRnaComplementTable() {
	std::array<char, 256> table = {};
	for (size_t i = 0; i < 256; i++) {
		table[i] = 0;
	}
	table['A'] = 'U';
	table['U'] = 'A';
	table['G'] = 'C';
	table['C'] = 'G';
	table['R'] = 'Y';
	table['Y'] = 'R';
	table['S'] = 'S';
	table['W'] = 'W';
	table['K'] = 'M';
	table['M'] = 'K';
	table['B'] = 'V';
	table['D'] = 'H';
	table['H'] = 'D';
	table['V'] = 'B';
	table['N'] = 'N';
	table['a'] = 'u';
	table['u'] = 'a';
	table['g'] = 'c';
	table['c'] = 'g';
	table['r'] = 'y';
	table['y'] = 'r';
	table['s'] = 's';
	table['w'] = 'w';
	table['k'] = 'm';
	table['m'] = 'k';
	table['b'] = 'v';
	table['d'] = 'h';
	table['h'] = 'd';
	table['v'] = 'b';
	table['n'] = 'n';
	table['-'] = '-';
	table['.'] = '.';
	return table;
}

static constexpr auto DNA_COMPLEMENT_TABLE = CreateDnaComplementTable();
static constexpr auto RNA_COMPLEMENT_TABLE = CreateRnaComplementTable();

// Reverse complement using a given complement table
template <const std::array<char, 256> &COMPLEMENT_TABLE>
inline std::string reverse_complement(const std::string &seq) {
	std::string rc(seq.size(), '\0');
	for (size_t i = 0; i < seq.size(); i++) {
		unsigned char base = static_cast<unsigned char>(seq[seq.size() - 1 - i]);
		char comp = COMPLEMENT_TABLE[base];
		if (comp == 0) {
			throw std::runtime_error(std::string("Invalid character: '") + static_cast<char>(base) + "'");
		}
		rc[i] = comp;
	}
	return rc;
}

inline std::string dna_reverse_complement(const std::string &seq) {
	return reverse_complement<DNA_COMPLEMENT_TABLE>(seq);
}

inline std::string rna_reverse_complement(const std::string &seq) {
	return reverse_complement<RNA_COMPLEMENT_TABLE>(seq);
}

} // namespace miint
