#pragma once

#include <cstdint>
#include <string>
#include <cctype>
#include <stdexcept>

namespace miint {

// Simple exception class for parsing errors
// This is used internally and can be caught/rethrown as DuckDB exceptions in the implementation
class InvalidInputException : public std::runtime_error {
public:
	explicit InvalidInputException(const std::string &msg) : std::runtime_error(msg) {
	}
};

// CIGAR parsing result structure
struct CigarStats {
	int64_t matches = 0;           // M, =, X operations
	int64_t match_ops = 0;         // = operations only
	int64_t mismatch_ops = 0;      // X operations only
	int64_t insertions = 0;        // I operations
	int64_t deletions = 0;         // D operations
	int64_t gap_opens = 0;         // Number of gap opening events
	int64_t alignment_columns = 0; // M + I + D operations
	int64_t soft_clips = 0;        // S operations
	int64_t hard_clips = 0;        // H operations
};

// MD parsing result structure
struct MdStats {
	int64_t matches = 0;
	int64_t mismatches = 0;
};

// Parse CIGAR string and extract statistics
static inline CigarStats ParseCigar(const std::string &cigar_str) {
	CigarStats stats;
	const char *cigar = cigar_str.data();
	size_t len = cigar_str.size();

	if (len == 0 || (len == 1 && cigar[0] == '*')) {
		return stats; // Empty or unmapped
	}

	int64_t op_len = 0;
	char prev_op_type = '\0'; // Track previous operation type for gap opens

	for (size_t i = 0; i < len; i++) {
		char c = cigar[i];

		if (std::isdigit(c)) {
			// Check for integer overflow before multiplication
			if (op_len > (INT64_MAX - 9) / 10) {
				throw InvalidInputException("CIGAR operation length exceeds maximum");
			}
			op_len = op_len * 10 + (c - '0');
		} else {
			// Operation character
			if (op_len == 0) {
				throw InvalidInputException("Invalid CIGAR string: operation without length");
			}

			switch (c) {
			case 'M': // Match or mismatch (alignment match)
				stats.matches += op_len;
				stats.alignment_columns += op_len;
				break;
			case '=': // Sequence match
				stats.matches += op_len;
				stats.match_ops += op_len;
				stats.alignment_columns += op_len;
				break;
			case 'X': // Sequence mismatch
				stats.matches += op_len;
				stats.mismatch_ops += op_len;
				stats.alignment_columns += op_len;
				break;
			case 'I': // Insertion to the reference
				stats.insertions += op_len;
				stats.alignment_columns += op_len;
				// Gap open counting: consecutive I operations count as a single gap event.
				// For example, "5I3I" = 1 gap open, but "5I3M2I" = 2 gap opens.
				// An insertion following a deletion (or vice versa) counts as separate events.
				if (prev_op_type != 'I') {
					stats.gap_opens++;
				}
				break;
			case 'D': // Deletion from the reference
				stats.deletions += op_len;
				stats.alignment_columns += op_len;
				// Gap open counting: consecutive D operations count as a single gap event.
				// For example, "5D3D" = 1 gap open, but "5D3M2D" = 2 gap opens.
				if (prev_op_type != 'D') {
					stats.gap_opens++;
				}
				break;
			case 'N': // Skipped region (ref skip, e.g., intron)
			case 'P': // Padding (silent deletion from padded reference)
				// These operations don't consume query or contribute to alignment columns.
				// N is used for spliced alignments (RNA-seq), P is for padded alignments.
				// Per SAM spec, both are ignored in identity and length calculations.
				break;
			case 'S': // Soft clipping
				stats.soft_clips += op_len;
				break;
			case 'H': // Hard clipping
				stats.hard_clips += op_len;
				break;
			default:
				throw InvalidInputException("Invalid CIGAR operation: " + std::string(1, c));
			}

			prev_op_type = c;
			op_len = 0;
		}
	}

	// Check for incomplete CIGAR (digits without operation at end)
	if (op_len > 0) {
		throw InvalidInputException("Invalid CIGAR string: incomplete operation (missing operation character)");
	}

	return stats;
}

// Parse MD tag string and extract match/mismatch counts
static inline MdStats ParseMd(const std::string &md_str) {
	MdStats stats;
	const char *md = md_str.data();
	size_t len = md_str.size();

	if (len == 0) {
		return stats; // Empty MD tag
	}

	int64_t match_len = 0;

	for (size_t i = 0; i < len; i++) {
		char c = md[i];

		if (std::isdigit(c)) {
			// Check for integer overflow before multiplication
			if (match_len > (INT64_MAX - 9) / 10) {
				throw InvalidInputException("MD tag match length exceeds maximum");
			}
			match_len = match_len * 10 + (c - '0');
		} else if (c == '^') {
			// Deletion marker: skip following deleted bases
			if (match_len > 0) {
				stats.matches += match_len;
				match_len = 0;
			}
			i++; // Skip the '^'
			// Skip deletion bases
			while (i < len && std::isalpha(md[i])) {
				i++;
			}
			i--; // Back up one since loop will increment
		} else if (std::isalpha(c)) {
			// Mismatch base
			if (match_len > 0) {
				stats.matches += match_len;
				match_len = 0;
			}
			stats.mismatches++;
		}
	}

	// Add remaining matches
	if (match_len > 0) {
		stats.matches += match_len;
	}

	return stats;
}

// Compute query length from CIGAR statistics
// When include_hard_clips=false, result matches HTSlib's bam_cigar2qlen
// (which counts M, I, S, =, X operations but excludes H, D, N, P)
static inline int64_t ComputeQueryLength(const CigarStats &stats, bool include_hard_clips) {
	// Query-consuming operations: M, I, S, =, X
	// matches field contains M + = + X
	int64_t length = stats.matches + stats.insertions + stats.soft_clips;

	if (include_hard_clips) {
		length += stats.hard_clips;
	}

	return length;
}

// Compute query coverage from CIGAR statistics
// Returns the proportion of query bases covered by the reference
static inline double ComputeQueryCoverage(const CigarStats &stats, const std::string &type) {
	// Total query length (always includes hard clips)
	int64_t query_length = stats.matches + stats.insertions + stats.soft_clips + stats.hard_clips;

	if (query_length == 0) {
		return 0.0; // Avoid division by zero, return 0% coverage
	}

	int64_t covered_bases = 0;
	if (type == "aligned") {
		// Only bases that align to reference (M, =, X)
		covered_bases = stats.matches;
	} else if (type == "mapped") {
		// Bases that are mapped (not clipped): M, I, =, X
		covered_bases = stats.matches + stats.insertions;
	} else {
		throw InvalidInputException("Invalid coverage type: " + type +
		                                    ". Must be 'aligned' or 'mapped'.");
	}

	return static_cast<double>(covered_bases) / static_cast<double>(query_length);
}

} // namespace miint
