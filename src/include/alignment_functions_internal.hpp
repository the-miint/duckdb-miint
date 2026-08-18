#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>
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

// Individual CIGAR operation for positional walking (e.g., coverage depth calculation).
// Deliberately separate from ParseCigar() which returns aggregate CigarStats — that function
// is hot-path for sequence identity and adding vector allocation there would regress performance.
struct CigarOperation {
	char op;
	int64_t length;
};

// Parse CIGAR string into per-operation type+length pairs for positional walking.
// Uses same validation as ParseCigar() but returns individual operations instead of aggregates.
static inline std::vector<CigarOperation> ParseCigarOperations(const std::string &cigar_str) {
	std::vector<CigarOperation> ops;
	const char *cigar = cigar_str.data();
	size_t len = cigar_str.size();

	if (len == 0 || (len == 1 && cigar[0] == '*')) {
		return ops; // Empty or unmapped -- no allocation at all, and unmapped records are common
	}

	// Every operation needs at least two characters (a length digit and an op code), so this
	// is an upper bound -- loose on real CIGARs, whose lengths run to several digits, but it
	// spares the growth-by-doubling reallocations that otherwise dominate this function's
	// cost on long ones.
	ops.reserve(len / 2 + 1);

	int64_t op_len = 0;

	for (size_t i = 0; i < len; i++) {
		char c = cigar[i];

		if (std::isdigit(c)) {
			if (op_len > (INT64_MAX - 9) / 10) {
				throw InvalidInputException("CIGAR operation length exceeds maximum");
			}
			op_len = op_len * 10 + (c - '0');
		} else {
			if (op_len == 0) {
				throw InvalidInputException("Invalid CIGAR string: operation without length");
			}

			switch (c) {
			case 'M':
			case '=':
			case 'X':
			case 'I':
			case 'D':
			case 'N':
			case 'P':
			case 'S':
			case 'H':
				ops.push_back({c, op_len});
				break;
			default:
				throw InvalidInputException("Invalid CIGAR operation: " + std::string(1, c));
			}

			op_len = 0;
		}
	}

	if (op_len > 0) {
		throw InvalidInputException("Invalid CIGAR string: incomplete operation (missing operation character)");
	}

	return ops;
}

// A half-open [start, stop) run of query positions, 0-based, in the ORIGINAL READ's
// orientation (not the reference orientation the CIGAR is written in).
struct QueryInterval {
	int64_t start;
	int64_t stop;
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

// Decode the coverage `type` vocabulary shared by ComputeQueryCoverage and
// ComputeQueryIntervals: both count M/=/X, and "mapped" additionally counts I.
// Single source of truth for the accepted values and the rejection message, so the two
// cannot drift apart on what they accept.
static inline bool CoverageTypeCountsInsertions(const std::string &type) {
	if (type == "aligned") {
		return false;
	}
	if (type == "mapped") {
		return true;
	}
	throw InvalidInputException("Invalid coverage type: " + type + ". Must be 'aligned' or 'mapped'.");
}

// Compute query coverage from CIGAR statistics
// Returns the proportion of query bases covered by the reference
static inline double ComputeQueryCoverage(const CigarStats &stats, const std::string &type) {
	// Total query length (always includes hard clips)
	int64_t query_length = stats.matches + stats.insertions + stats.soft_clips + stats.hard_clips;

	if (query_length == 0) {
		return 0.0; // Avoid division by zero, return 0% coverage
	}

	// "aligned" counts only bases that align to the reference (M, =, X); "mapped" also
	// counts insertions, which are unclipped but align to nothing.
	int64_t covered_bases = stats.matches + (CoverageTypeCountsInsertions(type) ? stats.insertions : 0);

	return static_cast<double>(covered_bases) / static_cast<double>(query_length);
}

// Compute the query positions this alignment covers, as half-open intervals in the
// ORIGINAL READ's orientation.
//
// `type` selects what counts as covered, using the same vocabulary as
// ComputeQueryCoverage: "aligned" counts M/=/X, "mapped" additionally counts I.
// S and H advance the query cursor but are never covered; D/N/P consume no query, so
// runs they separate are contiguous on the query axis and are emitted merged.
//
// SAM writes CIGARs in reference orientation, so on a reverse-strand alignment the
// leading clip sits at the read's 3' end. `reverse` (FLAG 0x10) puts the returned intervals
// on the read's own axis regardless of orientation, so intervals from fragments of the same
// read are directly comparable -- which is the whole point of returning intervals rather
// than a count. Intervals are always returned ascending and non-overlapping.
//
// Empty for an unmapped or empty CIGAR. Throws on an unrecognised type.
static inline std::vector<QueryInterval> ComputeQueryIntervals(const std::vector<CigarOperation> &ops, bool reverse,
                                                               const std::string &type) {
	const bool count_insertions = CoverageTypeCountsInsertions(type);

	std::vector<QueryInterval> intervals;
	intervals.reserve(ops.size());
	int64_t cursor = 0; // offset from the read's own 5' end

	// SAM writes the CIGAR in reference orientation, so on a reverse-strand alignment the
	// last operation is the one at the read's 5' end. Walking the operations backwards
	// therefore makes the cursor a read-axis offset directly, and intervals come out
	// ascending with the run-merging below unchanged -- no mirror pass afterwards.
	for (size_t i = 0; i < ops.size(); i++) {
		const auto &op = ops[reverse ? ops.size() - 1 - i : i];
		bool consumes_query = true;
		bool covered = false;

		switch (op.op) {
		case 'M':
		case '=':
		case 'X':
			covered = true;
			break;
		case 'I':
			covered = count_insertions;
			break;
		case 'S':
		case 'H':
			break;
		case 'D':
		case 'N':
		case 'P':
			consumes_query = false;
			break;
		default:
			// ParseCigarOperations admits exactly the nine SAM ops, so reaching here means
			// the parser gained an operation this consumer was not taught about. Listed
			// exhaustively rather than defaulted to reference-only so that drift fails
			// loudly instead of silently mis-placing intervals.
			throw InvalidInputException("Unhandled CIGAR operation in query interval computation: " +
			                            std::string(1, op.op));
		}

		// ParseCigarOperations bounds each operation's length but not their sum. Without this
		// a hand-written CIGAR summing past INT64_MAX yields intervals with stop < start,
		// because the interval below is built from the cursor directly.
		if (consumes_query && cursor > INT64_MAX - op.length) {
			throw InvalidInputException("CIGAR query length exceeds maximum");
		}

		if (covered) {
			// Runs separated only by D/N/P are contiguous on the query axis; extend
			// rather than emitting an interval that abuts its predecessor.
			if (!intervals.empty() && intervals.back().stop == cursor) {
				intervals.back().stop = cursor + op.length;
			} else {
				intervals.push_back({cursor, cursor + op.length});
			}
		}
		if (consumes_query) {
			cursor += op.length;
		}
	}

	return intervals;
}

// The counters sequence identity is a function of, isolated from the rest of CigarStats so
// that identity can be pooled over the several alignment records one read was split into.
//
// Only these four are additive across records. gap_opens is not: an operation ending one
// record and one beginning the next are not adjacent in any alignment, so summing gap opens
// would count events that never happened. The clip counters describe an individual record's
// accounting rather than the molecule. Carrying either here would be a number that looks
// poolable and is not.
struct IdentityCounts {
	int64_t matches = 0;           // M, =, X
	int64_t match_ops = 0;         // = only
	int64_t mismatch_ops = 0;      // X only
	int64_t alignment_columns = 0; // M + I + D

	// ParseCigar bounds each operation's length but not the sum of them, so a CIGAR can
	// already carry counters near INT64_MAX, and pooling widens that exposure from one CIGAR
	// to a whole read group. Checked because an overflowed sum does not fail visibly: signed
	// overflow is undefined behaviour, and the observable result is an identity outside
	// [0.0, 1.0] -- a number that looks like an answer.
	void Add(const IdentityCounts &other) {
		matches = AddCounts(matches, other.matches);
		match_ops = AddCounts(match_ops, other.match_ops);
		mismatch_ops = AddCounts(mismatch_ops, other.mismatch_ops);
		alignment_columns = AddCounts(alignment_columns, other.alignment_columns);
	}

private:
	// Both operands are operation-length totals, so both are non-negative.
	static int64_t AddCounts(int64_t a, int64_t b) {
		if (a > INT64_MAX - b) {
			throw InvalidInputException("CIGAR operation counts exceed maximum");
		}
		return a + b;
	}
};

// Project a parsed CIGAR onto the counters identity is computed from.
static inline IdentityCounts ToIdentityCounts(const CigarStats &stats) {
	IdentityCounts counts;
	counts.matches = stats.matches;
	counts.match_ops = stats.match_ops;
	counts.mismatch_ops = stats.mismatch_ops;
	counts.alignment_columns = stats.alignment_columns;
	return counts;
}

// Compute sequence identity from extended-CIGAR operations (= and X) alone.
//
// Formula: match_ops / alignment_columns, where alignment_columns = M + I + D
// (= and X count toward the matches field).
//
// Returns std::nullopt when identity cannot be determined:
//   - Legacy CIGAR with only M ops (can't distinguish matches from mismatches)
//   - Mixed M alongside = or X (inconsistent encoding)
//   - Degenerate CIGARs with no =/X ops at all (pure I/D/S/H/N/P)
// SQL wrappers materialize std::nullopt as NULL.
//
// Taking counts rather than a whole CigarStats is what lets one read's records be pooled
// through the identical rules: the caller sums the counters and asks once. Applied to sums,
// the second rule also rejects a read whose records disagree about the encoding -- one in
// legacy M and one extended -- where the =/X counts describe only part of the molecule.
static inline std::optional<double> ComputeCigarIdentity(const IdentityCounts &counts) {
	if (counts.match_ops + counts.mismatch_ops == 0) {
		// No = or X ops observed (M-only or degenerate)
		return std::nullopt;
	}

	// If M ops are present alongside =/X, the encoding is inconsistent
	int64_t m_only = counts.matches - counts.match_ops - counts.mismatch_ops;
	if (m_only > 0) {
		return std::nullopt;
	}

	// alignment_columns is guaranteed > 0 here because =/X ops exist
	return static_cast<double>(counts.match_ops) / static_cast<double>(counts.alignment_columns);
}

static inline std::optional<double> ComputeCigarIdentity(const CigarStats &stats) {
	return ComputeCigarIdentity(ToIdentityCounts(stats));
}

} // namespace miint
