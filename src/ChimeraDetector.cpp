#include "ChimeraDetector.hpp"

#include <cctype>
#include <stdexcept>

namespace miint {

static bool is_gap(char c) {
	return c == '-';
}

// Returns true for IUPAC ambiguity codes (N, R, Y, S, W, K, M, B, D, H, V)
// and any non-ACGT non-gap character. Gap ('-') is NOT ambiguous — it has
// its own handling in classify_diffs.
static bool is_ambiguous(char c) {
	char upper = std::toupper(static_cast<unsigned char>(c));
	switch (upper) {
	case 'A':
	case 'C':
	case 'G':
	case 'T':
		return false;
	default:
		return !is_gap(c);
	}
}

static char to_upper(char c) {
	return std::toupper(static_cast<unsigned char>(c));
}

StarAlignment build_star_alignment(const std::string &query_aligned_a, const std::string &subject_aligned_a,
                                   const std::string &query_aligned_b, const std::string &subject_aligned_b) {
	StarAlignment result;

	// Walk both pairwise alignments using ungapped query positions as the anchor.
	// Both alignments are global (alignEnd2End) on the same query, so they consume
	// exactly the same query bases. The only structural difference is where gaps
	// appear (insertions in one subject but not the other).
	//
	// When one alignment has a gap in the query (= insertion in that subject),
	// the other arm gets a gap column inserted to maintain positional correspondence.

	size_t ia = 0; // position in alignment A
	size_t ib = 0; // position in alignment B

	while (ia < query_aligned_a.size() && ib < query_aligned_b.size()) {
		bool a_query_gap = is_gap(query_aligned_a[ia]);
		bool b_query_gap = is_gap(query_aligned_b[ib]);

		if (a_query_gap && b_query_gap) {
			// Both have insertion at this position
			result.query_row += '-';
			result.parent_a_row += subject_aligned_a[ia];
			result.parent_b_row += subject_aligned_b[ib];
			ia++;
			ib++;
		} else if (a_query_gap) {
			// Insertion in A's subject only — gap in B
			result.query_row += '-';
			result.parent_a_row += subject_aligned_a[ia];
			result.parent_b_row += '-';
			ia++;
		} else if (b_query_gap) {
			// Insertion in B's subject only — gap in A
			result.query_row += '-';
			result.parent_a_row += '-';
			result.parent_b_row += subject_aligned_b[ib];
			ib++;
		} else {
			// Both have a query base — emit the column
			result.query_row += query_aligned_a[ia];
			result.parent_a_row += subject_aligned_a[ia];
			result.parent_b_row += subject_aligned_b[ib];
			ia++;
			ib++;
		}
	}

	// Drain remaining insertions from whichever alignment hasn't finished.
	// After the main loop, one alignment may have trailing insertions (query gaps).
	// These are subject-only insertions that must be represented.
	while (ia < query_aligned_a.size()) {
		// Remaining columns in A must be query gaps (subject insertions)
		if (!is_gap(query_aligned_a[ia])) {
			throw std::runtime_error("build_star_alignment: alignment A has unconsumed query bases after B finished");
		}
		result.query_row += '-';
		result.parent_a_row += subject_aligned_a[ia];
		result.parent_b_row += '-';
		ia++;
	}
	while (ib < query_aligned_b.size()) {
		if (!is_gap(query_aligned_b[ib])) {
			throw std::runtime_error("build_star_alignment: alignment B has unconsumed query bases after A finished");
		}
		result.query_row += '-';
		result.parent_a_row += '-';
		result.parent_b_row += subject_aligned_b[ib];
		ib++;
	}

	return result;
}

void classify_diffs(StarAlignment &star) {
	size_t len = star.query_row.size();
	star.diffs.resize(len, DiffType::IGNORE);

	for (size_t i = 0; i < len; i++) {
		char q = to_upper(star.query_row[i]);
		char a = to_upper(star.parent_a_row[i]);
		char b = to_upper(star.parent_b_row[i]);

		// Skip if any position has a gap
		if (is_gap(q) || is_gap(a) || is_gap(b)) {
			star.diffs[i] = DiffType::IGNORE;
			continue;
		}

		// Skip if adjacent to a gap (i-1 or i+1 has a gap in any row).
		// Gap chars ('-') are case-insensitive so raw string access is safe here.
		bool adjacent_gap = false;
		if (i > 0) {
			adjacent_gap =
			    is_gap(star.query_row[i - 1]) || is_gap(star.parent_a_row[i - 1]) || is_gap(star.parent_b_row[i - 1]);
		}
		if (!adjacent_gap && i + 1 < len) {
			adjacent_gap =
			    is_gap(star.query_row[i + 1]) || is_gap(star.parent_a_row[i + 1]) || is_gap(star.parent_b_row[i + 1]);
		}
		if (adjacent_gap) {
			star.diffs[i] = DiffType::IGNORE;
			continue;
		}

		// Skip if any base is ambiguous (N, R, Y, etc.)
		if (is_ambiguous(q) || is_ambiguous(a) || is_ambiguous(b)) {
			star.diffs[i] = DiffType::IGNORE;
			continue;
		}

		// Classify based on which parent the query matches
		if (a == b) {
			if (q == a) {
				star.diffs[i] = DiffType::IGNORE; // All identical
			} else {
				star.diffs[i] = DiffType::NO_VOTE; // Both parents agree, query differs
			}
		} else {
			// Parents differ
			if (q == a) {
				star.diffs[i] = DiffType::MATCH_A;
			} else if (q == b) {
				star.diffs[i] = DiffType::MATCH_B;
			} else {
				star.diffs[i] = DiffType::ABSTAIN;
			}
		}
	}
}

} // namespace miint
