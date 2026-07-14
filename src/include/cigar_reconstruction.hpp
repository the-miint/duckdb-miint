#pragma once

#include <stdexcept>
#include <string>

namespace miint {

// Parse a run-length CIGAR string, invoking `fn(count, op)` for each operation in order.
// Validates the grammar: a run of digits forms a positive length; a bare operation
// character (no leading digits) means length 1. Shared by reconstruct_aligned_from_cigar
// and eqx_split_cigar so the parse (and its error handling) lives in exactly one place.
template <typename Fn>
inline void for_each_cigar_op(const std::string &cigar, Fn &&fn) {
	size_t ci = 0;
	while (ci < cigar.size()) {
		size_t num_start = ci;
		while (ci < cigar.size() && cigar[ci] >= '0' && cigar[ci] <= '9') {
			ci++;
		}
		if (ci >= cigar.size()) {
			throw std::runtime_error("CIGAR string ends with digits but no operation character");
		}
		int count = 1;
		if (ci > num_start) {
			try {
				count = std::stoi(cigar.substr(num_start, ci - num_start));
			} catch (const std::out_of_range &) {
				throw std::runtime_error("CIGAR operation length overflows integer");
			}
			if (count <= 0) {
				throw std::runtime_error("CIGAR operation length must be positive");
			}
		}
		char op = cigar[ci++];
		fn(count, op);
	}
}

// Reconstruct gapped alignment strings from a CIGAR string and the original sequences.
// Accepts the union of CIGAR operations used by WFA2 (extended: '=' / 'X') and KSW2
// (standard: 'M' for match-or-mismatch; 'N' for intron skip emitted by ksw_exts2_sse).
// 'I' = insertion in query (gap in subject), 'D' = deletion from query (gap in subject).
// 'N' is treated the same as 'D' for the purposes of alignment display: subject base shown,
// query rendered as '-'. The op-name distinction (intron vs deletion) lives in the CIGAR;
// the reconstructed strings collapse them visually.
inline void reconstruct_aligned_from_cigar(const std::string &query, const std::string &subject,
                                           const std::string &cigar, std::string &query_aligned,
                                           std::string &subject_aligned) {
	query_aligned.clear();
	subject_aligned.clear();
	query_aligned.reserve(query.size() + subject.size());
	subject_aligned.reserve(query.size() + subject.size());

	size_t qi = 0;
	size_t si = 0;

	for_each_cigar_op(cigar, [&](int count, char op) {
		for (int k = 0; k < count; k++) {
			switch (op) {
			case '=':
			case 'X':
			case 'M':
				if (qi >= query.size() || si >= subject.size()) {
					throw std::runtime_error("CIGAR consumes more bases than available in sequences");
				}
				query_aligned += query[qi++];
				subject_aligned += subject[si++];
				break;
			case 'I':
				if (qi >= query.size()) {
					throw std::runtime_error("CIGAR consumes more query bases than available");
				}
				query_aligned += query[qi++];
				subject_aligned += '-';
				break;
			case 'D':
			case 'N':
				if (si >= subject.size()) {
					throw std::runtime_error("CIGAR consumes more subject bases than available");
				}
				query_aligned += '-';
				subject_aligned += subject[si++];
				break;
			default:
				throw std::runtime_error(std::string("Unknown CIGAR operation: ") + op);
			}
		}
	});

	if (qi != query.size()) {
		throw std::runtime_error("CIGAR did not consume all query bases: consumed " + std::to_string(qi) + " of " +
		                         std::to_string(query.size()));
	}
	if (si != subject.size()) {
		throw std::runtime_error("CIGAR did not consume all subject bases: consumed " + std::to_string(si) + " of " +
		                         std::to_string(subject.size()));
	}
}

// ASCII-uppercase one character (locale-independent; safe for signed char). Used so base
// comparison in eqx_split_cigar is case-insensitive.
inline char cigar_ascii_upper(char c) {
	return (c >= 'a' && c <= 'z') ? static_cast<char>(c - ('a' - 'A')) : c;
}

// Rewrite the 'M' (match-or-mismatch) operations of a CIGAR into extended '=' (match) /
// 'X' (mismatch) operations by comparing the aligned bases of the original sequences. This
// is the "eqx" post-pass, matching the repo-wide convention (minimap2's MM_F_EQX flag, the
// `eqx=true` default on the align_minimap2*/save_minimap2_index table functions, and WFA2's
// native extended CIGAR). It exists because the pairwise KSW2 functions call ksw_extz2_sse
// directly and so never reach minimap2's internal eqx pass. Splitting here lets sequence
// identity be read straight off a KSW2 CIGAR (e.g. cigar_sequence_identity), which a bare
// 'M'-only CIGAR cannot support. 'I', 'D', and 'N' operations pass through unchanged.
//
// Base comparison is case-insensitive (as minimap2's own eqx effectively is, working on
// case-folded encoded bases), so soft-masked (lowercase) bases are matches, not mismatches;
// two distinct letters (including different IUPAC ambiguity codes) are mismatches. '=' / 'X'
// inputs are re-derived from the bases (idempotent for honest input); KSW2 never emits them.
//
// Unlike reconstruct_aligned_from_cigar, this does NOT require the CIGAR to consume every base
// of both sequences: a z-drop-truncated alignment yields a CIGAR shorter than the inputs, and
// we faithfully re-encode whatever ops are present (the inner bounds checks still reject a
// CIGAR that consumes MORE bases than exist).
inline std::string eqx_split_cigar(const std::string &cigar, const std::string &query, const std::string &subject) {
	std::string out;
	out.reserve(cigar.size() + cigar.size() / 2);

	size_t qi = 0;
	size_t si = 0;

	auto append_run = [&out](int count, char op) {
		if (count > 0) {
			out += std::to_string(count);
			out += op;
		}
	};

	for_each_cigar_op(cigar, [&](int count, char op) {
		switch (op) {
		case 'M':
		case '=':
		case 'X': {
			// Walk the aligned columns, coalescing consecutive same-result columns into
			// a single '=' or 'X' run.
			int run = 0;
			bool run_is_match = false;
			for (int k = 0; k < count; k++) {
				if (qi >= query.size() || si >= subject.size()) {
					throw std::runtime_error("CIGAR consumes more bases than available in sequences");
				}
				bool is_match = cigar_ascii_upper(query[qi]) == cigar_ascii_upper(subject[si]);
				qi++;
				si++;
				if (run > 0 && is_match != run_is_match) {
					append_run(run, run_is_match ? '=' : 'X');
					run = 0;
				}
				run_is_match = is_match;
				run++;
			}
			append_run(run, run_is_match ? '=' : 'X');
			break;
		}
		case 'I':
			if (qi + static_cast<size_t>(count) > query.size()) {
				throw std::runtime_error("CIGAR consumes more query bases than available");
			}
			qi += count;
			append_run(count, 'I');
			break;
		case 'D':
		case 'N':
			if (si + static_cast<size_t>(count) > subject.size()) {
				throw std::runtime_error("CIGAR consumes more subject bases than available");
			}
			si += count;
			append_run(count, op);
			break;
		default:
			throw std::runtime_error(std::string("Unknown CIGAR operation: ") + op);
		}
	});

	return out;
}

} // namespace miint
