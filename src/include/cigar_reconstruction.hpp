#pragma once

#include <stdexcept>
#include <string>

namespace miint {

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
	}

	if (qi != query.size()) {
		throw std::runtime_error("CIGAR did not consume all query bases: consumed " + std::to_string(qi) + " of " +
		                         std::to_string(query.size()));
	}
	if (si != subject.size()) {
		throw std::runtime_error("CIGAR did not consume all subject bases: consumed " + std::to_string(si) + " of " +
		                         std::to_string(subject.size()));
	}
}

} // namespace miint
