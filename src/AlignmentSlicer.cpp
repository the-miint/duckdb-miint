#include "AlignmentSlicer.hpp"
#include "alignment_functions_internal.hpp"

#include <algorithm>
#include <sstream>

namespace miint {

AlignmentSlicer::AlignmentSlicer(int64_t region_start, int64_t region_stop, bool include_deletions)
    : region_start(region_start), region_stop(region_stop), include_deletions(include_deletions) {
	if (region_start < 1) {
		throw InvalidInputException("alignment_slice: region start must be >= 1, got " + std::to_string(region_start));
	}
	if (region_start > region_stop) {
		throw InvalidInputException("alignment_slice: region start (" + std::to_string(region_start) +
		                            ") must be <= region stop (" + std::to_string(region_stop) + ")");
	}
}

SliceResult AlignmentSlicer::Slice(const SliceInput &input) const {
	SliceResult result;
	result.overlaps = false;
	result.position = 0;
	result.stop_position = 0;

	// Guard against invalid position (0 means unmapped in SAM)
	if (input.position < 1) {
		return result;
	}

	// Parse CIGAR
	auto ops = ParseCigarOperations(input.cigar);
	if (ops.empty()) {
		return result; // Unmapped or empty CIGAR
	}

	// Empty region — nothing can overlap
	if (region_start == region_stop) {
		return result;
	}

	// Walk CIGAR operations, tracking reference position and query position
	int64_t ref_pos = input.position; // 1-based current reference position
	int64_t left_query_trim = 0;      // query bases trimmed from left
	int64_t right_query_trim = 0;     // query bases trimmed from right
	int64_t left_hard_clip = 0;       // existing H clips on left
	int64_t right_hard_clip = 0;      // existing H clips on right

	bool has_query_overlap = false; // M/=/X/I within region
	bool has_del_overlap = false;   // D within region

	// Phase tracking: BEFORE we've entered the region, WITHIN, AFTER we've passed it
	enum Phase { BEFORE, WITHIN, AFTER };
	Phase phase = BEFORE;

	// Leading S ops are deferred until the first ref-consuming op resolves whether
	// left trimming is needed. If no M/=/X trimming happens on the left, leading S
	// is preserved; otherwise it becomes part of the left hard clip.
	std::vector<CigarOperation> pending_leading_s;
	int64_t pending_leading_s_query = 0;
	bool leading_s_resolved = false;

	// Collect trimmed CIGAR ops (excluding H clips, which we track separately)
	std::vector<CigarOperation> trimmed_ops;

	for (const auto &op : ops) {
		// Determine how many reference and query bases this op consumes
		int64_t ref_consumed = 0;
		int64_t query_consumed = 0;

		switch (op.op) {
		case 'M':
		case '=':
		case 'X':
			ref_consumed = op.length;
			query_consumed = op.length;
			break;
		case 'I':
			query_consumed = op.length;
			break;
		case 'D':
		case 'N':
			ref_consumed = op.length;
			break;
		case 'S':
			query_consumed = op.length;
			break;
		case 'H':
			// Track existing hard clips
			if (!leading_s_resolved && trimmed_ops.empty() && phase == BEFORE) {
				left_hard_clip += op.length;
			} else {
				right_hard_clip += op.length;
			}
			continue; // Don't advance ref_pos or query_pos
		case 'P':
			// Padding: doesn't consume ref or query
			// Keep it if we're in the WITHIN phase
			if (phase == WITHIN) {
				trimmed_ops.push_back(op);
			}
			continue;
		default:
			continue;
		}

		// Handle S operations (no reference consumption)
		if (op.op == 'S') {
			if (!leading_s_resolved) {
				// Defer leading S decision
				pending_leading_s.push_back(op);
				pending_leading_s_query += op.length;
			} else if (phase == AFTER) {
				right_query_trim += op.length;
			} else {
				// WITHIN: keep soft clip
				trimmed_ops.push_back(op);
			}
			continue;
		}

		// Resolve pending leading S when we hit the first non-S, non-H, non-P op.
		// If the first ref-consuming op starts at or after region_start, no left
		// trimming is needed and leading S is kept. Otherwise, leading S is trimmed.
		// We do NOT set phase here — the actual op handler below determines phase
		// based on whether the op overlaps the region.
		if (!leading_s_resolved) {
			leading_s_resolved = true;
			if (ref_pos >= region_start) {
				// No left trimming — keep leading S
				for (const auto &s : pending_leading_s) {
					trimmed_ops.push_back(s);
				}
			} else {
				// Left trimming will happen — trim leading S
				left_query_trim += pending_leading_s_query;
			}
		}

		// Handle I operations (no reference consumption)
		if (op.op == 'I') {
			if (ref_pos >= region_start && ref_pos < region_stop) {
				// Insertion is at a ref position within the region
				trimmed_ops.push_back(op);
				has_query_overlap = true;
				phase = WITHIN;
			} else if (ref_pos < region_start) {
				left_query_trim += op.length;
			} else {
				right_query_trim += op.length;
				if (phase == WITHIN) {
					phase = AFTER;
				}
			}
			continue;
		}

		// Reference-consuming operations: M, =, X, D, N
		int64_t op_ref_start = ref_pos;
		int64_t op_ref_end = ref_pos + ref_consumed; // exclusive

		if (op_ref_end <= region_start) {
			// Entirely before region
			left_query_trim += query_consumed;
			// phase stays BEFORE
		} else if (op_ref_start >= region_stop) {
			// Entirely after region
			right_query_trim += query_consumed;
			if (phase == WITHIN) {
				phase = AFTER;
			}
		} else {
			// Overlaps region — may need splitting
			int64_t bases_before = std::max(int64_t(0), region_start - op_ref_start);
			int64_t bases_after = std::max(int64_t(0), op_ref_end - region_stop);
			int64_t bases_within = ref_consumed - bases_before - bases_after;

			if (bases_before > 0 && query_consumed > 0) {
				// For M/=/X: query and ref consumed 1:1
				left_query_trim += bases_before;
			}

			if (bases_within > 0) {
				trimmed_ops.push_back({op.op, bases_within});
				if (op.op == 'M' || op.op == '=' || op.op == 'X') {
					has_query_overlap = true;
				} else if (op.op == 'D') {
					has_del_overlap = true;
				}
				// N within region is kept in CIGAR but doesn't count as overlap
				phase = WITHIN;
			}

			if (bases_after > 0) {
				if (query_consumed > 0) {
					right_query_trim += bases_after;
				}
				if (phase == WITHIN) {
					phase = AFTER;
				}
			}
		}

		ref_pos += ref_consumed;
	}

	// Determine overlap
	result.overlaps = has_query_overlap || (include_deletions && has_del_overlap);
	if (!result.overlaps) {
		return result;
	}

	// Build final CIGAR string
	int64_t total_left_H = left_hard_clip + left_query_trim;
	int64_t total_right_H = right_hard_clip + right_query_trim;

	std::vector<CigarOperation> final_ops;
	if (total_left_H > 0) {
		final_ops.push_back({'H', total_left_H});
	}
	final_ops.insert(final_ops.end(), trimmed_ops.begin(), trimmed_ops.end());
	if (total_right_H > 0) {
		final_ops.push_back({'H', total_right_H});
	}

	// Merge adjacent same-type operations
	std::vector<CigarOperation> merged;
	for (const auto &op : final_ops) {
		if (!merged.empty() && merged.back().op == op.op) {
			merged.back().length += op.length;
		} else {
			merged.push_back(op);
		}
	}

	// Convert to CIGAR string
	std::ostringstream cigar_ss;
	for (const auto &op : merged) {
		cigar_ss << op.length << op.op;
	}
	result.cigar = cigar_ss.str();

	// Compute adjusted positions.
	// position: if the alignment starts before the region, it was trimmed to region_start.
	// stop_position: if the alignment extends past region_stop, it was trimmed to region_stop.
	// When no trimming occurs on a side, the original value is preserved (and is already
	// within the region, since we only get here for overlapping alignments).
	result.position = std::max(input.position, region_start);
	result.stop_position = std::min(input.stop_position, region_stop);

	// Trim sequence and quality
	if (!input.sequence.empty()) {
		int64_t seq_len = static_cast<int64_t>(input.sequence.size());
		int64_t start = left_query_trim;
		int64_t end = seq_len - right_query_trim;
		if (start < end && start >= 0 && end <= seq_len) {
			result.sequence = input.sequence.substr(start, end - start);
		}
	}

	if (!input.quality.empty()) {
		int64_t qual_len = static_cast<int64_t>(input.quality.size());
		int64_t start = left_query_trim;
		int64_t end = qual_len - right_query_trim;
		if (start < end && start >= 0 && end <= qual_len) {
			result.quality = std::vector<uint8_t>(input.quality.begin() + start, input.quality.begin() + end);
		}
	}

	return result;
}

} // namespace miint
