#include "CoverageDepthCalculator.hpp"
#include "alignment_functions_internal.hpp"

namespace miint {

// Memory note: allocates reference_length * 4 bytes. For a human chromosome (~250M positions)
// that's ~1GB per group. Callers should be aware of this when using with large references.
CoverageDepthCalculator::CoverageDepthCalculator(int64_t reference_length, bool include_deletions)
    : depths(static_cast<size_t>(reference_length), 0), reference_length(reference_length),
      include_deletions(include_deletions), has_reads(false) {
}

// Increment depths[start..end) with bounds clamping to [0, reference_length).
// start and end are 0-based indices.
void CoverageDepthCalculator::IncrementRange(int64_t start, int64_t end) {
	if (start < 0) {
		start = 0;
	}
	if (end < 0) {
		end = 0;
	}
	if (end > reference_length) {
		end = reference_length;
	}
	if (start >= end) {
		return; // Empty range or entirely out of bounds, nothing to do
	}
	// After clamping above, 0 <= start < end <= reference_length is guaranteed

	for (int64_t i = start; i < end; i++) {
		auto idx = static_cast<size_t>(i);
		if (depths[idx] == UINT32_MAX) {
			throw InvalidInputException("Coverage depth overflow at position " + std::to_string(i + 1) +
			                            ": depth exceeds maximum (" + std::to_string(UINT32_MAX) + ")");
		}
		depths[idx]++;
	}
}

void CoverageDepthCalculator::AddRead(int64_t position, int64_t stop_position, const std::string &cigar) {
	auto ops = ParseCigarOperations(cigar);
	if (ops.empty()) {
		return; // Unmapped or empty CIGAR
	}

	has_reads = true;

	// Check if we need CIGAR walking: only when there are N ops, or D ops in exclude_deletions mode
	bool needs_walking = false;
	for (const auto &op : ops) {
		if (op.op == 'N' || (op.op == 'D' && !include_deletions)) {
			needs_walking = true;
			break;
		}
	}

	if (!needs_walking) {
		// Fast path: just increment [position, stop_position) — all ref-consuming ops count
		IncrementRange(position - 1, stop_position - 1);
		return;
	}

	// Slow path: walk CIGAR ops to handle N and/or D exclusion
	int64_t ref_pos = position - 1; // 0-based index into depths

	for (const auto &op : ops) {
		switch (op.op) {
		case 'M':
		case '=':
		case 'X':
			// Always count as coverage
			IncrementRange(ref_pos, ref_pos + op.length);
			ref_pos += op.length;
			break;
		case 'D':
			// Only count in include_deletions mode (we're in the slow path
			// because of N ops elsewhere in this CIGAR)
			if (include_deletions) {
				IncrementRange(ref_pos, ref_pos + op.length);
			}
			ref_pos += op.length;
			break;
		case 'N':
			// Never counts as coverage — skip reference positions
			ref_pos += op.length;
			break;
		case 'I':
		case 'S':
		case 'P':
		case 'H':
			// These don't consume reference positions
			break;
		default:
			break;
		}
	}
}

// Element-wise addition for parallel GROUP BY merging.
void CoverageDepthCalculator::Combine(const CoverageDepthCalculator &other) {
	for (size_t i = 0; i < depths.size() && i < other.depths.size(); i++) {
		if (depths[i] > UINT32_MAX - other.depths[i]) {
			throw InvalidInputException(
			    "Coverage depth overflow at position " + std::to_string(static_cast<int64_t>(i) + 1) +
			    " during parallel merge: depth exceeds maximum (" + std::to_string(UINT32_MAX) + ")");
		}
		depths[i] += other.depths[i];
	}
	if (other.has_reads) {
		has_reads = true;
	}
}

const std::vector<uint32_t> &CoverageDepthCalculator::GetDepths() const {
	return depths;
}

bool CoverageDepthCalculator::Empty() const {
	return !has_reads;
}

} // namespace miint
