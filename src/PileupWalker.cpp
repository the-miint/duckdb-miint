#include "PileupWalker.hpp"
#include "alignment_functions_internal.hpp"

namespace miint {

static constexpr const char *FN_NAME = "compute_pileup";

void PileupWalker::Walk(const std::string &cigar, const std::string &ref_id, const std::string &ref_seq,
                        std::int64_t align_start_1based, const std::string &read_id, const std::string &seq,
                        const std::uint8_t *qual_data, std::size_t qual_len, bool qual_is_null,
                        std::vector<PileupRow> &out) {
	auto ops = ParseCigarOperations(cigar);
	if (ops.empty()) {
		return;
	}

	// Validate: sum of query-consuming ops (M=XIS) must equal seq.size().
	// Hard clips (H) don't consume query because the SAM seq is already
	// pre-trimmed per spec.
	std::int64_t expected_query_len = 0;
	for (const auto &op : ops) {
		switch (op.op) {
		case 'M':
		case '=':
		case 'X':
		case 'I':
		case 'S':
			expected_query_len += op.length;
			break;
		default:
			break;
		}
	}
	if (expected_query_len != static_cast<std::int64_t>(seq.size())) {
		throw InvalidInputException(std::string(FN_NAME) + ": read '" + read_id + "' CIGAR consumes " +
		                            std::to_string(expected_query_len) + " query bases but sequence has " +
		                            std::to_string(seq.size()) + " (CIGAR=" + cigar + ")");
	}
	if (!qual_is_null && qual_len != seq.size()) {
		throw InvalidInputException(std::string(FN_NAME) + ": read '" + read_id + "' sequence length (" +
		                            std::to_string(seq.size()) + ") does not match qual length (" +
		                            std::to_string(qual_len) + ")");
	}

	std::int64_t ref_pos = align_start_1based; // 1-based
	std::size_t query_pos = 0;                 // 0-based offset into seq/qual
	const std::int64_t ref_len = static_cast<std::int64_t>(ref_seq.size());

	auto emit_row = [&](std::int64_t length, bool query_present) {
		for (std::int64_t k = 0; k < length; ++k) {
			if (ref_pos + k - 1 >= ref_len) {
				throw InvalidInputException(std::string(FN_NAME) + ": read '" + read_id + "' alignment at ref_pos " +
				                            std::to_string(ref_pos + k) + " extends beyond reference '" + ref_id +
				                            "' (length " + std::to_string(ref_len) + ")");
			}
			PileupRow r;
			r.ref_id = ref_id;
			r.ref_pos = ref_pos + k;
			r.read_id = read_id;
			r.ref_base = ref_seq[ref_pos + k - 1];
			r.ref_base_is_null = false;
			r.insert_pos = 0;
			if (query_present) {
				r.query_base = seq[query_pos + k];
				r.query_qual = qual_is_null ? 0 : qual_data[query_pos + k];
				r.query_is_null = false;
				r.qual_is_null = qual_is_null;
			} else {
				r.query_base = 0;
				r.query_qual = 0;
				r.query_is_null = true;
				r.qual_is_null = true;
			}
			out.push_back(std::move(r));
		}
	};

	for (const auto &op : ops) {
		switch (op.op) {
		case 'M':
		case '=':
		case 'X':
			emit_row(op.length, /*query_present=*/true);
			ref_pos += op.length;
			query_pos += op.length;
			break;
		case 'I': {
			std::int64_t ins_ref_pos = ref_pos - 1;
			for (std::int64_t k = 0; k < op.length; ++k) {
				PileupRow r;
				r.ref_id = ref_id;
				r.ref_pos = ins_ref_pos;
				r.read_id = read_id;
				r.ref_base = 0;
				r.ref_base_is_null = true;
				r.insert_pos = static_cast<std::int32_t>(k + 1);
				r.query_base = seq[query_pos + k];
				r.query_qual = qual_is_null ? 0 : qual_data[query_pos + k];
				r.query_is_null = false;
				r.qual_is_null = qual_is_null;
				out.push_back(std::move(r));
			}
			query_pos += op.length;
			break;
		}
		case 'D':
		case 'N':
			emit_row(op.length, /*query_present=*/false);
			ref_pos += op.length;
			break;
		case 'S':
			query_pos += op.length;
			break;
		case 'H':
		case 'P':
			break;
		default:
			throw InvalidInputException(std::string(FN_NAME) + ": unsupported CIGAR op '" + std::string(1, op.op) +
			                            "' in read '" + read_id + "'");
		}
	}
}

} // namespace miint
