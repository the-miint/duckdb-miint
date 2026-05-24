#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

// One per-base row emitted by PileupWalker::Walk. `ref_pos` is 1-based
// (matches SAM POS). For insertion rows, `ref_pos` is the preceding
// reference position (SAM convention) and `insert_pos` is 1-based within
// the insertion event (0 for reference-aligned rows). For a leading
// insertion (before any ref-consuming op), `ref_pos` = align_start - 1,
// which can be 0. This does not collide with the SAM unmapped sentinel
// because unmapped alignments (CIGAR '*') produce zero rows.
struct PileupRow {
	std::string ref_id;
	std::int64_t ref_pos;
	std::string read_id;
	char ref_base;
	char query_base;
	std::uint8_t query_qual;
	bool query_is_null;
	bool qual_is_null;
	bool ref_base_is_null;
	std::int32_t insert_pos;
};

// Walks a single alignment's CIGAR against a reference sequence and emits
// per-base PileupRow entries by appending to `out`.
//
//   M/=/X : emit one row per ref pos; advance both cursors
//   D/N   : emit one row per ref pos with query_base/query_qual NULL
//   I     : emit one row per inserted base with ref_base NULL,
//           ref_pos = preceding ref position, insert_pos = 1..N
//   S     : advance query only (clip drops from output)
//   H/P   : no-op
class PileupWalker {
public:
	static void Walk(const std::string &cigar, const std::string &ref_id, const std::string &ref_seq,
	                 std::int64_t align_start_1based, const std::string &read_id, const std::string &seq,
	                 const std::uint8_t *qual_data, std::size_t qual_len, bool qual_is_null,
	                 std::vector<PileupRow> &out);
};

} // namespace miint
