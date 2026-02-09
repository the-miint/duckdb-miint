#pragma once

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include "SequenceRecord.hpp"

namespace miint {

// Maximum sequence length we'll accept from an SFF file (10 MB).
// 454 reads are typically <1000 bp; anything beyond this is file corruption.
static constexpr uint32_t SFF_MAX_SEQ_LEN = 10 * 1024 * 1024;

class SFFReader {
public:
	// trim parameter controls whether clip values are applied (Phase 4).
	// No default - callers must be explicit.
	SFFReader(const std::string &path, bool trim);

	// Header accessors
	uint32_t number_of_reads() const {
		return number_of_reads_;
	}
	uint16_t number_of_flows() const {
		return number_of_flows_;
	}
	const std::string &flow_chars() const {
		return flow_chars_;
	}
	const std::string &key_sequence() const {
		return key_sequence_;
	}
	uint64_t index_offset() const {
		return index_offset_;
	}
	uint32_t index_length() const {
		return index_length_;
	}

	// Read up to n records, returning a batch. Returns empty batch at EOF.
	SequenceRecordBatch read(int n);

private:
	std::string path_;
	std::ifstream file_;
	bool trim_;

	// File header fields
	uint32_t number_of_reads_ = 0;
	uint16_t number_of_flows_ = 0;
	std::string flow_chars_;
	std::string key_sequence_;
	uint64_t index_offset_ = 0;
	uint32_t index_length_ = 0;

	// Precomputed from index_length_ (padded to 8-byte boundary)
	size_t padded_index_length_ = 0;

	// Read state
	uint32_t reads_parsed_ = 0;
	bool index_skipped_ = false;

	// Internal parsing helpers
	struct ReadHeader {
		std::string name;
		uint32_t seq_len;
		uint16_t clip_qual_left;
		uint16_t clip_qual_right;
		uint16_t clip_adapter_left;
		uint16_t clip_adapter_right;
	};

	ReadHeader parse_read_header();
	void parse_read_data(const ReadHeader &rh, SequenceRecordBatch &batch);
	void apply_trim(const ReadHeader &rh, std::string &bases, std::vector<uint8_t> &quals);
	void skip_index_if_needed();

	// SFF is always big-endian per spec, regardless of host endianness
	static uint16_t read_uint16_be(const uint8_t *buf);
	static uint32_t read_uint32_be(const uint8_t *buf);
	static uint64_t read_uint64_be(const uint8_t *buf);
};

} // namespace miint
