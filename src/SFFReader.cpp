// SFF (Standard Flowgram Format) binary parser.
// Format specification by 454 Life Sciences/Roche, Whitehead Institute,
// and Wellcome Trust Sanger Institute.
// Format details referenced from Biopython's Bio.SeqIO.SffIO module
// (https://github.com/biopython/biopython), licensed under the
// Biopython License Agreement. See THIRD_PARTY_LICENSES.md.

#include "SFFReader.hpp"
#include <algorithm>
#include <stdexcept>

namespace miint {

// SFF magic number: ".sff" as big-endian uint32
static constexpr uint32_t SFF_MAGIC = 0x2E736666;
// SFF version 1
static constexpr uint8_t SFF_VERSION[4] = {0, 0, 0, 1};
// Only supported flowgram format code
static constexpr uint8_t SFF_FLOWGRAM_FORMAT = 1;
// Minimum size of the fixed portion of the file header
static constexpr size_t SFF_HEADER_FIXED_SIZE = 31;
// Per-read header: 16 bytes fixed + name + padding to 8-byte boundary
static constexpr size_t SFF_READ_HEADER_FIXED_SIZE = 16;
// Per-read data sizes
static constexpr size_t SFF_FLOWGRAM_VALUE_SIZE = 2; // uint16 per flow
static constexpr size_t SFF_PER_BASE_FIELDS = 3;     // flow_index + base + quality

static size_t pad_to_8(size_t len) {
	size_t remainder = len % 8;
	return remainder == 0 ? len : len + (8 - remainder);
}

SFFReader::SFFReader(const std::string &path, bool trim) : path_(path), trim_(trim) {
	file_.open(path, std::ios::binary);
	if (!file_.is_open()) {
		throw std::runtime_error("Failed to open SFF file: " + path);
	}

	// Determine file size
	file_.seekg(0, std::ios::end);
	auto file_size = file_.tellg();
	file_.seekg(0, std::ios::beg);

	if (file_size < static_cast<std::streamoff>(SFF_HEADER_FIXED_SIZE)) {
		throw std::runtime_error("SFF file too small for header (" + std::to_string(file_size) + " bytes): " + path);
	}

	// Read the 31-byte fixed header
	uint8_t fixed_header[SFF_HEADER_FIXED_SIZE];
	file_.read(reinterpret_cast<char *>(fixed_header), SFF_HEADER_FIXED_SIZE);
	if (!file_) {
		throw std::runtime_error("Failed to read SFF file header: " + path);
	}

	// Validate magic number
	uint32_t magic = read_uint32_be(fixed_header);
	if (magic != SFF_MAGIC) {
		throw std::runtime_error("Invalid SFF magic number (expected 0x2E736666, got 0x" + ([&]() {
			                         char hex[9];
			                         snprintf(hex, sizeof(hex), "%08X", magic);
			                         return std::string(hex);
		                         })() +
		                         "): " + path);
	}

	// Validate version
	if (fixed_header[4] != SFF_VERSION[0] || fixed_header[5] != SFF_VERSION[1] || fixed_header[6] != SFF_VERSION[2] ||
	    fixed_header[7] != SFF_VERSION[3]) {
		throw std::runtime_error("Unsupported SFF version (expected 0.0.0.1, got " + std::to_string(fixed_header[4]) +
		                         "." + std::to_string(fixed_header[5]) + "." + std::to_string(fixed_header[6]) + "." +
		                         std::to_string(fixed_header[7]) + "): " + path);
	}

	// Parse header fields
	index_offset_ = read_uint64_be(fixed_header + 8);
	index_length_ = read_uint32_be(fixed_header + 16);
	number_of_reads_ = read_uint32_be(fixed_header + 20);
	uint16_t header_length = read_uint16_be(fixed_header + 24);
	uint16_t key_length = read_uint16_be(fixed_header + 26);
	number_of_flows_ = read_uint16_be(fixed_header + 28);
	uint8_t flowgram_format = fixed_header[30];

	// Validate flowgram format
	if (flowgram_format != SFF_FLOWGRAM_FORMAT) {
		throw std::runtime_error("Unsupported SFF flowgram format code (expected 1, got " +
		                         std::to_string(flowgram_format) + "): " + path);
	}

	// Read flow_chars directly into string (#6: avoid temp vector copy)
	flow_chars_.resize(number_of_flows_);
	file_.read(flow_chars_.data(), number_of_flows_);
	if (!file_) {
		throw std::runtime_error("Failed to read flow_chars from SFF file: " + path);
	}

	// Read key_sequence directly into string (#6: avoid temp vector copy)
	key_sequence_.resize(key_length);
	file_.read(key_sequence_.data(), key_length);
	if (!file_) {
		throw std::runtime_error("Failed to read key_sequence from SFF file: " + path);
	}

	// #8: precompute padded index length once
	padded_index_length_ = (index_length_ > 0) ? pad_to_8(index_length_) : 0;

	// Seek to end of header (padded to 8-byte boundary)
	file_.seekg(header_length, std::ios::beg);
}

SequenceRecordBatch SFFReader::read(int n) {
	SequenceRecordBatch batch(/*paired=*/false);
	if (n <= 0 || reads_parsed_ >= number_of_reads_) {
		return batch;
	}

	uint32_t remaining = number_of_reads_ - reads_parsed_;
	uint32_t to_read = std::min(static_cast<uint32_t>(n), remaining);
	batch.reserve(to_read);

	for (uint32_t i = 0; i < to_read; i++) {
		skip_index_if_needed();
		auto rh = parse_read_header();
		parse_read_data(rh, batch);
		reads_parsed_++;
	}

	return batch;
}

SFFReader::ReadHeader SFFReader::parse_read_header() {
	// #5: cache position before read for error messages
	auto header_start = file_.tellg();

	uint8_t buf[SFF_READ_HEADER_FIXED_SIZE];
	file_.read(reinterpret_cast<char *>(buf), SFF_READ_HEADER_FIXED_SIZE);
	if (!file_) {
		throw std::runtime_error("SFF file truncated reading read header at position " + std::to_string(header_start) +
		                         ": " + path_);
	}

	ReadHeader rh;
	uint16_t read_header_length = read_uint16_be(buf);
	uint16_t name_length = read_uint16_be(buf + 2);
	rh.seq_len = read_uint32_be(buf + 4);
	rh.clip_qual_left = read_uint16_be(buf + 8);
	rh.clip_qual_right = read_uint16_be(buf + 10);
	rh.clip_adapter_left = read_uint16_be(buf + 12);
	rh.clip_adapter_right = read_uint16_be(buf + 14);

	// #2: bounds check seq_len
	if (rh.seq_len > SFF_MAX_SEQ_LEN) {
		throw std::runtime_error("SFF read at position " + std::to_string(header_start) +
		                         " has seq_len=" + std::to_string(rh.seq_len) + " which exceeds maximum (" +
		                         std::to_string(SFF_MAX_SEQ_LEN) + "), file is likely corrupted: " + path_);
	}

	// #7: validate clip coordinates against seq_len
	if (rh.clip_qual_left > rh.seq_len) {
		throw std::runtime_error("SFF read at position " + std::to_string(header_start) +
		                         " has clip_qual_left=" + std::to_string(rh.clip_qual_left) +
		                         " > seq_len=" + std::to_string(rh.seq_len) + ": " + path_);
	}
	if (rh.clip_qual_right > rh.seq_len) {
		throw std::runtime_error("SFF read at position " + std::to_string(header_start) +
		                         " has clip_qual_right=" + std::to_string(rh.clip_qual_right) +
		                         " > seq_len=" + std::to_string(rh.seq_len) + ": " + path_);
	}
	if (rh.clip_adapter_left > rh.seq_len) {
		throw std::runtime_error("SFF read at position " + std::to_string(header_start) +
		                         " has clip_adapter_left=" + std::to_string(rh.clip_adapter_left) +
		                         " > seq_len=" + std::to_string(rh.seq_len) + ": " + path_);
	}
	if (rh.clip_adapter_right > rh.seq_len) {
		throw std::runtime_error("SFF read at position " + std::to_string(header_start) +
		                         " has clip_adapter_right=" + std::to_string(rh.clip_adapter_right) +
		                         " > seq_len=" + std::to_string(rh.seq_len) + ": " + path_);
	}

	// #6: read name directly into string (avoid temp vector copy)
	rh.name.resize(name_length);
	file_.read(rh.name.data(), name_length);
	if (!file_) {
		throw std::runtime_error("SFF file truncated reading read name at position " + std::to_string(header_start) +
		                         ": " + path_);
	}

	// Seek to end of read header (padded to 8-byte boundary)
	file_.seekg(header_start + static_cast<std::streamoff>(read_header_length));

	return rh;
}

void SFFReader::parse_read_data(const ReadHeader &rh, SequenceRecordBatch &batch) {
	// #5: cache position before reads for error messages
	auto data_start = file_.tellg();

	// #9: use named constants for data section layout
	// Skip flowgram values: number_of_flows * SFF_FLOWGRAM_VALUE_SIZE bytes
	file_.seekg(static_cast<std::streamoff>(number_of_flows_) * SFF_FLOWGRAM_VALUE_SIZE, std::ios::cur);

	// Skip flow_index: seq_len bytes
	file_.seekg(rh.seq_len, std::ios::cur);

	// Read bases: seq_len bytes
	std::string bases(rh.seq_len, '\0');
	file_.read(bases.data(), rh.seq_len);
	if (!file_) {
		throw std::runtime_error("SFF file truncated reading sequence data at position " + std::to_string(data_start) +
		                         ": " + path_);
	}

	// Read quality scores: seq_len bytes (raw Phred, 0-93)
	std::vector<uint8_t> raw_quals(rh.seq_len);
	file_.read(reinterpret_cast<char *>(raw_quals.data()), rh.seq_len);
	if (!file_) {
		throw std::runtime_error("SFF file truncated reading quality data at position " + std::to_string(data_start) +
		                         ": " + path_);
	}

	// Seek past padding to 8-byte boundary
	// #9: data section = flowgram + flow_index + bases + quals
	size_t data_raw_len = static_cast<size_t>(number_of_flows_) * SFF_FLOWGRAM_VALUE_SIZE +
	                      static_cast<size_t>(rh.seq_len) * SFF_PER_BASE_FIELDS;
	file_.seekg(data_start + static_cast<std::streamoff>(pad_to_8(data_raw_len)));

	// Apply trimming if enabled
	if (trim_) {
		apply_trim(rh, bases, raw_quals);
	}

	// Store in batch
	batch.read_ids.push_back(std::move(rh.name));
	batch.comments.emplace_back(); // SFF has no comment field
	batch.sequences1.push_back(std::move(bases));
	batch.quals1.emplace_back(raw_quals, 33);
}

void SFFReader::apply_trim(const ReadHeader &rh, std::string &bases, std::vector<uint8_t> &quals) {
	uint32_t seq_len = rh.seq_len;

	// Compute effective trim boundaries per SFF spec:
	// clip value of 0 means "no clip" (use full extent)
	// Left clips are 1-based; convert to 0-based after taking max
	uint32_t left_1based = std::max(rh.clip_qual_left > 0 ? rh.clip_qual_left : static_cast<uint16_t>(1),
	                                rh.clip_adapter_left > 0 ? rh.clip_adapter_left : static_cast<uint16_t>(1));
	uint32_t effective_left = left_1based - 1; // 0-based

	// Right clips are 1-based inclusive
	uint32_t effective_right =
	    std::min(rh.clip_qual_right > 0 ? static_cast<uint32_t>(rh.clip_qual_right) : seq_len,
	             rh.clip_adapter_right > 0 ? static_cast<uint32_t>(rh.clip_adapter_right) : seq_len);

	// Edge case: overlapping clips produce empty result
	if (effective_left >= effective_right) {
		bases.clear();
		quals.clear();
		return;
	}

	uint32_t trim_len = effective_right - effective_left;
	bases = bases.substr(effective_left, trim_len);
	quals.erase(quals.begin() + effective_right, quals.end());
	quals.erase(quals.begin(), quals.begin() + effective_left);
}

void SFFReader::skip_index_if_needed() {
	// #8: early return if no index or already skipped
	if (index_offset_ == 0 || index_skipped_) {
		return;
	}

	// #4: check tellg() for error before casting to unsigned
	auto raw_pos = file_.tellg();
	if (raw_pos < 0) {
		return; // stream error - let the next read() catch it
	}

	if (static_cast<uint64_t>(raw_pos) >= index_offset_) {
		file_.seekg(static_cast<std::streamoff>(index_offset_ + padded_index_length_), std::ios::beg);
		index_skipped_ = true;
	}
}

// SFF is always big-endian per spec, regardless of host endianness
uint16_t SFFReader::read_uint16_be(const uint8_t *buf) {
	return static_cast<uint16_t>(buf[0]) << 8 | static_cast<uint16_t>(buf[1]);
}

uint32_t SFFReader::read_uint32_be(const uint8_t *buf) {
	return static_cast<uint32_t>(buf[0]) << 24 | static_cast<uint32_t>(buf[1]) << 16 |
	       static_cast<uint32_t>(buf[2]) << 8 | static_cast<uint32_t>(buf[3]);
}

uint64_t SFFReader::read_uint64_be(const uint8_t *buf) {
	return static_cast<uint64_t>(buf[0]) << 56 | static_cast<uint64_t>(buf[1]) << 48 |
	       static_cast<uint64_t>(buf[2]) << 40 | static_cast<uint64_t>(buf[3]) << 32 |
	       static_cast<uint64_t>(buf[4]) << 24 | static_cast<uint64_t>(buf[5]) << 16 |
	       static_cast<uint64_t>(buf[6]) << 8 | static_cast<uint64_t>(buf[7]);
}

} // namespace miint
