// SPDX-License-Identifier: MIT
//
// FASTQ record encoder — emits the four-line @id\nseq\n+\nqual\n form via a
// caller-supplied byte sink. Quality scores are encoded with the configurable
// offset (Phred+33 by default; +64 for legacy data) and validated against the
// printable ASCII range.
//
// The output is a plain `void(const char *, size_t)` callback rather than any
// DuckDB type, so the unit-test target can link this file directly without
// pulling in the duckdb library. Production wires the callback to a
// `MemoryStream::WriteData` call in the COPY path or to a streaming gzip+MD5
// pipeline in `ena_upload_reads`.

#pragma once

#include <cstdint>
#include <cstdlib>
#include <functional>
#include <string>

namespace duckdb {

class FastqEncoder {
public:
	// Bytes written are streamed to the caller via this signature. The encoder
	// does not retain the buffer between calls; safe to write straight into a
	// MemoryStream, gzip filter, or std::string.
	using Sink = std::function<void(const char *data, std::size_t size)>;

	// `qual_offset` is added to each raw quality byte before emission (33 for
	// Phred+33, 64 for legacy Solexa). The encoder rejects any encoded value
	// that exceeds 126 (printable-ASCII upper bound) with std::runtime_error.
	explicit FastqEncoder(uint8_t qual_offset = 33) : qual_offset(qual_offset) {
	}

	// Append one FASTQ record via the sink. `comment` may be empty (or nullptr
	// if `comment_size` is zero); when non-empty it's emitted on the header
	// line separated from `id` by a single space (no trailing space when
	// comment is empty). `quality_data` must be at least `quality_length`
	// bytes; sequence and quality lengths must match (caller enforced).
	void Encode(const Sink &sink, const char *id, std::size_t id_size, const char *comment, std::size_t comment_size,
	            const char *seq, std::size_t seq_size, const std::uint8_t *quality_data, std::size_t quality_length);

	uint8_t GetQualOffset() const {
		return qual_offset;
	}

private:
	uint8_t qual_offset;
	// Reused across records to amortize the per-call allocation; resized to
	// `quality_length` on each Encode call.
	std::string qual_encoded_buf;
};

} // namespace duckdb
