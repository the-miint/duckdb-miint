// SPDX-License-Identifier: MIT
//
// Implementation of FastqEncoder. See fastq_encoder.hpp.

#include "fastq_encoder.hpp"

#include <stdexcept>
#include <string>

namespace duckdb {

void FastqEncoder::Encode(const Sink &sink, const char *id, std::size_t id_size, const char *comment,
                          std::size_t comment_size, const char *seq, std::size_t seq_size,
                          const std::uint8_t *quality_data, std::size_t quality_length) {
	qual_encoded_buf.resize(quality_length);
	for (std::size_t k = 0; k < quality_length; k++) {
		const int encoded = static_cast<int>(quality_data[k]) + qual_offset;
		if (encoded > 126) {
			throw std::runtime_error("Quality score overflow: " + std::to_string(quality_data[k]) + " + " +
			                         std::to_string(qual_offset) + " = " + std::to_string(encoded) +
			                         " exceeds valid ASCII range (max 126)");
		}
		qual_encoded_buf[k] = static_cast<char>(encoded);
	}

	sink("@", 1);
	sink(id, id_size);
	if (comment_size > 0) {
		sink(" ", 1);
		sink(comment, comment_size);
	}
	sink("\n", 1);
	sink(seq, seq_size);
	sink("\n+\n", 3);
	sink(qual_encoded_buf.data(), quality_length);
	sink("\n", 1);
}

} // namespace duckdb
