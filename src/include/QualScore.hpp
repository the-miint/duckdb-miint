#pragma once
#include <string>
#include <cstdint>
#include <vector>
#include <span>
#include <stdexcept>
#include <cstdint>

namespace miint {
//! Object to represent Qual store data independent of offset
class QualScore {
private:
	//! Quality in ascii
	std::string qual_;

public:
	//! Constructor from string: store as-is
	explicit QualScore(const std::string &qual_str) noexcept : qual_(qual_str) {
	}

	//! Constructor from uint8 vector: interpret as quality scores, convert to characters
	explicit QualScore(const std::vector<uint8_t> &qual_vec, int offset = 33) {
		qual_.reserve(qual_vec.size());
		for (uint8_t q : qual_vec) {
			int offset_q = q + offset;

			if (offset_q > 126) {
				throw std::invalid_argument("Invalid quality score: " + std::to_string(q));
			}
			qual_.push_back(static_cast<char>(offset_q));
		}
	}

	//! Return as string (raw stored characters)
	const std::string &as_string() const noexcept {
		return qual_;
	}

	//! Number of quality scores
	size_t size() const noexcept {
		return qual_.size();
	}

	//! Write decoded quality scores directly to caller-supplied buffer (avoids heap allocation).
	//! Returns number of bytes written.
	size_t write_decoded(uint8_t *dest, int offset = 33) const {
		size_t n = qual_.size();
		for (size_t i = 0; i < n; i++) {
			int val = static_cast<int>(qual_[i]) - offset;
			if (val < 0 || val > 93) {
				throw std::runtime_error("Stored quality character out of range: " + std::string(1, qual_[i]));
			}
			dest[i] = static_cast<uint8_t>(val);
		}
		return n;
	}

	//! Return as vector<uint8_t> (quality scores: stored char minus offset)
	std::vector<uint8_t> as_vec(int offset = 33) const {
		std::vector<uint8_t> out;
		out.reserve(qual_.size());
		for (char c : qual_) {
			int val = static_cast<int>(c) - offset;
			if (val < 0 || val > 93) {
				throw std::runtime_error("Stored quality character out of range: " + std::string(1, c));
			}
			out.push_back(static_cast<uint8_t>(val));
		}
		return out;
	}
};

//! phred_scores: already decoded (ASCII–offset), one entry per base
//! min_quality: minimum acceptable score
//! window_length: number of consecutive low‐quality bases to flag
//! returns: zero‐based index of first low‐quality window, length + 1 if no window
std::size_t find_low_quality_window(const std::span<uint8_t> &phred_scores, uint8_t min_quality,
                                    std::size_t window_length);
} // namespace miint
