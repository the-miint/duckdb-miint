#pragma once

#include "SpectrumBatch.hpp"
#include <cstdint>
#include <memory>
#include <string>

namespace miint {

class MzXMLReader {
public:
	explicit MzXMLReader(const std::string &path);
	~MzXMLReader();

	// Non-copyable, moveable
	MzXMLReader(const MzXMLReader &) = delete;
	MzXMLReader &operator=(const MzXMLReader &) = delete;
	MzXMLReader(MzXMLReader &&) noexcept;
	MzXMLReader &operator=(MzXMLReader &&) noexcept;

	// Read up to n spectra. Returns batch with 0..n spectra (0 = end of file).
	[[nodiscard]] MzMLSpectrumBatch read_spectra(size_t n);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace miint
