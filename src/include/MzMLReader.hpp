#pragma once

#include "SpectrumBatch.hpp"
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// SOA batch for chromatogram data
struct MzMLChromatogramBatch {
	std::vector<int32_t> chromatogram_index;
	std::vector<std::string> chromatogram_id;
	std::vector<std::string> chromatogram_type; // "TIC" / "SIC" / "SRM" / "BPC" / ""
	std::vector<double> precursor_mz;
	std::vector<bool> precursor_mz_valid;
	std::vector<double> product_mz;
	std::vector<bool> product_mz_valid;
	std::vector<std::vector<double>> time_array;
	std::vector<std::vector<double>> intensity_array;

	[[nodiscard]] size_t size() const {
		return chromatogram_index.size();
	}
	[[nodiscard]] bool empty() const {
		return chromatogram_index.empty();
	}
};

class MzMLReader {
public:
	explicit MzMLReader(const std::string &path);
	~MzMLReader();

	// Non-copyable, moveable
	MzMLReader(const MzMLReader &) = delete;
	MzMLReader &operator=(const MzMLReader &) = delete;
	MzMLReader(MzMLReader &&) noexcept;
	MzMLReader &operator=(MzMLReader &&) noexcept;

	// Read up to n spectra. Returns batch with 0..n spectra (0 = end of spectra).
	[[nodiscard]] MzMLSpectrumBatch read_spectra(size_t n);

	// Read up to n chromatograms. Returns batch with 0..n chromatograms (0 = end).
	// Must be called after all spectra have been consumed (or after spectra are exhausted).
	[[nodiscard]] MzMLChromatogramBatch read_chromatograms(size_t n);

	[[nodiscard]] bool has_spectra() const;
	[[nodiscard]] bool has_chromatograms() const;

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace miint
