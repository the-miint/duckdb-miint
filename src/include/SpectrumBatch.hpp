#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

// SOA batch for spectrum data — all vectors are parallel (same index = same spectrum).
// Shared output type for MzMLReader and MzXMLReader.
struct SpectrumBatch {
	std::vector<int32_t> spectrum_index;
	std::vector<std::string> spectrum_id;
	std::vector<int32_t> scan_number;
	std::vector<bool> scan_number_valid;
	std::vector<int32_t> ms_level;
	std::vector<double> retention_time;
	std::vector<bool> retention_time_valid;
	std::vector<std::string> spectrum_type; // "centroid" / "profile" / ""
	std::vector<std::string> polarity;      // "positive" / "negative" / ""
	std::vector<double> base_peak_mz;
	std::vector<bool> base_peak_mz_valid;
	std::vector<double> base_peak_intensity;
	std::vector<bool> base_peak_intensity_valid;
	std::vector<double> total_ion_current;
	std::vector<bool> total_ion_current_valid;
	std::vector<double> lowest_mz;
	std::vector<bool> lowest_mz_valid;
	std::vector<double> highest_mz;
	std::vector<bool> highest_mz_valid;
	std::vector<int32_t> default_array_length;
	std::vector<double> precursor_mz;
	std::vector<bool> precursor_mz_valid;
	std::vector<int32_t> precursor_charge;
	std::vector<bool> precursor_charge_valid;
	std::vector<double> precursor_intensity;
	std::vector<bool> precursor_intensity_valid;
	std::vector<double> isolation_window_target;
	std::vector<bool> isolation_window_target_valid;
	std::vector<double> isolation_window_lower;
	std::vector<bool> isolation_window_lower_valid;
	std::vector<double> isolation_window_upper;
	std::vector<bool> isolation_window_upper_valid;
	std::vector<std::string> activation_method; // "CID" / "HCD" / "ETD" / ""
	std::vector<double> collision_energy;
	std::vector<bool> collision_energy_valid;
	std::vector<std::vector<double>> mz_array;
	std::vector<std::vector<double>> intensity_array;
	std::vector<std::string> filter_string;
	std::vector<double> scan_window_lower;
	std::vector<bool> scan_window_lower_valid;
	std::vector<double> scan_window_upper;
	std::vector<bool> scan_window_upper_valid;
	std::vector<int32_t> ms1_scan_index;
	std::vector<bool> ms1_scan_index_valid;

	[[nodiscard]] size_t size() const {
		return spectrum_index.size();
	}
	[[nodiscard]] bool empty() const {
		return spectrum_index.empty();
	}
};

// Backward compatibility
using MzMLSpectrumBatch = SpectrumBatch;

} // namespace miint
