#include "MzXMLReader.hpp"
#include "MsBinaryDecoder.hpp"
#include "mz_parse_utils.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <expat.h>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace miint {

struct MzXMLScanState {
	int32_t spectrum_index = 0;
	int32_t scan_num = 0;
	int32_t ms_level = 0;
	std::string spectrum_type;
	std::string polarity;
	double base_peak_mz = 0;
	bool base_peak_mz_valid = false;
	double base_peak_intensity = 0;
	bool base_peak_intensity_valid = false;
	double total_ion_current = 0;
	bool total_ion_current_valid = false;
	double lowest_mz = 0;
	bool lowest_mz_valid = false;
	double highest_mz = 0;
	bool highest_mz_valid = false;
	double retention_time = 0;
	bool retention_time_valid = false;
	std::string filter_string;
	double scan_window_lower = 0;
	bool scan_window_lower_valid = false;
	double scan_window_upper = 0;
	bool scan_window_upper_valid = false;

	double precursor_mz = 0;
	bool precursor_mz_valid = false;
	int32_t precursor_charge = 0;
	bool precursor_charge_valid = false;
	double precursor_intensity = 0;
	bool precursor_intensity_valid = false;
	std::string activation_method;
	double collision_energy = 0;
	bool collision_energy_valid = false;
	double window_wideness = 0;
	bool window_wideness_valid = false;
	int32_t precursor_scan_num = 0;
	bool precursor_scan_num_valid = false;

	bool peaks_compressed = false;
	bool peaks_64bit = true;
	std::string peaks_text;
	std::string precursor_mz_text;

	std::vector<double> mz_array;
	std::vector<double> intensity_array;

	int32_t ms1_scan_index = 0;
	bool ms1_scan_index_valid = false;

	// Index into completed_spectra when this scan opened (for parent backfill)
	size_t children_start = 0;
};

// Parse ISO 8601 duration: PT{h}H{m}M{s}S (any subset) -> minutes
static bool parse_retention_time(const char *rt_str, double &minutes) {
	if (!rt_str || rt_str[0] != 'P' || rt_str[1] != 'T') {
		return false;
	}
	const char *p = rt_str + 2;
	double total_minutes = 0;
	bool found_any = false;

	while (*p) {
		char *end = nullptr;
		double val = std::strtod(p, &end);
		if (end == p) {
			return false;
		}
		if (*end == 'H') {
			total_minutes += val * 60.0;
			found_any = true;
		} else if (*end == 'M') {
			total_minutes += val;
			found_any = true;
		} else if (*end == 'S') {
			total_minutes += val / 60.0;
			found_any = true;
		} else {
			return false;
		}
		p = end + 1;
	}
	if (found_any) {
		minutes = total_minutes;
	}
	return found_any;
}

enum class MzXMLContext {
	NONE,
	MZXML,
	MSRUN,
	SCAN,
	PRECURSOR_MZ,
	PEAKS,
};

struct MzXMLReader::Impl {
	std::string filepath;
	FILE *file = nullptr;
	std::function<size_t(void *, size_t)> read_fn;
	XML_Parser parser = nullptr;

	std::vector<MzXMLContext> context_stack;
	std::vector<MzXMLScanState> scan_stack;
	std::deque<MzXMLScanState> completed_spectra;

	int32_t next_spectrum_index = 0;
	int32_t last_ms1_index = -1;
	std::unordered_map<int32_t, int32_t> scan_num_to_index;

	bool parse_done = false;
	std::exception_ptr pending_exception;

	Impl(const std::string &path) : filepath(path) {
		file = std::fopen(path.c_str(), "rb");
		if (!file) {
			throw std::runtime_error("MzXMLReader: cannot open file: " + path);
		}
		read_fn = [this](void *buf, size_t n) {
			return std::fread(buf, 1, n, file);
		};
		init_parser();
	}

	// Constructor for DuckDB FileSystem access (supports local + remote).
	// Guarded because the test binary does not link DuckDB.
#ifdef MIINT_STATIC_BUILD
	Impl(duckdb::FileSystem &fs, const std::string &path) : filepath(path) {
		auto handle = fs.OpenFile(path, duckdb::FileOpenFlags(duckdb::FileOpenFlags::FILE_FLAGS_READ));
		if (!handle) {
			throw std::runtime_error("MzXMLReader: cannot open file: " + path);
		}
		auto shared_handle = std::shared_ptr<duckdb::FileHandle>(handle.release());
		read_fn = [shared_handle](void *buf, size_t n) -> size_t {
			auto result = shared_handle->Read(buf, n);
			if (result < 0) {
				throw std::runtime_error("MzXMLReader: I/O error reading file");
			}
			return static_cast<size_t>(result);
		};
		init_parser();
	}
#endif

	~Impl() {
		if (parser) {
			XML_ParserFree(parser);
		}
		if (file) {
			std::fclose(file);
		}
	}

	void init_parser() {
		parser = XML_ParserCreateNS(nullptr, '|');
		if (!parser) {
			throw std::runtime_error("MzXMLReader: failed to create XML parser");
		}
		XML_SetUserData(parser, this);
		XML_SetElementHandler(parser, start_element_handler, end_element_handler);
		XML_SetCharacterDataHandler(parser, char_data_handler);
		context_stack.push_back(MzXMLContext::NONE);
	}

	MzXMLContext current_context() const {
		return context_stack.back();
	}

	void check_pending_exception() {
		if (pending_exception) {
			std::rethrow_exception(pending_exception);
		}
	}

	void parse_until(size_t spectra_needed) {
		if (parse_done) {
			return;
		}

		char buffer[65536];
		while (true) {
			if (completed_spectra.size() >= spectra_needed) {
				return;
			}

			size_t bytes_read = read_fn(buffer, sizeof(buffer));
			bool is_final = (bytes_read == 0);

			auto status = XML_Parse(parser, buffer, static_cast<int>(bytes_read), is_final);
			check_pending_exception();

			if (status == XML_STATUS_ERROR) {
				auto line = XML_GetCurrentLineNumber(parser);
				auto err = XML_ErrorString(XML_GetErrorCode(parser));
				throw std::runtime_error("MzXMLReader: XML parse error at line " + std::to_string(line) + " in " +
				                         filepath + ": " + std::string(err));
			}

			if (is_final) {
				parse_done = true;
				return;
			}

			if (completed_spectra.size() >= spectra_needed) {
				return;
			}
		}
	}

	static const char *get_attr(const char **attrs, const char *name) {
		for (int i = 0; attrs[i]; i += 2) {
			if (std::strcmp(attrs[i], name) == 0) {
				return attrs[i + 1];
			}
		}
		return nullptr;
	}

	static const char *local_name(const char *full_name) {
		const char *pipe = std::strrchr(full_name, '|');
		return pipe ? pipe + 1 : full_name;
	}

	static void XMLCALL start_element_handler(void *user_data, const char *name, const char **attrs) {
		auto *self = static_cast<Impl *>(user_data);
		if (self->pending_exception) {
			return;
		}
		try {
			self->on_start_element(local_name(name), attrs);
		} catch (...) {
			self->pending_exception = std::current_exception();
			XML_StopParser(self->parser, XML_FALSE);
		}
	}

	static void XMLCALL end_element_handler(void *user_data, const char *name) {
		auto *self = static_cast<Impl *>(user_data);
		if (self->pending_exception) {
			return;
		}
		try {
			self->on_end_element(local_name(name));
		} catch (...) {
			self->pending_exception = std::current_exception();
			XML_StopParser(self->parser, XML_FALSE);
		}
	}

	static void XMLCALL char_data_handler(void *user_data, const char *s, int len) {
		auto *self = static_cast<Impl *>(user_data);
		if (self->pending_exception) {
			return;
		}
		try {
			self->on_char_data(s, len);
		} catch (...) {
			self->pending_exception = std::current_exception();
			XML_StopParser(self->parser, XML_FALSE);
		}
	}

	static constexpr size_t MAX_CONTEXT_DEPTH = 64;

	void on_start_element(const char *name, const char **attrs) {
		if (context_stack.size() >= MAX_CONTEXT_DEPTH) {
			throw std::runtime_error("MzXMLReader: XML nesting too deep (>" + std::to_string(MAX_CONTEXT_DEPTH) +
			                         " levels) in " + filepath);
		}
		auto ctx = current_context();

		if (std::strcmp(name, "mzXML") == 0 && ctx == MzXMLContext::NONE) {
			context_stack.push_back(MzXMLContext::MZXML);
			return;
		}

		if (std::strcmp(name, "msRun") == 0 && ctx == MzXMLContext::MZXML) {
			context_stack.push_back(MzXMLContext::MSRUN);
			return;
		}

		if (std::strcmp(name, "scan") == 0 && (ctx == MzXMLContext::MSRUN || ctx == MzXMLContext::SCAN)) {
			context_stack.push_back(MzXMLContext::SCAN);

			MzXMLScanState scan;
			scan.children_start = completed_spectra.size();

			safe_stoi(get_attr(attrs, "num"), scan.scan_num);
			safe_stoi(get_attr(attrs, "msLevel"), scan.ms_level);

			auto *rt = get_attr(attrs, "retentionTime");
			if (rt) {
				scan.retention_time_valid = parse_retention_time(rt, scan.retention_time);
			}

			auto *centroided = get_attr(attrs, "centroided");
			if (centroided) {
				if (std::strcmp(centroided, "1") == 0) {
					scan.spectrum_type = "centroid";
				} else if (std::strcmp(centroided, "0") == 0) {
					scan.spectrum_type = "profile";
				}
			}

			auto *polarity = get_attr(attrs, "polarity");
			if (polarity) {
				if (std::strcmp(polarity, "+") == 0) {
					scan.polarity = "positive";
				} else if (std::strcmp(polarity, "-") == 0) {
					scan.polarity = "negative";
				}
			}

			scan.base_peak_mz_valid = safe_stod(get_attr(attrs, "basePeakMz"), scan.base_peak_mz);
			scan.base_peak_intensity_valid = safe_stod(get_attr(attrs, "basePeakIntensity"), scan.base_peak_intensity);
			scan.total_ion_current_valid = safe_stod(get_attr(attrs, "totIonCurrent"), scan.total_ion_current);
			scan.lowest_mz_valid = safe_stod(get_attr(attrs, "lowMz"), scan.lowest_mz);
			scan.highest_mz_valid = safe_stod(get_attr(attrs, "highMz"), scan.highest_mz);

			auto *filter = get_attr(attrs, "filterLine");
			if (filter) {
				scan.filter_string = filter;
			}

			scan.scan_window_lower_valid = safe_stod(get_attr(attrs, "startMz"), scan.scan_window_lower);
			scan.scan_window_upper_valid = safe_stod(get_attr(attrs, "endMz"), scan.scan_window_upper);
			scan.collision_energy_valid = safe_stod(get_attr(attrs, "collisionEnergy"), scan.collision_energy);

			scan_stack.push_back(std::move(scan));
			return;
		}

		// mzXML spec allows at most one <precursorMz> per scan; if multiple appear, last wins.
		if (std::strcmp(name, "precursorMz") == 0 && ctx == MzXMLContext::SCAN && !scan_stack.empty()) {
			context_stack.push_back(MzXMLContext::PRECURSOR_MZ);
			auto &scan = scan_stack.back();
			scan.precursor_mz_text.clear();

			scan.precursor_charge_valid = safe_stoi(get_attr(attrs, "precursorCharge"), scan.precursor_charge);
			scan.precursor_intensity_valid = safe_stod(get_attr(attrs, "precursorIntensity"), scan.precursor_intensity);
			scan.window_wideness_valid = safe_stod(get_attr(attrs, "windowWideness"), scan.window_wideness);
			scan.precursor_scan_num_valid = safe_stoi(get_attr(attrs, "precursorScanNum"), scan.precursor_scan_num);

			auto *method = get_attr(attrs, "activationMethod");
			if (method) {
				scan.activation_method = method;
			}
			return;
		}

		if (std::strcmp(name, "peaks") == 0 && ctx == MzXMLContext::SCAN && !scan_stack.empty()) {
			context_stack.push_back(MzXMLContext::PEAKS);
			auto &scan = scan_stack.back();
			scan.peaks_text.clear();

			auto *precision = get_attr(attrs, "precision");
			if (precision) {
				scan.peaks_64bit = (std::strcmp(precision, "64") == 0);
			}

			// compressionType: absent or "none" → uncompressed, "zlib" → zlib
			auto *compression = get_attr(attrs, "compressionType");
			if (compression) {
				if (std::strcmp(compression, "zlib") == 0) {
					scan.peaks_compressed = true;
				} else if (std::strcmp(compression, "none") == 0) {
					scan.peaks_compressed = false;
				} else {
					throw std::runtime_error("MzXMLReader: unsupported compressionType '" + std::string(compression) +
					                         "' in " + filepath);
				}
			}

			auto *byte_order = get_attr(attrs, "byteOrder");
			if (byte_order && std::strcmp(byte_order, "network") != 0) {
				throw std::runtime_error("MzXMLReader: unsupported byteOrder '" + std::string(byte_order) +
				                         "' (expected 'network') in " + filepath);
			}

			auto *pair_order = get_attr(attrs, "pairOrder");
			if (pair_order && std::strcmp(pair_order, "m/z-int") != 0) {
				throw std::runtime_error("MzXMLReader: unsupported pairOrder '" + std::string(pair_order) +
				                         "' (expected 'm/z-int') in " + filepath);
			}

			auto *content_type = get_attr(attrs, "contentType");
			if (content_type && std::strcmp(content_type, "m/z-int") != 0) {
				throw std::runtime_error("MzXMLReader: unsupported contentType '" + std::string(content_type) +
				                         "' (expected 'm/z-int') in " + filepath);
			}
			return;
		}
	}

	void complete_scan(MzXMLScanState &scan) {
		scan.spectrum_index = next_spectrum_index++;
		scan_num_to_index[scan.scan_num] = scan.spectrum_index;

		if (scan.ms_level == 1) {
			// Backfill: nested children emitted while this scan was open
			for (size_t i = scan.children_start; i < completed_spectra.size(); i++) {
				auto &child = completed_spectra[i];
				if (child.ms_level > 1 && !child.ms1_scan_index_valid) {
					child.ms1_scan_index = scan.spectrum_index;
					child.ms1_scan_index_valid = true;
				}
			}
			last_ms1_index = scan.spectrum_index;
		} else if (scan.ms_level > 1) {
			// If any ancestor on the scan stack is MS1, it will backfill us when it completes.
			// This handles MS3 nested inside MS2 nested inside MS1.
			bool ancestor_will_backfill = false;
			for (const auto &ancestor : scan_stack) {
				if (ancestor.ms_level == 1) {
					ancestor_will_backfill = true;
					break;
				}
			}
			if (!ancestor_will_backfill) {
				// Priority 2: precursorScanNum
				if (scan.precursor_scan_num_valid) {
					auto it = scan_num_to_index.find(scan.precursor_scan_num);
					if (it != scan_num_to_index.end()) {
						scan.ms1_scan_index = it->second;
						scan.ms1_scan_index_valid = true;
					}
				}
				// Priority 3: last_ms1_index
				if (!scan.ms1_scan_index_valid && last_ms1_index >= 0) {
					scan.ms1_scan_index = last_ms1_index;
					scan.ms1_scan_index_valid = true;
				}
			}
		}

		completed_spectra.push_back(std::move(scan));
	}

	void on_end_element(const char *name) {
		auto ctx = current_context();

		if (std::strcmp(name, "mzXML") == 0 && ctx == MzXMLContext::MZXML) {
			context_stack.pop_back();
			return;
		}

		if (std::strcmp(name, "msRun") == 0 && ctx == MzXMLContext::MSRUN) {
			context_stack.pop_back();
			return;
		}

		if (std::strcmp(name, "precursorMz") == 0 && ctx == MzXMLContext::PRECURSOR_MZ) {
			if (!scan_stack.empty()) {
				auto &scan = scan_stack.back();
				scan.precursor_mz_valid = safe_stod(scan.precursor_mz_text.c_str(), scan.precursor_mz);
			}
			context_stack.pop_back();
			return;
		}

		if (std::strcmp(name, "peaks") == 0 && ctx == MzXMLContext::PEAKS) {
			if (!scan_stack.empty()) {
				auto &scan = scan_stack.back();
				auto &text = scan.peaks_text;
				while (!text.empty() &&
				       (text.back() == ' ' || text.back() == '\n' || text.back() == '\r' || text.back() == '\t')) {
					text.pop_back();
				}
				if (!text.empty()) {
					auto [mz, intensity] = MsBinaryDecoder::decode_mzxml(text, scan.peaks_compressed, scan.peaks_64bit);
					scan.mz_array = std::move(mz);
					scan.intensity_array = std::move(intensity);
				}
				scan.peaks_text.clear();
				scan.peaks_text.shrink_to_fit();
			}
			context_stack.pop_back();
			return;
		}

		if (std::strcmp(name, "scan") == 0 && ctx == MzXMLContext::SCAN) {
			if (!scan_stack.empty()) {
				auto scan = std::move(scan_stack.back());
				scan_stack.pop_back();
				complete_scan(scan);
			}
			context_stack.pop_back();
			return;
		}
	}

	void on_char_data(const char *s, int len) {
		auto ctx = current_context();
		if (!scan_stack.empty()) {
			if (ctx == MzXMLContext::PEAKS) {
				scan_stack.back().peaks_text.append(s, len);
			} else if (ctx == MzXMLContext::PRECURSOR_MZ) {
				scan_stack.back().precursor_mz_text.append(s, len);
			}
		}
	}
};

// Public API

MzXMLReader::MzXMLReader(const std::string &path) : impl_(std::make_unique<Impl>(path)) {
}

#ifdef MIINT_STATIC_BUILD
MzXMLReader::MzXMLReader(duckdb::FileSystem &fs, const std::string &path) : impl_(std::make_unique<Impl>(fs, path)) {
}
#endif

MzXMLReader::~MzXMLReader() = default;
MzXMLReader::MzXMLReader(MzXMLReader &&) noexcept = default;
MzXMLReader &MzXMLReader::operator=(MzXMLReader &&) noexcept = default;

static void drain_scan_to_batch(MzMLSpectrumBatch &batch, MzXMLScanState &&scan) {
	batch.spectrum_index.push_back(scan.spectrum_index);
	batch.spectrum_id.push_back("scan=" + std::to_string(scan.scan_num));
	batch.scan_number.push_back(scan.scan_num);
	batch.scan_number_valid.push_back(true);
	batch.ms_level.push_back(scan.ms_level);
	batch.retention_time.push_back(scan.retention_time);
	batch.retention_time_valid.push_back(scan.retention_time_valid);
	batch.spectrum_type.push_back(std::move(scan.spectrum_type));
	batch.polarity.push_back(std::move(scan.polarity));
	batch.base_peak_mz.push_back(scan.base_peak_mz);
	batch.base_peak_mz_valid.push_back(scan.base_peak_mz_valid);
	batch.base_peak_intensity.push_back(scan.base_peak_intensity);
	batch.base_peak_intensity_valid.push_back(scan.base_peak_intensity_valid);
	batch.total_ion_current.push_back(scan.total_ion_current);
	batch.total_ion_current_valid.push_back(scan.total_ion_current_valid);
	batch.lowest_mz.push_back(scan.lowest_mz);
	batch.lowest_mz_valid.push_back(scan.lowest_mz_valid);
	batch.highest_mz.push_back(scan.highest_mz);
	batch.highest_mz_valid.push_back(scan.highest_mz_valid);
	batch.default_array_length.push_back(static_cast<int32_t>(scan.mz_array.size()));
	batch.precursor_mz.push_back(scan.precursor_mz);
	batch.precursor_mz_valid.push_back(scan.precursor_mz_valid);
	batch.precursor_charge.push_back(scan.precursor_charge);
	batch.precursor_charge_valid.push_back(scan.precursor_charge_valid);
	batch.precursor_intensity.push_back(scan.precursor_intensity);
	batch.precursor_intensity_valid.push_back(scan.precursor_intensity_valid);

	if (scan.precursor_mz_valid) {
		batch.isolation_window_target.push_back(scan.precursor_mz);
		batch.isolation_window_target_valid.push_back(true);
	} else {
		batch.isolation_window_target.push_back(0);
		batch.isolation_window_target_valid.push_back(false);
	}
	if (scan.window_wideness_valid) {
		double half = scan.window_wideness / 2.0;
		batch.isolation_window_lower.push_back(half);
		batch.isolation_window_lower_valid.push_back(true);
		batch.isolation_window_upper.push_back(half);
		batch.isolation_window_upper_valid.push_back(true);
	} else {
		batch.isolation_window_lower.push_back(0);
		batch.isolation_window_lower_valid.push_back(false);
		batch.isolation_window_upper.push_back(0);
		batch.isolation_window_upper_valid.push_back(false);
	}

	batch.activation_method.push_back(std::move(scan.activation_method));
	batch.collision_energy.push_back(scan.collision_energy);
	batch.collision_energy_valid.push_back(scan.collision_energy_valid);
	batch.mz_array.push_back(std::move(scan.mz_array));
	batch.intensity_array.push_back(std::move(scan.intensity_array));
	batch.filter_string.push_back(std::move(scan.filter_string));
	batch.scan_window_lower.push_back(scan.scan_window_lower);
	batch.scan_window_lower_valid.push_back(scan.scan_window_lower_valid);
	batch.scan_window_upper.push_back(scan.scan_window_upper);
	batch.scan_window_upper_valid.push_back(scan.scan_window_upper_valid);
	batch.ms1_scan_index.push_back(scan.ms1_scan_index);
	batch.ms1_scan_index_valid.push_back(scan.ms1_scan_index_valid);
}

MzMLSpectrumBatch MzXMLReader::read_spectra(size_t n) {
	MzMLSpectrumBatch batch;
	if (n == 0) {
		return batch;
	}

	if (impl_->completed_spectra.size() < n && !impl_->parse_done) {
		impl_->parse_until(n);
	}

	size_t count = std::min(n, impl_->completed_spectra.size());

	batch.spectrum_index.reserve(count);
	batch.spectrum_id.reserve(count);
	batch.scan_number.reserve(count);
	batch.scan_number_valid.reserve(count);
	batch.ms_level.reserve(count);
	batch.retention_time.reserve(count);
	batch.retention_time_valid.reserve(count);
	batch.spectrum_type.reserve(count);
	batch.polarity.reserve(count);
	batch.base_peak_mz.reserve(count);
	batch.base_peak_mz_valid.reserve(count);
	batch.base_peak_intensity.reserve(count);
	batch.base_peak_intensity_valid.reserve(count);
	batch.total_ion_current.reserve(count);
	batch.total_ion_current_valid.reserve(count);
	batch.lowest_mz.reserve(count);
	batch.lowest_mz_valid.reserve(count);
	batch.highest_mz.reserve(count);
	batch.highest_mz_valid.reserve(count);
	batch.default_array_length.reserve(count);
	batch.precursor_mz.reserve(count);
	batch.precursor_mz_valid.reserve(count);
	batch.precursor_charge.reserve(count);
	batch.precursor_charge_valid.reserve(count);
	batch.precursor_intensity.reserve(count);
	batch.precursor_intensity_valid.reserve(count);
	batch.isolation_window_target.reserve(count);
	batch.isolation_window_target_valid.reserve(count);
	batch.isolation_window_lower.reserve(count);
	batch.isolation_window_lower_valid.reserve(count);
	batch.isolation_window_upper.reserve(count);
	batch.isolation_window_upper_valid.reserve(count);
	batch.activation_method.reserve(count);
	batch.collision_energy.reserve(count);
	batch.collision_energy_valid.reserve(count);
	batch.mz_array.reserve(count);
	batch.intensity_array.reserve(count);
	batch.filter_string.reserve(count);
	batch.scan_window_lower.reserve(count);
	batch.scan_window_lower_valid.reserve(count);
	batch.scan_window_upper.reserve(count);
	batch.scan_window_upper_valid.reserve(count);
	batch.ms1_scan_index.reserve(count);
	batch.ms1_scan_index_valid.reserve(count);

	for (size_t i = 0; i < count; i++) {
		drain_scan_to_batch(batch, std::move(impl_->completed_spectra.front()));
		impl_->completed_spectra.pop_front();
	}

	return batch;
}

} // namespace miint
