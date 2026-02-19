#include "MzMLReader.hpp"
#include "MzMLBinaryDecoder.hpp"
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <expat.h>
#include <stdexcept>
#include <unordered_map>

namespace miint {

// cvParam attributes parsed from XML
struct CvParam {
	std::string accession;
	std::string value;
	std::string unit_accession;
};

// Per-binaryDataArray state accumulated during parsing
struct BinaryDataArrayState {
	bool is_compressed = false;
	bool is_64bit = true;
	bool is_mz = false;
	bool is_intensity = false;
	bool is_time = false;
	std::string binary_text;
};

// Per-spectrum state accumulated during parsing
struct SpectrumState {
	int32_t index = 0;
	std::string id;
	int32_t default_array_length = 0;
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
	double precursor_mz = 0;
	bool precursor_mz_valid = false;
	int32_t precursor_charge = 0;
	bool precursor_charge_valid = false;
	double precursor_intensity = 0;
	bool precursor_intensity_valid = false;
	double isolation_window_target = 0;
	bool isolation_window_target_valid = false;
	double isolation_window_lower = 0;
	bool isolation_window_lower_valid = false;
	double isolation_window_upper = 0;
	bool isolation_window_upper_valid = false;
	std::string activation_method;
	double collision_energy = 0;
	bool collision_energy_valid = false;
	std::string filter_string;
	double scan_window_lower = 0;
	bool scan_window_lower_valid = false;
	double scan_window_upper = 0;
	bool scan_window_upper_valid = false;
	double retention_time = 0;
	bool retention_time_valid = false;
	std::vector<double> mz_array;
	std::vector<double> intensity_array;
	// ms1_scan_index: recorded at parse time (not batch-drain time)
	int32_t ms1_scan_index = 0;
	bool ms1_scan_index_valid = false;
};

// Per-chromatogram state accumulated during parsing
struct ChromatogramState {
	int32_t index = 0;
	std::string id;
	std::string chromatogram_type;
	double precursor_mz = 0;
	bool precursor_mz_valid = false;
	double product_mz = 0;
	bool product_mz_valid = false;
	std::vector<double> time_array;
	std::vector<double> intensity_array;
};

// XML element path tracking
enum class ParseContext {
	NONE,
	MZML,
	REFERENCEABLE_PARAM_GROUP_LIST,
	REFERENCEABLE_PARAM_GROUP,
	RUN,
	SPECTRUM_LIST,
	SPECTRUM,
	SPECTRUM_SCAN_LIST,
	SPECTRUM_SCAN,
	SPECTRUM_SCAN_WINDOW,
	SPECTRUM_PRECURSOR_LIST,
	SPECTRUM_PRECURSOR,
	SPECTRUM_ISOLATION_WINDOW,
	SPECTRUM_SELECTED_ION_LIST,
	SPECTRUM_SELECTED_ION,
	SPECTRUM_ACTIVATION,
	SPECTRUM_BINARY_DATA_ARRAY_LIST,
	SPECTRUM_BINARY_DATA_ARRAY,
	SPECTRUM_BINARY,
	CHROMATOGRAM_LIST,
	CHROMATOGRAM,
	CHROMATOGRAM_PRECURSOR,
	CHROMATOGRAM_PRECURSOR_ISOLATION_WINDOW,
	CHROMATOGRAM_PRODUCT,
	CHROMATOGRAM_PRODUCT_ISOLATION_WINDOW,
	CHROMATOGRAM_BINARY_DATA_ARRAY_LIST,
	CHROMATOGRAM_BINARY_DATA_ARRAY,
	CHROMATOGRAM_BINARY,
};

// Safe numeric conversion helpers: return false on failure instead of throwing
static bool safe_stoi(const std::string &s, int32_t &out) {
	if (s.empty()) {
		return false;
	}
	try {
		size_t pos = 0;
		long val = std::stol(s, &pos);
		if (pos != s.size() || val < INT32_MIN || val > INT32_MAX) {
			return false;
		}
		out = static_cast<int32_t>(val);
		return true;
	} catch (...) {
		return false;
	}
}

static bool safe_stod(const std::string &s, double &out) {
	if (s.empty()) {
		return false;
	}
	try {
		size_t pos = 0;
		out = std::stod(s, &pos);
		return pos == s.size();
	} catch (...) {
		return false;
	}
}

struct MzMLReader::Impl {
	std::string filepath;
	FILE *file = nullptr;
	XML_Parser parser = nullptr;

	// Referenceable param groups: id -> list of CvParams
	std::unordered_map<std::string, std::vector<CvParam>> param_groups;

	// Parse state
	std::vector<ParseContext> context_stack;
	std::string current_param_group_id;

	// Current spectrum/chromatogram being built
	SpectrumState current_spectrum;
	BinaryDataArrayState current_binary_array;
	ChromatogramState current_chromatogram;

	// ms1_scan_index tracking: updated at parse time
	int32_t last_ms1_index = -1;

	// Completed items waiting to be returned (deque for O(1) pop_front)
	std::deque<SpectrumState> completed_spectra;
	std::deque<ChromatogramState> completed_chromatograms;

	// Parsing progress flags
	bool spectra_done = false;
	bool chromatograms_done = false;
	bool parse_done = false;
	bool in_chromatogram_mode = false;

	// Did we see a spectrumList / chromatogramList?
	bool has_spectrum_list = false;
	bool has_chromatogram_list = false;

	// Exception captured from SAX callback (expat is C — cannot throw through it)
	std::exception_ptr pending_exception;

	Impl(const std::string &path) : filepath(path) {
		file = std::fopen(path.c_str(), "rb");
		if (!file) {
			throw std::runtime_error("MzMLReader: cannot open file: " + path);
		}

		parser = XML_ParserCreateNS(nullptr, '|');
		if (!parser) {
			std::fclose(file);
			throw std::runtime_error("MzMLReader: failed to create XML parser");
		}

		XML_SetUserData(parser, this);
		XML_SetElementHandler(parser, start_element_handler, end_element_handler);
		XML_SetCharacterDataHandler(parser, char_data_handler);

		context_stack.push_back(ParseContext::NONE);
	}

	~Impl() {
		if (parser) {
			XML_ParserFree(parser);
		}
		if (file) {
			std::fclose(file);
		}
	}

	ParseContext current_context() const {
		return context_stack.back();
	}

	// Rethrow any exception captured from a SAX callback
	void check_pending_exception() {
		if (pending_exception) {
			std::rethrow_exception(pending_exception);
		}
	}

	// Parse until we have enough completed items or EOF
	void parse_until(size_t spectra_needed, size_t chromatograms_needed) {
		if (parse_done) {
			return;
		}

		char buffer[65536];
		while (true) {
			// Check if we have enough data
			if (!in_chromatogram_mode && completed_spectra.size() >= spectra_needed) {
				return;
			}
			if (in_chromatogram_mode && completed_chromatograms.size() >= chromatograms_needed) {
				return;
			}
			// When fast-forwarding spectra, stop as soon as spectra_done is set
			if (!in_chromatogram_mode && spectra_done) {
				return;
			}

			size_t bytes_read = std::fread(buffer, 1, sizeof(buffer), file);
			bool is_final = (bytes_read == 0);

			auto status = XML_Parse(parser, buffer, static_cast<int>(bytes_read), is_final);

			// Check for exceptions captured from SAX callbacks first
			check_pending_exception();

			if (status == XML_STATUS_ERROR) {
				auto line = XML_GetCurrentLineNumber(parser);
				auto err = XML_ErrorString(XML_GetErrorCode(parser));
				throw std::runtime_error("MzMLReader: XML parse error at line " + std::to_string(line) + " in " +
				                         filepath + ": " + std::string(err));
			}

			if (is_final) {
				parse_done = true;
				spectra_done = true;
				chromatograms_done = true;
				return;
			}

			// Check again after parsing
			if (!in_chromatogram_mode && completed_spectra.size() >= spectra_needed) {
				return;
			}
			if (in_chromatogram_mode && completed_chromatograms.size() >= chromatograms_needed) {
				return;
			}
			if (!in_chromatogram_mode && spectra_done) {
				return;
			}
		}
	}

	// Apply resolved param group cvParams to the current binary data array
	void apply_cv_param_to_binary_array(const CvParam &cv) {
		if (cv.accession == "MS:1000574") {
			current_binary_array.is_compressed = true;
		} else if (cv.accession == "MS:1000576") {
			current_binary_array.is_compressed = false;
		} else if (cv.accession == "MS:1000523") {
			current_binary_array.is_64bit = true;
		} else if (cv.accession == "MS:1000521" || cv.accession == "MS:1000519") {
			current_binary_array.is_64bit = false;
		} else if (cv.accession == "MS:1000514") {
			current_binary_array.is_mz = true;
		} else if (cv.accession == "MS:1000515") {
			current_binary_array.is_intensity = true;
		} else if (cv.accession == "MS:1000595") {
			current_binary_array.is_time = true;
		}
	}

	void apply_cv_param_to_spectrum(const CvParam &cv) {
		const auto &acc = cv.accession;
		const auto &val = cv.value;

		if (acc == "MS:1000511") {
			safe_stoi(val, current_spectrum.ms_level);
		} else if (acc == "MS:1000127") {
			current_spectrum.spectrum_type = "centroid";
		} else if (acc == "MS:1000128") {
			current_spectrum.spectrum_type = "profile";
		} else if (acc == "MS:1000130") {
			current_spectrum.polarity = "positive";
		} else if (acc == "MS:1000129") {
			current_spectrum.polarity = "negative";
		} else if (acc == "MS:1000504") {
			current_spectrum.base_peak_mz_valid = safe_stod(val, current_spectrum.base_peak_mz);
		} else if (acc == "MS:1000505") {
			current_spectrum.base_peak_intensity_valid = safe_stod(val, current_spectrum.base_peak_intensity);
		} else if (acc == "MS:1000285") {
			current_spectrum.total_ion_current_valid = safe_stod(val, current_spectrum.total_ion_current);
		} else if (acc == "MS:1000528") {
			current_spectrum.lowest_mz_valid = safe_stod(val, current_spectrum.lowest_mz);
		} else if (acc == "MS:1000527") {
			current_spectrum.highest_mz_valid = safe_stod(val, current_spectrum.highest_mz);
		} else if (acc == "MS:1000512") {
			current_spectrum.filter_string = val;
		}
	}

	void apply_cv_param_to_scan(const CvParam &cv) {
		if (cv.accession == "MS:1000016") {
			current_spectrum.retention_time_valid = safe_stod(cv.value, current_spectrum.retention_time);
			// Convert seconds to minutes if unit is seconds
			if (current_spectrum.retention_time_valid && cv.unit_accession == "UO:0000010") {
				current_spectrum.retention_time /= 60.0;
			}
		}
	}

	void apply_cv_param_to_scan_window(const CvParam &cv) {
		// Note: if multiple scan windows exist, later values overwrite earlier ones.
		// The mzML spec allows multiple scan windows, but we only report the last one.
		if (cv.accession == "MS:1000501") {
			current_spectrum.scan_window_lower_valid = safe_stod(cv.value, current_spectrum.scan_window_lower);
		} else if (cv.accession == "MS:1000500") {
			current_spectrum.scan_window_upper_valid = safe_stod(cv.value, current_spectrum.scan_window_upper);
		}
	}

	void apply_cv_param_to_selected_ion(const CvParam &cv) {
		if (cv.accession == "MS:1000744") {
			current_spectrum.precursor_mz_valid = safe_stod(cv.value, current_spectrum.precursor_mz);
		} else if (cv.accession == "MS:1000041") {
			current_spectrum.precursor_charge_valid = safe_stoi(cv.value, current_spectrum.precursor_charge);
		} else if (cv.accession == "MS:1000042") {
			current_spectrum.precursor_intensity_valid = safe_stod(cv.value, current_spectrum.precursor_intensity);
		}
	}

	void apply_cv_param_to_isolation_window(const CvParam &cv) {
		if (cv.accession == "MS:1000827") {
			current_spectrum.isolation_window_target_valid =
			    safe_stod(cv.value, current_spectrum.isolation_window_target);
		} else if (cv.accession == "MS:1000828") {
			current_spectrum.isolation_window_lower_valid =
			    safe_stod(cv.value, current_spectrum.isolation_window_lower);
		} else if (cv.accession == "MS:1000829") {
			current_spectrum.isolation_window_upper_valid =
			    safe_stod(cv.value, current_spectrum.isolation_window_upper);
		}
	}

	void apply_cv_param_to_activation(const CvParam &cv) {
		if (cv.accession == "MS:1000133") {
			current_spectrum.activation_method = "CID";
		} else if (cv.accession == "MS:1000422") {
			current_spectrum.activation_method = "HCD";
		} else if (cv.accession == "MS:1000598") {
			current_spectrum.activation_method = "ETD";
		} else if (cv.accession == "MS:1000045") {
			current_spectrum.collision_energy_valid = safe_stod(cv.value, current_spectrum.collision_energy);
		}
	}

	void apply_cv_param_to_chromatogram(const CvParam &cv) {
		if (cv.accession == "MS:1000235") {
			current_chromatogram.chromatogram_type = "TIC";
		} else if (cv.accession == "MS:1000628") {
			current_chromatogram.chromatogram_type = "BPC";
		} else if (cv.accession == "MS:1001473" || cv.accession == "MS:1000789") {
			current_chromatogram.chromatogram_type = "SRM";
		} else if (cv.accession == "MS:1000627") {
			current_chromatogram.chromatogram_type = "SIC";
		}
	}

	void apply_cv_param_to_chrom_precursor_isolation(const CvParam &cv) {
		if (cv.accession == "MS:1000827") {
			current_chromatogram.precursor_mz_valid = safe_stod(cv.value, current_chromatogram.precursor_mz);
		}
	}

	void apply_cv_param_to_chrom_product_isolation(const CvParam &cv) {
		if (cv.accession == "MS:1000827") {
			current_chromatogram.product_mz_valid = safe_stod(cv.value, current_chromatogram.product_mz);
		}
	}

	// Helper to extract an attribute value from expat's attrs array
	static const char *get_attr(const char **attrs, const char *name) {
		for (int i = 0; attrs[i]; i += 2) {
			if (std::strcmp(attrs[i], name) == 0) {
				return attrs[i + 1];
			}
		}
		return nullptr;
	}

	// Strip namespace prefix from element name (expat with NS uses "ns|local")
	static const char *local_name(const char *full_name) {
		const char *pipe = std::strrchr(full_name, '|');
		return pipe ? pipe + 1 : full_name;
	}

	CvParam parse_cv_param(const char **attrs) {
		CvParam cv;
		for (int i = 0; attrs[i]; i += 2) {
			if (std::strcmp(attrs[i], "accession") == 0) {
				cv.accession = attrs[i + 1];
			} else if (std::strcmp(attrs[i], "value") == 0) {
				cv.value = attrs[i + 1];
			} else if (std::strcmp(attrs[i], "unitAccession") == 0) {
				cv.unit_accession = attrs[i + 1];
			}
		}
		return cv;
	}

	// Static SAX callbacks — catch all C++ exceptions to avoid UB propagating through expat's C frames
	static void XMLCALL start_element_handler(void *user_data, const char *name, const char **attrs) {
		auto *self = static_cast<Impl *>(user_data);
		if (self->pending_exception) {
			return; // Already have an error, skip processing
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
			throw std::runtime_error("MzMLReader: XML nesting too deep (>" + std::to_string(MAX_CONTEXT_DEPTH) +
			                         " levels) in " + filepath);
		}
		auto ctx = current_context();

		if (std::strcmp(name, "mzML") == 0 && (ctx == ParseContext::NONE || ctx == ParseContext::MZML)) {
			context_stack.push_back(ParseContext::MZML);
			return;
		}

		if (std::strcmp(name, "referenceableParamGroupList") == 0 && ctx == ParseContext::MZML) {
			context_stack.push_back(ParseContext::REFERENCEABLE_PARAM_GROUP_LIST);
			return;
		}

		if (std::strcmp(name, "referenceableParamGroup") == 0 && ctx == ParseContext::REFERENCEABLE_PARAM_GROUP_LIST) {
			context_stack.push_back(ParseContext::REFERENCEABLE_PARAM_GROUP);
			auto *id = get_attr(attrs, "id");
			current_param_group_id = id ? id : "";
			param_groups[current_param_group_id] = {};
			return;
		}

		if (std::strcmp(name, "cvParam") == 0 && ctx == ParseContext::REFERENCEABLE_PARAM_GROUP) {
			auto cv = parse_cv_param(attrs);
			param_groups[current_param_group_id].push_back(std::move(cv));
			return;
		}

		if (std::strcmp(name, "run") == 0 && ctx == ParseContext::MZML) {
			context_stack.push_back(ParseContext::RUN);
			return;
		}

		if (std::strcmp(name, "spectrumList") == 0 && ctx == ParseContext::RUN) {
			context_stack.push_back(ParseContext::SPECTRUM_LIST);
			has_spectrum_list = true;
			return;
		}

		if (std::strcmp(name, "spectrum") == 0 && ctx == ParseContext::SPECTRUM_LIST) {
			context_stack.push_back(ParseContext::SPECTRUM);
			current_spectrum = SpectrumState {};
			auto *idx = get_attr(attrs, "index");
			if (idx) {
				safe_stoi(std::string(idx), current_spectrum.index);
			}
			auto *id = get_attr(attrs, "id");
			if (id) {
				current_spectrum.id = id;
			}
			auto *dal = get_attr(attrs, "defaultArrayLength");
			if (dal) {
				safe_stoi(std::string(dal), current_spectrum.default_array_length);
			}
			return;
		}

		// Spectrum child elements
		if (ctx == ParseContext::SPECTRUM) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_spectrum(parse_cv_param(attrs));
				return;
			}
			if (std::strcmp(name, "scanList") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_SCAN_LIST);
				return;
			}
			if (std::strcmp(name, "precursorList") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_PRECURSOR_LIST);
				return;
			}
			if (std::strcmp(name, "binaryDataArrayList") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_BINARY_DATA_ARRAY_LIST);
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_SCAN_LIST) {
			if (std::strcmp(name, "scan") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_SCAN);
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_SCAN) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_scan(parse_cv_param(attrs));
				return;
			}
			if (std::strcmp(name, "scanWindow") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_SCAN_WINDOW);
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_SCAN_WINDOW) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_scan_window(parse_cv_param(attrs));
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_PRECURSOR_LIST) {
			if (std::strcmp(name, "precursor") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_PRECURSOR);
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_PRECURSOR) {
			if (std::strcmp(name, "isolationWindow") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_ISOLATION_WINDOW);
				return;
			}
			if (std::strcmp(name, "selectedIonList") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_SELECTED_ION_LIST);
				return;
			}
			if (std::strcmp(name, "activation") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_ACTIVATION);
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_ISOLATION_WINDOW) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_isolation_window(parse_cv_param(attrs));
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_SELECTED_ION_LIST) {
			if (std::strcmp(name, "selectedIon") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_SELECTED_ION);
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_SELECTED_ION) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_selected_ion(parse_cv_param(attrs));
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_ACTIVATION) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_activation(parse_cv_param(attrs));
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_BINARY_DATA_ARRAY_LIST) {
			if (std::strcmp(name, "binaryDataArray") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_BINARY_DATA_ARRAY);
				current_binary_array = BinaryDataArrayState {};
				return;
			}
		}

		if (ctx == ParseContext::SPECTRUM_BINARY_DATA_ARRAY) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_binary_array(parse_cv_param(attrs));
				return;
			}
			if (std::strcmp(name, "referenceableParamGroupRef") == 0) {
				auto *ref = get_attr(attrs, "ref");
				if (ref) {
					auto it = param_groups.find(ref);
					if (it != param_groups.end()) {
						for (const auto &cv : it->second) {
							apply_cv_param_to_binary_array(cv);
						}
					}
				}
				return;
			}
			if (std::strcmp(name, "binary") == 0) {
				context_stack.push_back(ParseContext::SPECTRUM_BINARY);
				current_binary_array.binary_text.clear();
				return;
			}
		}

		// Chromatogram elements
		if (std::strcmp(name, "chromatogramList") == 0 && ctx == ParseContext::RUN) {
			context_stack.push_back(ParseContext::CHROMATOGRAM_LIST);
			has_chromatogram_list = true;
			// Spectra are done when we enter chromatogram list
			spectra_done = true;
			return;
		}

		if (std::strcmp(name, "chromatogram") == 0 && ctx == ParseContext::CHROMATOGRAM_LIST) {
			context_stack.push_back(ParseContext::CHROMATOGRAM);
			current_chromatogram = ChromatogramState {};
			auto *idx = get_attr(attrs, "index");
			if (idx) {
				safe_stoi(std::string(idx), current_chromatogram.index);
			}
			auto *id = get_attr(attrs, "id");
			if (id) {
				current_chromatogram.id = id;
			}
			return;
		}

		if (ctx == ParseContext::CHROMATOGRAM) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_chromatogram(parse_cv_param(attrs));
				return;
			}
			if (std::strcmp(name, "precursor") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_PRECURSOR);
				return;
			}
			if (std::strcmp(name, "product") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_PRODUCT);
				return;
			}
			if (std::strcmp(name, "binaryDataArrayList") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_BINARY_DATA_ARRAY_LIST);
				return;
			}
		}

		if (ctx == ParseContext::CHROMATOGRAM_PRECURSOR) {
			if (std::strcmp(name, "isolationWindow") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_PRECURSOR_ISOLATION_WINDOW);
				return;
			}
		}

		if (ctx == ParseContext::CHROMATOGRAM_PRECURSOR_ISOLATION_WINDOW) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_chrom_precursor_isolation(parse_cv_param(attrs));
				return;
			}
		}

		if (ctx == ParseContext::CHROMATOGRAM_PRODUCT) {
			if (std::strcmp(name, "isolationWindow") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_PRODUCT_ISOLATION_WINDOW);
				return;
			}
		}

		if (ctx == ParseContext::CHROMATOGRAM_PRODUCT_ISOLATION_WINDOW) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_chrom_product_isolation(parse_cv_param(attrs));
				return;
			}
		}

		if (ctx == ParseContext::CHROMATOGRAM_BINARY_DATA_ARRAY_LIST) {
			if (std::strcmp(name, "binaryDataArray") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_BINARY_DATA_ARRAY);
				current_binary_array = BinaryDataArrayState {};
				return;
			}
		}

		if (ctx == ParseContext::CHROMATOGRAM_BINARY_DATA_ARRAY) {
			if (std::strcmp(name, "cvParam") == 0) {
				apply_cv_param_to_binary_array(parse_cv_param(attrs));
				return;
			}
			if (std::strcmp(name, "referenceableParamGroupRef") == 0) {
				auto *ref = get_attr(attrs, "ref");
				if (ref) {
					auto it = param_groups.find(ref);
					if (it != param_groups.end()) {
						for (const auto &cv : it->second) {
							apply_cv_param_to_binary_array(cv);
						}
					}
				}
				return;
			}
			if (std::strcmp(name, "binary") == 0) {
				context_stack.push_back(ParseContext::CHROMATOGRAM_BINARY);
				current_binary_array.binary_text.clear();
				return;
			}
		}
	}

	void on_end_element(const char *name) {
		auto ctx = current_context();

		// Pop context for matching end elements
		if (std::strcmp(name, "mzML") == 0 && ctx == ParseContext::MZML) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "referenceableParamGroupList") == 0 &&
		    ctx == ParseContext::REFERENCEABLE_PARAM_GROUP_LIST) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "referenceableParamGroup") == 0 && ctx == ParseContext::REFERENCEABLE_PARAM_GROUP) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "run") == 0 && ctx == ParseContext::RUN) {
			context_stack.pop_back();
			spectra_done = true;
			chromatograms_done = true;
			return;
		}
		if (std::strcmp(name, "spectrumList") == 0 && ctx == ParseContext::SPECTRUM_LIST) {
			context_stack.pop_back();
			spectra_done = true;
			return;
		}
		if (std::strcmp(name, "spectrum") == 0 && ctx == ParseContext::SPECTRUM) {
			// Record ms1_scan_index at parse time, before adding to completed queue
			if (current_spectrum.ms_level == 1) {
				last_ms1_index = current_spectrum.index;
				current_spectrum.ms1_scan_index = 0;
				current_spectrum.ms1_scan_index_valid = false;
			} else {
				if (last_ms1_index >= 0) {
					current_spectrum.ms1_scan_index = last_ms1_index;
					current_spectrum.ms1_scan_index_valid = true;
				} else {
					current_spectrum.ms1_scan_index = 0;
					current_spectrum.ms1_scan_index_valid = false;
				}
			}
			completed_spectra.push_back(std::move(current_spectrum));
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "scanList") == 0 && ctx == ParseContext::SPECTRUM_SCAN_LIST) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "scan") == 0 && ctx == ParseContext::SPECTRUM_SCAN) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "scanWindow") == 0 && ctx == ParseContext::SPECTRUM_SCAN_WINDOW) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "precursorList") == 0 && ctx == ParseContext::SPECTRUM_PRECURSOR_LIST) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "precursor") == 0 && ctx == ParseContext::SPECTRUM_PRECURSOR) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "isolationWindow") == 0 && ctx == ParseContext::SPECTRUM_ISOLATION_WINDOW) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "selectedIonList") == 0 && ctx == ParseContext::SPECTRUM_SELECTED_ION_LIST) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "selectedIon") == 0 && ctx == ParseContext::SPECTRUM_SELECTED_ION) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "activation") == 0 && ctx == ParseContext::SPECTRUM_ACTIVATION) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "binaryDataArrayList") == 0 && ctx == ParseContext::SPECTRUM_BINARY_DATA_ARRAY_LIST) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "binaryDataArray") == 0 && ctx == ParseContext::SPECTRUM_BINARY_DATA_ARRAY) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "binary") == 0 && ctx == ParseContext::SPECTRUM_BINARY) {
			// Only decode arrays we care about; skip non-standard arrays (e.g. "ms level")
			if (current_binary_array.is_mz || current_binary_array.is_intensity) {
				auto decoded =
				    MzMLBinaryDecoder::decode(current_binary_array.binary_text, current_binary_array.is_compressed,
				                              current_binary_array.is_64bit);
				if (current_binary_array.is_mz) {
					current_spectrum.mz_array = std::move(decoded);
				} else {
					current_spectrum.intensity_array = std::move(decoded);
				}
			}
			context_stack.pop_back();
			return;
		}

		// Chromatogram end elements
		if (std::strcmp(name, "chromatogramList") == 0 && ctx == ParseContext::CHROMATOGRAM_LIST) {
			context_stack.pop_back();
			chromatograms_done = true;
			return;
		}
		if (std::strcmp(name, "chromatogram") == 0 && ctx == ParseContext::CHROMATOGRAM) {
			completed_chromatograms.push_back(std::move(current_chromatogram));
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "precursor") == 0 && ctx == ParseContext::CHROMATOGRAM_PRECURSOR) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "isolationWindow") == 0 && ctx == ParseContext::CHROMATOGRAM_PRECURSOR_ISOLATION_WINDOW) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "product") == 0 && ctx == ParseContext::CHROMATOGRAM_PRODUCT) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "isolationWindow") == 0 && ctx == ParseContext::CHROMATOGRAM_PRODUCT_ISOLATION_WINDOW) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "binaryDataArrayList") == 0 && ctx == ParseContext::CHROMATOGRAM_BINARY_DATA_ARRAY_LIST) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "binaryDataArray") == 0 && ctx == ParseContext::CHROMATOGRAM_BINARY_DATA_ARRAY) {
			context_stack.pop_back();
			return;
		}
		if (std::strcmp(name, "binary") == 0 && ctx == ParseContext::CHROMATOGRAM_BINARY) {
			// Only decode arrays we care about; skip non-standard arrays (e.g. "ms level")
			if (current_binary_array.is_time || current_binary_array.is_intensity) {
				auto decoded =
				    MzMLBinaryDecoder::decode(current_binary_array.binary_text, current_binary_array.is_compressed,
				                              current_binary_array.is_64bit);
				if (current_binary_array.is_time) {
					current_chromatogram.time_array = std::move(decoded);
				} else {
					current_chromatogram.intensity_array = std::move(decoded);
				}
			}
			context_stack.pop_back();
			return;
		}
	}

	void on_char_data(const char *s, int len) {
		auto ctx = current_context();
		if (ctx == ParseContext::SPECTRUM_BINARY || ctx == ParseContext::CHROMATOGRAM_BINARY) {
			current_binary_array.binary_text.append(s, len);
		}
	}
};

// MzMLReader public API

MzMLReader::MzMLReader(const std::string &path) : impl_(std::make_unique<Impl>(path)) {
}

MzMLReader::~MzMLReader() = default;

MzMLReader::MzMLReader(MzMLReader &&) noexcept = default;
MzMLReader &MzMLReader::operator=(MzMLReader &&) noexcept = default;

MzMLSpectrumBatch MzMLReader::read_spectra(size_t n) {
	MzMLSpectrumBatch batch;
	if (n == 0) {
		return batch;
	}

	// Parse more if needed and not done
	if (impl_->completed_spectra.size() < n && !impl_->spectra_done) {
		impl_->in_chromatogram_mode = false;
		impl_->parse_until(n, 0);
	}

	// Move completed spectra into batch (up to n)
	size_t count = std::min(n, impl_->completed_spectra.size());

	// Pre-allocate all batch vectors to avoid geometric reallocation
	batch.spectrum_index.reserve(count);
	batch.spectrum_id.reserve(count);
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
		auto &spec = impl_->completed_spectra.front();
		batch.spectrum_index.push_back(spec.index);
		batch.spectrum_id.push_back(std::move(spec.id));
		batch.ms_level.push_back(spec.ms_level);
		batch.retention_time.push_back(spec.retention_time);
		batch.retention_time_valid.push_back(spec.retention_time_valid);
		batch.spectrum_type.push_back(std::move(spec.spectrum_type));
		batch.polarity.push_back(std::move(spec.polarity));
		batch.base_peak_mz.push_back(spec.base_peak_mz);
		batch.base_peak_mz_valid.push_back(spec.base_peak_mz_valid);
		batch.base_peak_intensity.push_back(spec.base_peak_intensity);
		batch.base_peak_intensity_valid.push_back(spec.base_peak_intensity_valid);
		batch.total_ion_current.push_back(spec.total_ion_current);
		batch.total_ion_current_valid.push_back(spec.total_ion_current_valid);
		batch.lowest_mz.push_back(spec.lowest_mz);
		batch.lowest_mz_valid.push_back(spec.lowest_mz_valid);
		batch.highest_mz.push_back(spec.highest_mz);
		batch.highest_mz_valid.push_back(spec.highest_mz_valid);
		batch.default_array_length.push_back(spec.default_array_length);
		batch.precursor_mz.push_back(spec.precursor_mz);
		batch.precursor_mz_valid.push_back(spec.precursor_mz_valid);
		batch.precursor_charge.push_back(spec.precursor_charge);
		batch.precursor_charge_valid.push_back(spec.precursor_charge_valid);
		batch.precursor_intensity.push_back(spec.precursor_intensity);
		batch.precursor_intensity_valid.push_back(spec.precursor_intensity_valid);
		batch.isolation_window_target.push_back(spec.isolation_window_target);
		batch.isolation_window_target_valid.push_back(spec.isolation_window_target_valid);
		batch.isolation_window_lower.push_back(spec.isolation_window_lower);
		batch.isolation_window_lower_valid.push_back(spec.isolation_window_lower_valid);
		batch.isolation_window_upper.push_back(spec.isolation_window_upper);
		batch.isolation_window_upper_valid.push_back(spec.isolation_window_upper_valid);
		batch.activation_method.push_back(std::move(spec.activation_method));
		batch.collision_energy.push_back(spec.collision_energy);
		batch.collision_energy_valid.push_back(spec.collision_energy_valid);
		batch.mz_array.push_back(std::move(spec.mz_array));
		batch.intensity_array.push_back(std::move(spec.intensity_array));
		batch.filter_string.push_back(std::move(spec.filter_string));
		batch.scan_window_lower.push_back(spec.scan_window_lower);
		batch.scan_window_lower_valid.push_back(spec.scan_window_lower_valid);
		batch.scan_window_upper.push_back(spec.scan_window_upper);
		batch.scan_window_upper_valid.push_back(spec.scan_window_upper_valid);
		batch.ms1_scan_index.push_back(spec.ms1_scan_index);
		batch.ms1_scan_index_valid.push_back(spec.ms1_scan_index_valid);

		impl_->completed_spectra.pop_front();
	}

	return batch;
}

MzMLChromatogramBatch MzMLReader::read_chromatograms(size_t n) {
	MzMLChromatogramBatch batch;
	if (n == 0) {
		return batch;
	}

	// Fast-forward past all spectra. Clear after each parse chunk to avoid
	// accumulating unbounded spectrum data (binary arrays can be huge).
	if (!impl_->spectra_done) {
		impl_->in_chromatogram_mode = false;
		while (!impl_->spectra_done && !impl_->parse_done) {
			impl_->parse_until(SIZE_MAX, 0);
			impl_->completed_spectra.clear();
		}
	}

	// Parse more chromatograms if needed
	if (impl_->completed_chromatograms.size() < n && !impl_->chromatograms_done) {
		impl_->in_chromatogram_mode = true;
		impl_->parse_until(0, n);
	}

	size_t count = std::min(n, impl_->completed_chromatograms.size());

	// Pre-allocate all batch vectors to avoid geometric reallocation
	batch.chromatogram_index.reserve(count);
	batch.chromatogram_id.reserve(count);
	batch.chromatogram_type.reserve(count);
	batch.precursor_mz.reserve(count);
	batch.precursor_mz_valid.reserve(count);
	batch.product_mz.reserve(count);
	batch.product_mz_valid.reserve(count);
	batch.time_array.reserve(count);
	batch.intensity_array.reserve(count);

	for (size_t i = 0; i < count; i++) {
		auto &chrom = impl_->completed_chromatograms.front();
		batch.chromatogram_index.push_back(chrom.index);
		batch.chromatogram_id.push_back(std::move(chrom.id));
		batch.chromatogram_type.push_back(std::move(chrom.chromatogram_type));
		batch.precursor_mz.push_back(chrom.precursor_mz);
		batch.precursor_mz_valid.push_back(chrom.precursor_mz_valid);
		batch.product_mz.push_back(chrom.product_mz);
		batch.product_mz_valid.push_back(chrom.product_mz_valid);
		batch.time_array.push_back(std::move(chrom.time_array));
		batch.intensity_array.push_back(std::move(chrom.intensity_array));

		impl_->completed_chromatograms.pop_front();
	}

	return batch;
}

bool MzMLReader::has_spectra() const {
	return impl_->has_spectrum_list;
}

bool MzMLReader::has_chromatograms() const {
	return impl_->has_chromatogram_list;
}

} // namespace miint
