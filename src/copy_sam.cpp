#include "copy_sam.hpp"
#include "reference_table_reader.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include <htslib-1.22.1/htslib/sam.h>
#include <htslib-1.22.1/htslib/hts.h>
#include <unordered_map>
#include <unordered_set>

namespace duckdb {

//===--------------------------------------------------------------------===//
// HTSlib Smart Pointers
//===--------------------------------------------------------------------===//
struct SAMFileDeleter {
	void operator()(samFile *fp) const {
		if (fp) {
			sam_close(fp);
		}
	}
};

struct SAMHeaderDeleter {
	void operator()(sam_hdr_t *hdr) const {
		if (hdr) {
			sam_hdr_destroy(hdr);
		}
	}
};

struct BAMRecordDeleter {
	void operator()(bam1_t *aln) const {
		if (aln) {
			bam_destroy1(aln);
		}
	}
};

using SAMFilePtr = std::unique_ptr<samFile, SAMFileDeleter>;
using SAMHeaderPtr = std::unique_ptr<sam_hdr_t, SAMHeaderDeleter>;
using BAMRecordPtr = std::unique_ptr<bam1_t, BAMRecordDeleter>;

//===--------------------------------------------------------------------===//
// Column Indices for SAM
//===--------------------------------------------------------------------===//
struct SAMColumnIndices {
	idx_t read_id_idx = DConstants::INVALID_INDEX;
	idx_t flags_idx = DConstants::INVALID_INDEX;
	idx_t reference_idx = DConstants::INVALID_INDEX;
	idx_t position_idx = DConstants::INVALID_INDEX;
	idx_t stop_position_idx = DConstants::INVALID_INDEX;
	idx_t mapq_idx = DConstants::INVALID_INDEX;
	idx_t cigar_idx = DConstants::INVALID_INDEX;
	idx_t mate_reference_idx = DConstants::INVALID_INDEX;
	idx_t mate_position_idx = DConstants::INVALID_INDEX;
	idx_t template_length_idx = DConstants::INVALID_INDEX;

	// Tag columns
	idx_t tag_as_idx = DConstants::INVALID_INDEX;
	idx_t tag_xs_idx = DConstants::INVALID_INDEX;
	idx_t tag_ys_idx = DConstants::INVALID_INDEX;
	idx_t tag_xn_idx = DConstants::INVALID_INDEX;
	idx_t tag_xm_idx = DConstants::INVALID_INDEX;
	idx_t tag_xo_idx = DConstants::INVALID_INDEX;
	idx_t tag_xg_idx = DConstants::INVALID_INDEX;
	idx_t tag_nm_idx = DConstants::INVALID_INDEX;
	idx_t tag_yt_idx = DConstants::INVALID_INDEX;
	idx_t tag_md_idx = DConstants::INVALID_INDEX;
	idx_t tag_sa_idx = DConstants::INVALID_INDEX;

	void FindIndices(const vector<string> &names) {
		for (idx_t i = 0; i < names.size(); i++) {
			auto &name = names[i];
			if (StringUtil::CIEquals(name, "read_id")) {
				read_id_idx = i;
			} else if (StringUtil::CIEquals(name, "flags")) {
				flags_idx = i;
			} else if (StringUtil::CIEquals(name, "reference")) {
				reference_idx = i;
			} else if (StringUtil::CIEquals(name, "position")) {
				position_idx = i;
			} else if (StringUtil::CIEquals(name, "stop_position")) {
				stop_position_idx = i;
			} else if (StringUtil::CIEquals(name, "mapq")) {
				mapq_idx = i;
			} else if (StringUtil::CIEquals(name, "cigar")) {
				cigar_idx = i;
			} else if (StringUtil::CIEquals(name, "mate_reference")) {
				mate_reference_idx = i;
			} else if (StringUtil::CIEquals(name, "mate_position")) {
				mate_position_idx = i;
			} else if (StringUtil::CIEquals(name, "template_length")) {
				template_length_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_as")) {
				tag_as_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_xs")) {
				tag_xs_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_ys")) {
				tag_ys_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_xn")) {
				tag_xn_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_xm")) {
				tag_xm_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_xo")) {
				tag_xo_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_xg")) {
				tag_xg_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_nm")) {
				tag_nm_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_yt")) {
				tag_yt_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_md")) {
				tag_md_idx = i;
			} else if (StringUtil::CIEquals(name, "tag_sa")) {
				tag_sa_idx = i;
			}
		}
	}
};

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
enum class SAMOutputFormat : uint8_t { SAM = 0, BAM = 1 };

struct SAMCopyBindData : public FunctionData {
	bool include_header = true;
	bool use_gzip = false;
	SAMOutputFormat format = SAMOutputFormat::SAM;
	int compression_level = -1; // -1 means use HTSlib default (6 for BAM)
	string file_path;
	vector<string> names;
	std::optional<string> reference_lengths_table;
	SAMColumnIndices indices; // Cached column indices

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<SAMCopyBindData>();
		result->include_header = include_header;
		result->use_gzip = use_gzip;
		result->format = format;
		result->compression_level = compression_level;
		result->file_path = file_path;
		result->names = names;
		result->reference_lengths_table = reference_lengths_table;
		result->indices = indices;
		return std::move(result);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<SAMCopyBindData>();
		return include_header == other.include_header && use_gzip == other.use_gzip && format == other.format &&
		       compression_level == other.compression_level && file_path == other.file_path &&
		       reference_lengths_table == other.reference_lengths_table;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> SAMCopyBindInternal(ClientContext &context, CopyFunctionBindInput &input,
                                                    const vector<string> &names, const vector<LogicalType> &sql_types,
                                                    SAMOutputFormat default_format) {
	auto result = make_uniq<SAMCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;
	result->format = default_format;

	// Detect and cache column indices
	result->indices.FindIndices(names);
	auto &indices = result->indices;

	// Validate required columns exist
	if (indices.read_id_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'read_id' column");
	}
	if (indices.flags_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'flags' column");
	}
	if (indices.reference_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'reference' column");
	}
	if (indices.position_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'position' column");
	}
	if (indices.mapq_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'mapq' column");
	}
	if (indices.cigar_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'cigar' column");
	}
	if (indices.mate_reference_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'mate_reference' column");
	}
	if (indices.mate_position_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'mate_position' column");
	}
	if (indices.template_length_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT SAM requires 'template_length' column");
	}

	// Validate column types
	if (sql_types[indices.read_id_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'read_id' must be VARCHAR");
	}
	if (sql_types[indices.flags_idx].id() != LogicalTypeId::USMALLINT) {
		throw BinderException("Column 'flags' must be USMALLINT");
	}
	if (sql_types[indices.reference_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'reference' must be VARCHAR");
	}
	if (sql_types[indices.position_idx].id() != LogicalTypeId::BIGINT) {
		throw BinderException("Column 'position' must be BIGINT");
	}
	if (sql_types[indices.mapq_idx].id() != LogicalTypeId::UTINYINT) {
		throw BinderException("Column 'mapq' must be UTINYINT");
	}
	if (sql_types[indices.cigar_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'cigar' must be VARCHAR");
	}
	if (sql_types[indices.mate_reference_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'mate_reference' must be VARCHAR");
	}
	if (sql_types[indices.mate_position_idx].id() != LogicalTypeId::BIGINT) {
		throw BinderException("Column 'mate_position' must be BIGINT");
	}
	if (sql_types[indices.template_length_idx].id() != LogicalTypeId::BIGINT) {
		throw BinderException("Column 'template_length' must be BIGINT");
	}

	// Parse options
	bool compression_specified = false;
	for (auto &option : input.info.options) {
		if (StringUtil::CIEquals(option.first, "include_header")) {
			result->include_header = option.second[0].GetValue<bool>();
		} else if (StringUtil::CIEquals(option.first, "compression")) {
			compression_specified = true;
			auto comp_value = option.second[0].ToString();
			if (StringUtil::CIEquals(comp_value, "gzip") || StringUtil::CIEquals(comp_value, "gz")) {
				result->use_gzip = true;
			} else if (StringUtil::CIEquals(comp_value, "none") || StringUtil::CIEquals(comp_value, "uncompressed")) {
				result->use_gzip = false;
			} else {
				throw BinderException("Unknown compression type for COPY FORMAT SAM: %s (supported: gzip, none)",
				                      comp_value);
			}
		} else if (StringUtil::CIEquals(option.first, "compression_level")) {
			result->compression_level = option.second[0].GetValue<int32_t>();
			if (result->compression_level < 0 || result->compression_level > 9) {
				throw BinderException("COMPRESSION_LEVEL must be between 0 and 9, got %d", result->compression_level);
			}
		} else if (StringUtil::CIEquals(option.first, "reference_lengths")) {
			const auto &table_value = option.second[0];
			if (table_value.type().id() != LogicalTypeId::VARCHAR) {
				throw BinderException("reference_lengths must be a VARCHAR (table or view name)");
			}

			result->reference_lengths_table = table_value.ToString();

			// Validate table or view exists (use TABLE_ENTRY lookup which returns either)
			EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, result->reference_lengths_table.value(), QueryErrorContext());
			auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);
			if (!entry) {
				throw BinderException("Table or view '%s' does not exist", result->reference_lengths_table.value());
			}
			if (entry->type != CatalogType::TABLE_ENTRY && entry->type != CatalogType::VIEW_ENTRY) {
				throw BinderException("'%s' is not a table or view", result->reference_lengths_table.value());
			}
		} else {
			throw BinderException("Unknown option for COPY FORMAT SAM: %s", option.first);
		}
	}

	// Validate COMPRESSION_LEVEL only makes sense for BAM
	if (result->compression_level != -1 && result->format != SAMOutputFormat::BAM) {
		throw BinderException("COMPRESSION_LEVEL parameter only applies to BAM format");
	}

	// BAM files must have headers (binary format requirement)
	if (result->format == SAMOutputFormat::BAM && !result->include_header) {
		throw BinderException("BAM format requires INCLUDE_HEADER=true (BAM is a binary format that requires headers)");
	}

	// Validate reference_lengths requirement
	// For BAM format, reference_lengths is always required (BAM needs complete header upfront)
	// For SAM format with INCLUDE_HEADER=true, reference_lengths is required to build the header
	// For SAM format with INCLUDE_HEADER=false, reference_lengths is optional
	if (result->format == SAMOutputFormat::BAM && !result->reference_lengths_table.has_value()) {
		throw BinderException("COPY FORMAT BAM requires REFERENCE_LENGTHS parameter "
		                      "(e.g., REFERENCE_LENGTHS='ref_table')");
	}
	if (result->include_header && !result->reference_lengths_table.has_value()) {
		throw BinderException("COPY FORMAT SAM with INCLUDE_HEADER=true requires REFERENCE_LENGTHS parameter "
		                      "(e.g., REFERENCE_LENGTHS='ref_table')");
	}

	// Auto-detect compression from file extension if not specified (for SAM format)
	if (!compression_specified && result->format == SAMOutputFormat::SAM) {
		if (result->file_path.size() >= 3 && result->file_path.substr(result->file_path.size() - 3) == ".gz") {
			result->use_gzip = true;
		}
	}

	return std::move(result);
}

// Bind function for SAM format
static unique_ptr<FunctionData> SAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                            const vector<string> &names, const vector<LogicalType> &sql_types) {
	return SAMCopyBindInternal(context, input, names, sql_types, SAMOutputFormat::SAM);
}

// Bind function for BAM format
static unique_ptr<FunctionData> BAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                            const vector<string> &names, const vector<LogicalType> &sql_types) {
	return SAMCopyBindInternal(context, input, names, sql_types, SAMOutputFormat::BAM);
}

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct SAMCopyGlobalState : public GlobalFunctionData {
	mutex write_lock;
	SAMFilePtr sam_file;
	SAMHeaderPtr header;
	std::atomic<bool> header_written{false};
	bool include_header = true;
};

static unique_ptr<GlobalFunctionData> SAMCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                              const string &file_path) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();

	auto gstate = make_uniq<SAMCopyGlobalState>();
	gstate->include_header = fdata.include_header;

	// Determine write mode based on format and compression
	string mode;
	if (fdata.format == SAMOutputFormat::BAM) {
		// BAM format
		if (fdata.compression_level >= 0) {
			// Explicit compression level (0-9)
			mode = "wb" + std::to_string(fdata.compression_level);
		} else {
			// Default BAM compression (HTSlib default level 6)
			mode = "wb";
		}
	} else {
		// SAM format
		mode = fdata.use_gzip ? "wz" : "w";
	}

	// Open file for writing
	gstate->sam_file = SAMFilePtr(sam_open(file_path.c_str(), mode.c_str()));
	if (!gstate->sam_file) {
		throw IOException("Failed to open file for writing: " + file_path);
	}

	// Create header from reference_lengths if provided, otherwise create empty header
	gstate->header = SAMHeaderPtr(sam_hdr_init());
	if (!gstate->header) {
		throw std::runtime_error("Failed to create SAM header");
	}

	if (fdata.reference_lengths_table.has_value()) {
		auto reference_lengths_map = ReadReferenceTable(context, fdata.reference_lengths_table.value());

		// Add each reference to header
		for (const auto &ref : reference_lengths_map) {
			if (sam_hdr_add_line(gstate->header.get(), "SQ", "SN", ref.first.c_str(), "LN",
			                     std::to_string(ref.second).c_str(), NULL) < 0) {
				throw std::runtime_error("Failed to add reference to SAM header: " + ref.first);
			}
		}
	}
	// If no reference_lengths provided, header remains empty
	// References will be added dynamically as we encounter them during writing

	return std::move(gstate);
}

//===--------------------------------------------------------------------===//
// Local State
//===--------------------------------------------------------------------===//
struct SAMCopyLocalState : public LocalFunctionData {
	BAMRecordPtr record;
	std::vector<uint32_t> cigar_buffer;
	size_t cigar_buffer_capacity = 0;

	SAMCopyLocalState() {
		record = BAMRecordPtr(bam_init1());
		if (!record) {
			throw std::runtime_error("Failed to allocate BAM record");
		}
	}
};

static unique_ptr<LocalFunctionData> SAMCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	return make_uniq<SAMCopyLocalState>();
}

//===--------------------------------------------------------------------===//
// Helper Functions
//===--------------------------------------------------------------------===//

// Parse CIGAR string and store in buffer
static bool ParseCIGAR(const string &cigar_str, std::vector<uint32_t> &cigar_buffer) {
	cigar_buffer.clear();

	// Handle "*" CIGAR (no alignment)
	if (cigar_str == "*") {
		return true;
	}

	// Use HTSlib's CIGAR parser
	char *end = nullptr;
	uint32_t *cigar_array = nullptr;
	size_t cigar_mem = 0;

	ssize_t n_cigar = sam_parse_cigar(cigar_str.c_str(), &end, &cigar_array, &cigar_mem);

	if (n_cigar < 0 || (end && *end != '\0')) {
		if (cigar_array) {
			free(cigar_array);
		}
		return false;
	}

	// Copy to buffer
	cigar_buffer.assign(cigar_array, cigar_array + n_cigar);
	free(cigar_array);

	return true;
}

// Get reference TID from header, adding it dynamically if not present
// When writing headerless SAM files, we add references with a sentinel length
static int32_t GetReferenceTID(sam_hdr_t *header, const string &reference) {
	if (reference == "*") {
		return -1;
	}
	int32_t tid = sam_hdr_name2tid(header, reference.c_str());

	// If reference not in header, add it with a sentinel length (2^31-1)
	// This allows writing headerless SAM files without pre-defined reference lengths
	if (tid < 0) {
		constexpr const char *SENTINEL_LENGTH = "2147483647"; // 2^31-1, unknown length
		if (sam_hdr_add_line(header, "SQ", "SN", reference.c_str(), "LN", SENTINEL_LENGTH, NULL) < 0) {
			throw std::runtime_error("Failed to dynamically add reference to SAM header: " + reference);
		}
		// Get the tid after adding
		tid = sam_hdr_name2tid(header, reference.c_str());
		if (tid < 0) {
			throw std::runtime_error("Failed to get TID for dynamically added reference: " + reference);
		}
	}

	return tid;
}

//===--------------------------------------------------------------------===//
// Tag Helpers
//===--------------------------------------------------------------------===//
static inline void AppendIntegerTag(bam1_t *record, const char *tag_name, const int64_t *data_ptr,
                                    const UnifiedVectorFormat &format, idx_t row_idx, const string &read_id) {
	if (!data_ptr) {
		return;
	}
	auto data_idx = format.sel->get_index(row_idx);
	if (format.validity.RowIsValid(data_idx)) {
		int32_t value = static_cast<int32_t>(data_ptr[data_idx]);
		if (bam_aux_append(record, tag_name, 'i', sizeof(int32_t), reinterpret_cast<uint8_t *>(&value)) < 0) {
			throw IOException("Failed to append %s tag to BAM record for read: %s", tag_name, read_id);
		}
	}
}

static inline void AppendStringTag(bam1_t *record, const char *tag_name, const string_t *data_ptr,
                                   const UnifiedVectorFormat &format, idx_t row_idx, const string &read_id) {
	if (!data_ptr) {
		return;
	}
	auto data_idx = format.sel->get_index(row_idx);
	if (format.validity.RowIsValid(data_idx)) {
		string value = data_ptr[data_idx].GetString();
		if (bam_aux_append(record, tag_name, 'Z', (int)value.size() + 1,
		                   reinterpret_cast<const uint8_t *>(value.c_str())) < 0) {
			throw IOException("Failed to append %s tag to BAM record for read: %s", tag_name, read_id);
		}
	}
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void SAMCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                        LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();
	auto &lstate = lstate_p.Cast<SAMCopyLocalState>();

	// Use cached column indices from bind data
	auto &indices = fdata.indices;

	// Extract data from input vectors
	UnifiedVectorFormat read_id_data, flags_data, reference_data, position_data, mapq_data;
	UnifiedVectorFormat cigar_data, mate_reference_data, mate_position_data, template_length_data;
	UnifiedVectorFormat tag_as_data, tag_xs_data, tag_ys_data, tag_xn_data, tag_xm_data, tag_xo_data, tag_xg_data,
	    tag_nm_data;
	UnifiedVectorFormat tag_yt_data, tag_md_data, tag_sa_data;

	input.data[indices.read_id_idx].ToUnifiedFormat(input.size(), read_id_data);
	input.data[indices.flags_idx].ToUnifiedFormat(input.size(), flags_data);
	input.data[indices.reference_idx].ToUnifiedFormat(input.size(), reference_data);
	input.data[indices.position_idx].ToUnifiedFormat(input.size(), position_data);
	input.data[indices.mapq_idx].ToUnifiedFormat(input.size(), mapq_data);
	input.data[indices.cigar_idx].ToUnifiedFormat(input.size(), cigar_data);
	input.data[indices.mate_reference_idx].ToUnifiedFormat(input.size(), mate_reference_data);
	input.data[indices.mate_position_idx].ToUnifiedFormat(input.size(), mate_position_data);
	input.data[indices.template_length_idx].ToUnifiedFormat(input.size(), template_length_data);

	if (indices.tag_as_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_as_idx].ToUnifiedFormat(input.size(), tag_as_data);
	}
	if (indices.tag_xs_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_xs_idx].ToUnifiedFormat(input.size(), tag_xs_data);
	}
	if (indices.tag_ys_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_ys_idx].ToUnifiedFormat(input.size(), tag_ys_data);
	}
	if (indices.tag_xn_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_xn_idx].ToUnifiedFormat(input.size(), tag_xn_data);
	}
	if (indices.tag_xm_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_xm_idx].ToUnifiedFormat(input.size(), tag_xm_data);
	}
	if (indices.tag_xo_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_xo_idx].ToUnifiedFormat(input.size(), tag_xo_data);
	}
	if (indices.tag_xg_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_xg_idx].ToUnifiedFormat(input.size(), tag_xg_data);
	}
	if (indices.tag_nm_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_nm_idx].ToUnifiedFormat(input.size(), tag_nm_data);
	}
	if (indices.tag_yt_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_yt_idx].ToUnifiedFormat(input.size(), tag_yt_data);
	}
	if (indices.tag_md_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_md_idx].ToUnifiedFormat(input.size(), tag_md_data);
	}
	if (indices.tag_sa_idx != DConstants::INVALID_INDEX) {
		input.data[indices.tag_sa_idx].ToUnifiedFormat(input.size(), tag_sa_data);
	}

	auto read_id_ptr = UnifiedVectorFormat::GetData<string_t>(read_id_data);
	auto flags_ptr = UnifiedVectorFormat::GetData<uint16_t>(flags_data);
	auto reference_ptr = UnifiedVectorFormat::GetData<string_t>(reference_data);
	auto position_ptr = UnifiedVectorFormat::GetData<int64_t>(position_data);
	auto mapq_ptr = UnifiedVectorFormat::GetData<uint8_t>(mapq_data);
	auto cigar_ptr = UnifiedVectorFormat::GetData<string_t>(cigar_data);
	auto mate_reference_ptr = UnifiedVectorFormat::GetData<string_t>(mate_reference_data);
	auto mate_position_ptr = UnifiedVectorFormat::GetData<int64_t>(mate_position_data);
	auto template_length_ptr = UnifiedVectorFormat::GetData<int64_t>(template_length_data);

	auto tag_as_ptr =
	    indices.tag_as_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_as_data) : nullptr;
	auto tag_xs_ptr =
	    indices.tag_xs_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_xs_data) : nullptr;
	auto tag_ys_ptr =
	    indices.tag_ys_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_ys_data) : nullptr;
	auto tag_xn_ptr =
	    indices.tag_xn_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_xn_data) : nullptr;
	auto tag_xm_ptr =
	    indices.tag_xm_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_xm_data) : nullptr;
	auto tag_xo_ptr =
	    indices.tag_xo_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_xo_data) : nullptr;
	auto tag_xg_ptr =
	    indices.tag_xg_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_xg_data) : nullptr;
	auto tag_nm_ptr =
	    indices.tag_nm_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<int64_t>(tag_nm_data) : nullptr;
	auto tag_yt_ptr =
	    indices.tag_yt_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<string_t>(tag_yt_data) : nullptr;
	auto tag_md_ptr =
	    indices.tag_md_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<string_t>(tag_md_data) : nullptr;
	auto tag_sa_ptr =
	    indices.tag_sa_idx != DConstants::INVALID_INDEX ? UnifiedVectorFormat::GetData<string_t>(tag_sa_data) : nullptr;

	// Write header once (first thread to arrive)
	if (!gstate.header_written.load()) {
		lock_guard<mutex> glock(gstate.write_lock);
		if (!gstate.header_written.load()) {
			if (fdata.include_header) {
				if (sam_hdr_write(gstate.sam_file.get(), gstate.header.get()) < 0) {
					throw IOException("Failed to write SAM header");
				}
			}
			gstate.header_written.store(true);
		}
	}

	// Process each row
	for (idx_t i = 0; i < input.size(); i++) {
		auto read_id_idx = read_id_data.sel->get_index(i);
		auto flags_idx = flags_data.sel->get_index(i);
		auto reference_idx = reference_data.sel->get_index(i);
		auto position_idx = position_data.sel->get_index(i);
		auto mapq_idx = mapq_data.sel->get_index(i);
		auto cigar_idx = cigar_data.sel->get_index(i);
		auto mate_reference_idx = mate_reference_data.sel->get_index(i);
		auto mate_position_idx = mate_position_data.sel->get_index(i);
		auto template_length_idx = template_length_data.sel->get_index(i);

		string read_id = read_id_ptr[read_id_idx].GetString();
		uint16_t flags = flags_ptr[flags_idx];
		string reference = reference_ptr[reference_idx].GetString();
		int64_t position = position_ptr[position_idx];
		uint8_t mapq = mapq_ptr[mapq_idx];
		string cigar_str = cigar_ptr[cigar_idx].GetString();
		string mate_reference = mate_reference_ptr[mate_reference_idx].GetString();
		int64_t mate_position = mate_position_ptr[mate_position_idx];
		int64_t template_length = template_length_ptr[template_length_idx];

		// Parse CIGAR
		if (!ParseCIGAR(cigar_str, lstate.cigar_buffer)) {
			throw InvalidInputException("Failed to parse CIGAR string '%s' for read: %s", cigar_str, read_id);
		}

		// Validate position values before casting
		if (position < 0) {
			throw InvalidInputException("Invalid position value %lld for read: %s (must be non-negative)", position, read_id);
		}
		if (position > 0 && (position - 1) > HTS_POS_MAX) {
			throw InvalidInputException("Position value %lld exceeds maximum allowed position for read: %s", position, read_id);
		}
		if (mate_position < 0) {
			throw InvalidInputException("Invalid mate_position value %lld for read: %s (must be non-negative)", mate_position, read_id);
		}
		if (mate_position > 0 && (mate_position - 1) > HTS_POS_MAX) {
			throw InvalidInputException("Mate position value %lld exceeds maximum allowed position for read: %s", mate_position, read_id);
		}

		// Get reference TIDs
		int32_t tid = GetReferenceTID(gstate.header.get(), reference);

		// Handle mate_reference: "=" means same as reference
		int32_t mtid;
		if (mate_reference == "=") {
			mtid = tid;
		} else {
			mtid = GetReferenceTID(gstate.header.get(), mate_reference);
		}

		// Build bam1_t record
		// Note: SEQ and QUAL are "*" (represented as NULL/0-length in bam_set1)
		// l_qname is the length WITHOUT null terminator (bam_set1 adds it internally)
		if (bam_set1(lstate.record.get(), read_id.length(), read_id.c_str(), flags, tid,
		             static_cast<hts_pos_t>(position - 1), // Convert to 0-based
		             mapq, lstate.cigar_buffer.size(), lstate.cigar_buffer.data(), mtid,
		             static_cast<hts_pos_t>(mate_position - 1), // Convert to 0-based
		             static_cast<hts_pos_t>(template_length), 0, "*", "*", 0) < 0) {
			throw IOException("Failed to build BAM record for read: " + read_id);
		}

		// Add optional tags
		AppendIntegerTag(lstate.record.get(), "AS", tag_as_ptr, tag_as_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "XS", tag_xs_ptr, tag_xs_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "YS", tag_ys_ptr, tag_ys_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "XN", tag_xn_ptr, tag_xn_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "XM", tag_xm_ptr, tag_xm_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "XO", tag_xo_ptr, tag_xo_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "XG", tag_xg_ptr, tag_xg_data, i, read_id);
		AppendIntegerTag(lstate.record.get(), "NM", tag_nm_ptr, tag_nm_data, i, read_id);
		AppendStringTag(lstate.record.get(), "YT", tag_yt_ptr, tag_yt_data, i, read_id);
		AppendStringTag(lstate.record.get(), "MD", tag_md_ptr, tag_md_data, i, read_id);
		AppendStringTag(lstate.record.get(), "SA", tag_sa_ptr, tag_sa_data, i, read_id);

		// Write record (with lock for thread safety on file I/O)
		{
			lock_guard<mutex> glock(gstate.write_lock);
			if (sam_write1(gstate.sam_file.get(), gstate.header.get(), lstate.record.get()) < 0) {
				throw IOException("Failed to write SAM record for read: " + read_id);
			}
		}
	}
}

//===--------------------------------------------------------------------===//
// Combine
//===--------------------------------------------------------------------===//
static void SAMCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                           LocalFunctionData &lstate_p) {
	// HTSlib handles buffering internally, nothing to combine
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void SAMCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();

	// Write header if we haven't written anything yet (empty file with header)
	if (!gstate.header_written.load() && fdata.include_header) {
		if (sam_hdr_write(gstate.sam_file.get(), gstate.header.get()) < 0) {
			throw IOException("Failed to write SAM header for empty file");
		}
	}

	// File is closed automatically by smart pointer destructor
}

//===--------------------------------------------------------------------===//
// Function Definition
//===--------------------------------------------------------------------===//
CopyFunction CopySAMFunction::GetFunction() {
	CopyFunction function("sam");
	function.copy_to_bind = SAMCopyBind;
	function.copy_to_initialize_global = SAMCopyInitializeGlobal;
	function.copy_to_initialize_local = SAMCopyInitializeLocal;
	function.copy_to_sink = SAMCopySink;
	function.copy_to_combine = SAMCopyCombine;
	function.copy_to_finalize = SAMCopyFinalize;
	function.extension = "sam";

	return function;
}

CopyFunction CopySAMFunction::GetBAMFunction() {
	CopyFunction function("bam");
	function.copy_to_bind = BAMCopyBind;
	function.copy_to_initialize_global = SAMCopyInitializeGlobal;
	function.copy_to_initialize_local = SAMCopyInitializeLocal;
	function.copy_to_sink = SAMCopySink;
	function.copy_to_combine = SAMCopyCombine;
	function.copy_to_finalize = SAMCopyFinalize;
	function.extension = "bam";

	return function;
}

void CopySAMFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
	loader.RegisterFunction(GetBAMFunction());
}

} // namespace duckdb
