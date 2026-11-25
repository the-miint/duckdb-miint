#include "copy_sam.hpp"
#include "copy_format_common.hpp"
#include "copy_option_parsing.hpp"
#include "htslib_raii.hpp"
#include "id_column_codec.hpp"
#include "id_column_utils.hpp"
#include "reference_table_reader.hpp"
#include "sequence_data_reader.hpp"
#include "sequence_utils.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/types/hash.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/client_context_state.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/storage/arena_allocator.hpp"
#include <htslib-1.22.1/htslib/sam.h>
#include <htslib-1.22.1/htslib/hts.h>
#include <htslib-1.22.1/htslib/bgzf.h>
#include <htslib-1.22.1/htslib/hfile.h>
#include <cerrno>
#include <cstring>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace duckdb {

// HTSlib RAII smart pointers (SAMFilePtr / SAMHeaderPtr / BAMRecordPtr) live in
// htslib_raii.hpp, shared with copy_ubam.cpp.

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

// The tag columns, in the order their tags are appended to a record: the order is in the output's bytes. An integer
// tag column must be BIGINT and a string one VARCHAR (checked at bind).
struct SAMTagColumn {
	const char *column;
	const char *tag;
	idx_t SAMColumnIndices::*index;
	bool integer;
};

static const SAMTagColumn SAM_TAG_COLUMNS[] = {
    {"tag_as", "AS", &SAMColumnIndices::tag_as_idx, true},  {"tag_xs", "XS", &SAMColumnIndices::tag_xs_idx, true},
    {"tag_ys", "YS", &SAMColumnIndices::tag_ys_idx, true},  {"tag_xn", "XN", &SAMColumnIndices::tag_xn_idx, true},
    {"tag_xm", "XM", &SAMColumnIndices::tag_xm_idx, true},  {"tag_xo", "XO", &SAMColumnIndices::tag_xo_idx, true},
    {"tag_xg", "XG", &SAMColumnIndices::tag_xg_idx, true},  {"tag_nm", "NM", &SAMColumnIndices::tag_nm_idx, true},
    {"tag_yt", "YT", &SAMColumnIndices::tag_yt_idx, false}, {"tag_md", "MD", &SAMColumnIndices::tag_md_idx, false},
    {"tag_sa", "SA", &SAMColumnIndices::tag_sa_idx, false},
};

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
enum class SAMOutputFormat : uint8_t { SAM = 0, BAM = 1 };

// "SAM" or "BAM", as errors name the format
static const char *SAMFormatName(SAMOutputFormat format) {
	return format == SAMOutputFormat::BAM ? "BAM" : "SAM";
}

struct SAMCopyBindData : public FunctionData {
	bool include_header = true;
	bool use_gzip = false;
	SAMOutputFormat format = SAMOutputFormat::SAM;
	int compression_level = -1; // -1 means use HTSlib default (6), for BAM and gzip SAM
	string file_path;
	vector<string> names;
	std::optional<string> reference_lengths_table;
	std::optional<string> sequence_data_table;
	SAMColumnIndices indices; // Cached column indices

	// Captured types of the three identifier columns. Each may be VARCHAR or
	// BIGINT (validated at bind). Default-constructed LogicalType is INVALID so
	// any path that forgets to set them fails loud at the record encoder's dispatch.
	LogicalType read_id_type;
	LogicalType reference_type;
	LogicalType mate_reference_type;

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<SAMCopyBindData>();
		result->include_header = include_header;
		result->use_gzip = use_gzip;
		result->format = format;
		result->compression_level = compression_level;
		result->file_path = file_path;
		result->names = names;
		result->reference_lengths_table = reference_lengths_table;
		result->sequence_data_table = sequence_data_table;
		result->indices = indices;
		result->read_id_type = read_id_type;
		result->reference_type = reference_type;
		result->mate_reference_type = mate_reference_type;
		return std::move(result);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<SAMCopyBindData>();
		return include_header == other.include_header && use_gzip == other.use_gzip && format == other.format &&
		       compression_level == other.compression_level && file_path == other.file_path &&
		       reference_lengths_table == other.reference_lengths_table &&
		       sequence_data_table == other.sequence_data_table;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> SAMCopyBindInternal(ClientContext &context, CopyFunctionBindInput &input,
                                                    const vector<Identifier> &names,
                                                    const vector<LogicalType> &sql_types,
                                                    SAMOutputFormat default_format) {
	auto result = make_uniq<SAMCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = IdentifiersToStrings(names);
	result->format = default_format;
	const char *format_name = SAMFormatName(default_format);
	RefuseRemoteCopyTarget(format_name, result->file_path);

	// Detect and cache column indices
	result->indices.FindIndices(IdentifiersToStrings(names));
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
	if (!IsAllowedIdType(sql_types[indices.read_id_idx])) {
		throw BinderException("Column 'read_id' must be %s", AllowedIdTypeList());
	}
	if (sql_types[indices.flags_idx].id() != LogicalTypeId::USMALLINT) {
		throw BinderException("Column 'flags' must be USMALLINT");
	}
	if (!IsAllowedIdType(sql_types[indices.reference_idx])) {
		throw BinderException("Column 'reference' must be %s", AllowedIdTypeList());
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
	if (!IsAllowedIdType(sql_types[indices.mate_reference_idx])) {
		throw BinderException("Column 'mate_reference' must be %s", AllowedIdTypeList());
	}

	// Capture id-column types for the record encoder's per-row dispatch.
	result->read_id_type = sql_types[indices.read_id_idx];
	result->reference_type = sql_types[indices.reference_idx];
	result->mate_reference_type = sql_types[indices.mate_reference_idx];
	if (sql_types[indices.mate_position_idx].id() != LogicalTypeId::BIGINT) {
		throw BinderException("Column 'mate_position' must be BIGINT");
	}
	if (sql_types[indices.template_length_idx].id() != LogicalTypeId::BIGINT) {
		throw BinderException("Column 'template_length' must be BIGINT");
	}
	// The tag columns are read with typed accessors too: a column of another type would raise an internal error, which
	// invalidates the database.
	for (auto &tag : SAM_TAG_COLUMNS) {
		auto index = indices.*tag.index;
		auto expected = tag.integer ? LogicalTypeId::BIGINT : LogicalTypeId::VARCHAR;
		if (index != DConstants::INVALID_INDEX && sql_types[index].id() != expected) {
			throw BinderException("Column '%s' must be %s", tag.column, tag.integer ? "BIGINT" : "VARCHAR");
		}
	}

	// Parse options
	Value compression_param;
	const vector<Value> *compression_level_values = nullptr;
	for (auto &option : input.info.options) {
		if (option.first == "include_header") {
			result->include_header = CopyOptionFlag(option.second);
		} else if (option.first == "compression") {
			compression_param = option.second[0];
		} else if (option.first == "compression_level") {
			compression_level_values = &option.second;
		} else if (option.first == "reference_lengths") {
			result->reference_lengths_table = option.second[0].ToString();

			// Validate table or view exists (use TABLE_ENTRY lookup which returns either)
			EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY,
			                            QualifiedName(Identifier(result->reference_lengths_table.value())),
			                            QueryErrorContext());
			auto entry = Catalog::GetEntry(context, lookup_info, OnEntryNotFound::RETURN_NULL);
			if (!entry) {
				throw BinderException("Table or view '%s' does not exist", result->reference_lengths_table.value());
			}
			if (entry->type != CatalogType::TABLE_ENTRY && entry->type != CatalogType::VIEW_ENTRY) {
				throw BinderException("'%s' is not a table or view", result->reference_lengths_table.value());
			}
		} else if (option.first == "sequence_data") {
			result->sequence_data_table = option.second[0].ToString();

			// Validate table or view exists
			EntryLookupInfo sd_lookup(CatalogType::TABLE_ENTRY,
			                          QualifiedName(Identifier(result->sequence_data_table.value())),
			                          QueryErrorContext());
			auto sd_entry = Catalog::GetEntry(context, sd_lookup, OnEntryNotFound::RETURN_NULL);
			if (!sd_entry) {
				throw BinderException("Table or view '%s' does not exist", result->sequence_data_table.value());
			}
			if (sd_entry->type != CatalogType::TABLE_ENTRY && sd_entry->type != CatalogType::VIEW_ENTRY) {
				throw BinderException("'%s' is not a table or view", result->sequence_data_table.value());
			}
		} else {
			ThrowUnreadCopyOption(format_name, option.first);
		}
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

	// BAM is always BGZF, so it takes any level; SAM takes one only for gzip output, as FASTQ does.
	if (result->format == SAMOutputFormat::SAM) {
		const auto compression = DetectCompressionType(result->file_path, compression_param);
		result->use_gzip = compression == FileCompressionType::GZIP;
		AddGzipFileExtension(input, compression);
		if (compression_level_values) {
			result->compression_level = GzipCompressionLevel(*compression_level_values, compression, result->file_path);
		}
	} else {
		// COMPRESSION 'gzip' or 'gz' only says what BAM always is; any other value is refused, not ignored.
		if (!compression_param.IsNull() && !IsGzipCompressionName(compression_param.ToString())) {
			throw BinderException("COPY FORMAT BAM is always BGZF-compressed, so COMPRESSION '%s' does not apply; "
			                      "COMPRESSION_LEVEL 0 writes its BGZF blocks uncompressed",
			                      compression_param.ToString());
		}
		if (compression_level_values) {
			result->compression_level = CopyOptionCompressionLevel(*compression_level_values);
		}
	}

	return std::move(result);
}

// Bind function for SAM format
static unique_ptr<FunctionData> SAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                            const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	return SAMCopyBindInternal(context, input, names, sql_types, SAMOutputFormat::SAM);
}

// Bind function for BAM format
static unique_ptr<FunctionData> BAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                            const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	return SAMCopyBindInternal(context, input, names, sql_types, SAMOutputFormat::BAM);
}

//===--------------------------------------------------------------------===//
// REFERENCE_LENGTHS and SEQUENCE_DATA, read once per COPY
//===--------------------------------------------------------------------===//
// DuckDB initializes each output file of a COPY on its own and gives a format no execution-time hook
// for the COPY as a whole. So the first file's InitializeGlobal reads the relations and the COPY's
// other files share what it read: every file gets the same header and the same sequences. What was
// read is kept in the client context until the query ends, so a prepared COPY reads the relations
// again each time it is executed.
struct SAMCopyRelations {
	std::vector<std::pair<std::string, uint64_t>> reference_lengths;
	// Each reference's TID, its index in reference_lengths: the order of the @SQ lines of headered output
	std::unordered_map<std::string, int32_t> reference_tids;
	std::optional<SequenceDataMap> sequence_data;
};

class SAMCopyRelationsCache : public ClientContextState {
public:
	shared_ptr<const SAMCopyRelations> GetOrRead(ClientContext &context, const SAMCopyBindData &fdata) {
		const CopyExecution copy {&fdata, context.transaction.GetActiveQuery()};
		// Held while reading, so the COPY's other files wait for this read instead of making their own
		lock_guard<mutex> guard(lock);
		// A context runs one statement at a time, so what another statement read is no longer needed
		for (auto it = by_copy.begin(); it != by_copy.end();) {
			it = it->first.second == copy.second ? std::next(it) : by_copy.erase(it);
		}
		auto &relations = by_copy[copy];
		if (!relations) {
			auto read = make_shared_ptr<SAMCopyRelations>();
			if (fdata.reference_lengths_table.has_value()) {
				// SORTED BY NAME, not the map's iteration order. A reference's index in the @SQ list IS
				// its TID, and a BAM is coordinate-sorted by TID rather than by reference name -- so
				// emitting @SQ in unordered_map hash order meant no ORDER BY the caller could write
				// produced a coordinate-sorted file, and htslib consumers rejected the output.
				// Name order is deterministic regardless of scan parallelism and makes the natural
				// `ORDER BY reference, position` correct.
				read->reference_lengths =
				    ReadReferenceTableSortedByName(context, fdata.reference_lengths_table.value());
				for (idx_t tid = 0; tid < read->reference_lengths.size(); tid++) {
					read->reference_tids.emplace(read->reference_lengths[tid].first, NumericCast<int32_t>(tid));
				}
			}
			if (fdata.sequence_data_table.has_value()) {
				read->sequence_data = ReadSequenceDataTable(context, fdata.sequence_data_table.value());
			}
			relations = std::move(read);
		}
		return relations;
	}

	void QueryEnd(ClientContext &context) override {
		lock_guard<mutex> guard(lock);
		by_copy.clear();
	}

private:
	// The COPY's bind data, one object for all of its files, and the statement running it. A query can
	// run several statements (IMPORT DATABASE runs its load.sql's), freeing each one's bind data
	// before the next is planned, so the address alone can repeat within a query.
	using CopyExecution = std::pair<const SAMCopyBindData *, transaction_t>;

	mutex lock;
	std::map<CopyExecution, shared_ptr<const SAMCopyRelations>> by_copy;
};

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
// One output file. What prepare_batch reads here, the header of headered output and the relations, is the same in every
// file of the COPY and is not changed after init (architecture.md, "Record writers on DuckDB's batch API",
// Invariants).
struct SAMCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	string file_path;
	SAMFilePtr sam_file;
	// Headered output: the @SQ lines of REFERENCE_LENGTHS, written at init and never changed after, so prepare_batch
	// formats with it on any thread (sam_format1 only reads it). Null for headerless SAM.
	SAMHeaderPtr header;
	// REFERENCE_LENGTHS and SEQUENCE_DATA as the COPY read them; null without either
	shared_ptr<const SAMCopyRelations> relations;
};

static unique_ptr<GlobalFunctionData> SAMCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                              const string &file_path) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();

	auto gstate = make_uniq<SAMCopyGlobalState>();
	gstate->file_path = file_path;

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
	} else if (fdata.use_gzip) {
		mode = fdata.compression_level >= 0 ? "wz" + std::to_string(fdata.compression_level) : "wz";
	} else {
		mode = "w";
	}

	// Open file for writing
	gstate->sam_file = SAMFilePtr(sam_open(file_path.c_str(), mode.c_str()));
	if (!gstate->sam_file) {
		throw IOException("Failed to open file for writing: " + file_path);
	}

	// BAM and gzip SAM compress their BGZF blocks on a worker pool, which does not change the bytes. bgzf_mt, because
	// hts_set_threads on SAM would also start a pool formatting records this writer has formatted (embedded-tools.md,
	// htslib).
	bool compressed_output = fdata.format == SAMOutputFormat::BAM || fdata.use_gzip;
	if (compressed_output) {
		auto n_threads = NumericCast<int>(TaskScheduler::GetScheduler(context).NumberOfThreads());
		if (n_threads > 1) {
			bgzf_mt(gstate->sam_file->fp.bgzf, n_threads, 256);
		}
	}

	if (fdata.reference_lengths_table.has_value() || fdata.sequence_data_table.has_value()) {
		gstate->relations = context.registered_state->GetOrCreate<SAMCopyRelationsCache>("miint_copy_sam_relations")
		                        ->GetOrRead(context, fdata);
	}

	// Headered output (bind requires REFERENCE_LENGTHS for it) writes its header now, so that every batch formats with
	// the header the file holds, and an empty COPY still writes one.
	if (fdata.include_header) {
		D_ASSERT(gstate->relations);
		gstate->header = SAMHeaderPtr(sam_hdr_init());
		if (!gstate->header) {
			throw IOException("Failed to create SAM header");
		}
		for (const auto &ref : gstate->relations->reference_lengths) {
			if (sam_hdr_add_line(gstate->header.get(), "SQ", "SN", ref.first.c_str(), "LN",
			                     std::to_string(ref.second).c_str(), NULL) < 0) {
				throw IOException("Failed to add reference to SAM header: %s", ref.first);
			}
		}
		errno = 0;
		if (sam_hdr_write(gstate->sam_file.get(), gstate->header.get()) < 0) {
			ThrowCopyWriteFailed(SAMFormatName(fdata.format), file_path, "sam_hdr_write");
		}
	}

	return std::move(gstate);
}

//===--------------------------------------------------------------------===//
// Helper Functions
//===--------------------------------------------------------------------===//

// Parse CIGAR string and store in buffer
static bool ParseCIGAR(const string &cigar_str, std::vector<uint32_t> &cigar_buffer) {
	cigar_buffer.clear();

	// "*" and "" are the two spellings of "no CIGAR": the SAM sentinel, and what
	// read_alignments/read_sam return for a record with zero CIGAR operations. Both parse to an
	// empty buffer. Whether a record is ALLOWED to have no CIGAR depends on its unmapped flag,
	// which the caller checks -- see the BAM_FUNMAP check in SAMRecordEncoder::Build.
	if (cigar_str == "*" || cigar_str.empty()) {
		return true;
	}

	// Use HTSlib's CIGAR parser
	char *end = nullptr;
	uint32_t *cigar_array = nullptr;
	size_t cigar_mem = 0;

	ssize_t n_cigar = sam_parse_cigar(cigar_str.c_str(), &end, &cigar_array, &cigar_mem);

	// A 0 return needs no separate test: both strings htslib returns 0 for -- "*" and one with
	// no operation letters, e.g. "123" -- either returned above or leave *end at the unconsumed
	// first byte, so the end check below already rejects them.
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

//===--------------------------------------------------------------------===//
// Sequence/Quality Preparation
//===--------------------------------------------------------------------===//

// Prepare SEQ and QUAL for a single alignment record by looking up the original
// FASTQ sequence, applying hard clipping and reverse complement as needed.
static void PrepareSeqQual(const SequenceDataMap &seq_data, const std::string &read_id, uint16_t flags,
                           const std::vector<uint32_t> &cigar, std::string &out_seq, std::vector<uint8_t> &out_qual) {
	auto it = seq_data.find(read_id);
	if (it == seq_data.end()) {
		throw InvalidInputException("Read '%s' not found in SEQUENCE_DATA table", read_id);
	}

	const auto &entry = it->second;

	// Pick read1 vs read2 based on paired-end flag
	bool is_read2 = (flags & 0x80) != 0;
	const std::string &raw_seq = is_read2 ? entry.sequence2 : entry.sequence1;
	const std::vector<uint8_t> &raw_qual = is_read2 ? entry.qual2 : entry.qual1;

	if (raw_seq.empty()) {
		throw InvalidInputException("Read '%s' has no %s sequence in SEQUENCE_DATA table", read_id,
		                            is_read2 ? "sequence2" : "sequence1");
	}

	// Validate sequence and quality lengths match (when quality is present)
	if (!raw_qual.empty() && raw_qual.size() != raw_seq.size()) {
		throw InvalidInputException("Quality score length (%llu) does not match sequence length (%llu) for read_id "
		                            "'%s'%s in the SEQUENCE_DATA table",
		                            static_cast<uint64_t>(raw_qual.size()), static_cast<uint64_t>(raw_seq.size()),
		                            read_id, is_read2 ? " (R2)" : "");
	}

	// Determine hard clipping from CIGAR
	int64_t leading_clip = 0;
	int64_t trailing_clip = 0;
	if (!cigar.empty()) {
		if (bam_cigar_op(cigar.front()) == BAM_CHARD_CLIP) {
			leading_clip = bam_cigar_oplen(cigar.front());
		}
		if (cigar.size() > 1 && bam_cigar_op(cigar.back()) == BAM_CHARD_CLIP) {
			trailing_clip = bam_cigar_oplen(cigar.back());
		}
	}

	int64_t full_len = static_cast<int64_t>(raw_seq.size());
	int64_t clipped_len = full_len - leading_clip - trailing_clip;

	if (clipped_len <= 0) {
		throw InvalidInputException("Hard clipping (%lld + %lld) exceeds sequence length (%lld) for read '%s'",
		                            leading_clip, trailing_clip, full_len, read_id);
	}

	out_seq = raw_seq.substr(leading_clip, clipped_len);
	if (!raw_qual.empty()) {
		out_qual.assign(raw_qual.begin() + leading_clip, raw_qual.begin() + leading_clip + clipped_len);
	} else {
		out_qual.clear();
	}

	// Reverse complement if on reverse strand
	if (flags & 0x10) {
		out_seq = miint::dna_reverse_complement(out_seq);
		std::reverse(out_qual.begin(), out_qual.end());
	}

	// Validate: SEQ length must match sum of query-consuming CIGAR ops (skip for unmapped)
	if (!cigar.empty()) {
		int64_t expected_len = 0;
		for (auto c : cigar) {
			int op = bam_cigar_op(c);
			// Query-consuming operations: M(0), I(1), S(4), =(7), X(8)
			if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
				expected_len += bam_cigar_oplen(c);
			}
		}
		if (static_cast<int64_t>(out_seq.size()) != expected_len) {
			throw InvalidInputException(
			    "After hard clipping, sequence length (%lld) does not match CIGAR query length (%lld) for read '%s'",
			    static_cast<int64_t>(out_seq.size()), expected_len, read_id);
		}
	}
}

//===--------------------------------------------------------------------===//
// Identifier-column dispatch (read_id, reference, mate_reference)
//===--------------------------------------------------------------------===//
// Returns the SAM textual representation of an identifier-column cell: the SAM '*' sentinel for a NULL, of every
// accepted type (reading an invalid slot gave an empty field: invalid SAM, and a BAM htslib refuses to read at all);
// otherwise the text the sequence writers give it (ResolveSequenceRecordId): BIGINT in decimal, UUID in its canonical
// 36-character lowercase form, VARCHAR verbatim.
static inline string ExtractSamIdField(const UnifiedVectorFormat &fmt, idx_t row_idx, const LogicalType &id_type) {
	if (!fmt.validity.RowIsValid(fmt.sel->get_index(row_idx))) {
		return "*";
	}
	string buf;
	const char *ptr;
	idx_t size;
	ResolveSequenceRecordId(fmt, row_idx, id_type, buf, ptr, size);
	return ptr == buf.data() ? std::move(buf) : string(ptr, size);
}

//===--------------------------------------------------------------------===//
// Required scalar columns
//===--------------------------------------------------------------------===//
// Reads a required fixed-width column, failing loud on NULL.
//
// Consulting validity is not optional: indexing the data array of an invalid slot is undefined,
// and in practice yields different garbage depending on the vector's representation. Silently
// accepting it is how a NULL flags could write a *mapped* record at an arbitrary position with the
// COPY reporting success -- the same defect class as #197's NULL identifiers.
//
// null_means_zero is set by the caller for the columns where SAM itself defines 0 as "absent" AND
// the record is unmapped, so nothing can be lost by substituting it.
template <class T>
static inline T RequiredScalar(const T *data, const UnifiedVectorFormat &fmt, idx_t row_idx, const char *column,
                               const string &read_id, bool null_means_zero = false) {
	auto idx = fmt.sel->get_index(row_idx);
	if (fmt.validity.RowIsValid(idx)) {
		return data[idx];
	}
	if (null_means_zero) {
		return static_cast<T>(0);
	}
	throw InvalidInputException("Column '%s' is NULL for read: %s", column, read_id);
}

//===--------------------------------------------------------------------===//
// Record Encoder
//===--------------------------------------------------------------------===//

// The references one batch of headerless SAM meets, in the order it meets them, as a header for sam_format1 to name
// them from. It is never written. Records name their reference and mate reference, and write '=' where the two TIDs are
// equal, so the bytes do not depend on which TIDs a batch gives. Only the target arrays are filled, as bam_hdr_read
// fills them: sam_format1 reads nothing else. It is built anew for each batch, which meets up to two references a row,
// so a reference costs no allocation of its own: its name goes to an arena, and is found again through an
// open-addressing table of TIDs; the arrays and the table start with room for a chunk's rows and double.
class BatchReferences {
public:
	BatchReferences() : names(Allocator::DefaultAllocator()), slots(INITIAL_SLOTS, EMPTY) {
		target_names.reserve(INITIAL_SLOTS / 2);
		target_lengths.reserve(INITIAL_SLOTS / 2);
		entries.reserve(INITIAL_SLOTS / 2);
	}
	BatchReferences(const BatchReferences &) = delete;
	BatchReferences &operator=(const BatchReferences &) = delete;

	int32_t TID(const string &name) {
		const hash_t hash = Hash(name.c_str(), name.size());
		idx_t slot = hash & (slots.size() - 1);
		for (; slots[slot] != EMPTY; slot = (slot + 1) & (slots.size() - 1)) {
			const auto tid = slots[slot];
			if (entries[tid].hash == hash && entries[tid].size == name.size() &&
			    memcmp(target_names[tid], name.c_str(), name.size()) == 0) {
				return tid;
			}
		}
		return Add(name, hash, slot);
	}

	const sam_hdr_t &Header() const {
		return header;
	}

private:
	// A full chunk's rows, one reference each, at half load
	static constexpr idx_t INITIAL_SLOTS = 2 * STANDARD_VECTOR_SIZE;
	static constexpr int32_t EMPTY = -1;

	struct Entry {
		hash_t hash;
		idx_t size;
	};

	// A name the batch has not met, at the empty `slot` its probe ended on. Kept out of TID, which runs twice a row
	// (reference and mate reference), so that TID stays small where it is inlined.
	[[gnu::noinline]] int32_t Add(const string &name, hash_t hash, idx_t slot) {
		const auto tid = NumericCast<int32_t>(target_names.size());
		auto copy = char_ptr_cast(names.Allocate(name.size() + 1));
		memcpy(copy, name.c_str(), name.size() + 1);
		target_names.push_back(copy);
		target_lengths.push_back(0); // unknown, and never written
		entries.push_back(Entry {hash, name.size()});
		slots[slot] = tid;
		if (2 * target_names.size() > slots.size()) {
			Grow();
		}
		// The vectors may have moved
		header.n_targets = tid + 1;
		header.target_name = target_names.data();
		header.target_len = target_lengths.data();
		return tid;
	}

	// Doubles the table, which is never more than half full
	void Grow() {
		vector<int32_t> grown(2 * slots.size(), EMPTY);
		for (idx_t tid = 0; tid < entries.size(); tid++) {
			idx_t slot = entries[tid].hash & (grown.size() - 1);
			while (grown[slot] != EMPTY) {
				slot = (slot + 1) & (grown.size() - 1);
			}
			grown[slot] = NumericCast<int32_t>(tid);
		}
		slots = std::move(grown);
	}

	ArenaAllocator names; // each name and the NUL sam_format1 reads it to
	vector<char *> target_names;
	vector<uint32_t> target_lengths;
	vector<Entry> entries; // each TID's name: its hash and length
	vector<int32_t> slots; // a TID or EMPTY; a power of two in size
	// Its arrays are the vectors above and its names are in the arena, so it is never passed to sam_hdr_destroy
	sam_hdr_t header {};
};

// The ENCODER RecordPrepareBatch and RecordFlushBatch take (copy_format_common.hpp). It builds a row into a bam1_t
// (every check the row must pass, SEQ and QUAL from SEQUENCE_DATA, the TIDs and the tags), then stores it as a record:
// for SAM, its sam_format1 line and '\n'; for BAM, as StoreBAMRecord stores it, which WriteSAMRecord writes in
// place. Made for one prepare_batch or flush_batch call, so one thread uses it. It reads only the global
// state's header and relations (see SAMCopyGlobalState).
class SAMRecordEncoder {
public:
	SAMRecordEncoder(const SAMCopyBindData &fdata_p, const SAMCopyGlobalState &gstate)
	    : fdata(fdata_p), relations(gstate.relations.get()), shared_header(gstate.header.get()), record(bam_init1()) {
		if (!record) {
			throw IOException("Failed to allocate BAM record");
		}
		if (!shared_header) {
			references = make_uniq<BatchReferences>();
		}
		for (auto &column : SAM_TAG_COLUMNS) {
			auto index = fdata.indices.*column.index;
			if (index != DConstants::INVALID_INDEX) {
				tags.push_back(Tag {&column, index, {}});
			}
		}
	}

	SAMRecordEncoder(const SAMRecordEncoder &) = delete;
	SAMRecordEncoder &operator=(const SAMRecordEncoder &) = delete;
	~SAMRecordEncoder() {
		ks_free(&line);
	}

	void SetChunk(DataChunk &input) {
		auto &indices = fdata.indices;
		input.data[indices.read_id_idx].ToUnifiedFormat(read_id_data);
		input.data[indices.flags_idx].ToUnifiedFormat(flags_data);
		input.data[indices.reference_idx].ToUnifiedFormat(reference_data);
		input.data[indices.position_idx].ToUnifiedFormat(position_data);
		input.data[indices.mapq_idx].ToUnifiedFormat(mapq_data);
		input.data[indices.cigar_idx].ToUnifiedFormat(cigar_data);
		input.data[indices.mate_reference_idx].ToUnifiedFormat(mate_reference_data);
		input.data[indices.mate_position_idx].ToUnifiedFormat(mate_position_data);
		input.data[indices.template_length_idx].ToUnifiedFormat(template_length_data);
		for (auto &tag : tags) {
			input.data[tag.index].ToUnifiedFormat(tag.data);
		}
	}

	// Checks row `row` of the chunk and appends its record to `out`
	void Encode(idx_t row, RecordBuffer &out);

	// Checks row `row` of the chunk, raising every data error Encode would
	void Check(idx_t row) {
		Build(row);
	}

private:
	// Checks row `i` of the chunk and builds it into `record`. Every data error a row can have is raised here.
	void Build(idx_t i);
	int32_t ReferenceTID(const string &reference, const char *column);
	void AppendTags(idx_t i);
	void CheckBAMLimits() const;

	struct Tag {
		const SAMTagColumn *column;
		idx_t index;
		UnifiedVectorFormat data;
	};

	const SAMCopyBindData &fdata;
	const SAMCopyRelations *relations;
	const sam_hdr_t *shared_header;
	unique_ptr<BatchReferences> references; // headerless SAM only
	BAMRecordPtr record;
	std::vector<uint32_t> cigar_buffer;
	std::string seq_buffer;
	std::vector<uint8_t> qual_buffer;
	string read_id;
	UnifiedVectorFormat read_id_data, flags_data, reference_data, position_data, mapq_data;
	UnifiedVectorFormat cigar_data, mate_reference_data, mate_position_data, template_length_data;
	vector<Tag> tags;               // the tag columns the input has, in SAM_TAG_COLUMNS' order
	kstring_t line = KS_INITIALIZE; // SAM: the record's sam_format1 line
};

// A reference's TID: its index in REFERENCE_LENGTHS (sorted by name) for headered output, or this batch's own for
// headerless SAM, which names references only by name.
//
// Headered output has written its @SQ lines, so a reference not among them cannot be written: SAM read it back as
// '*', BAM dropped the record outright, neither with any error. Supplying REFERENCE_LENGTHS asserts the table
// is the complete reference set, so a gap is a caller error.
int32_t SAMRecordEncoder::ReferenceTID(const string &reference, const char *column) {
	if (reference == "*") {
		return -1;
	}
	if (references) {
		return references->TID(reference);
	}
	auto entry = relations->reference_tids.find(reference);
	if (entry == relations->reference_tids.end()) {
		throw InvalidInputException(
		    "Column '%s' is '%s' for read '%s', which is not in the REFERENCE_LENGTHS table. The header has "
		    "already been written and cannot gain new @SQ lines, so this record could not be represented. "
		    "REFERENCE_LENGTHS must cover every reference the records use (or use INCLUDE_HEADER=false to write "
		    "headerless SAM).",
		    column, reference, read_id);
	}
	return entry->second;
}

// Appends the row's tags, in SAM_TAG_COLUMNS' order; a NULL is no tag. An integer tag is a 32-bit signed value in SAM
// and BAM alike, so one outside that range is an error rather than a value cut to 32 bits.
void SAMRecordEncoder::AppendTags(idx_t i) {
	for (auto &tag : tags) {
		auto data_idx = tag.data.sel->get_index(i);
		if (!tag.data.validity.RowIsValid(data_idx)) {
			continue;
		}
		int appended;
		if (tag.column->integer) {
			auto value = UnifiedVectorFormat::GetData<int64_t>(tag.data)[data_idx];
			if (value < NumericLimits<int32_t>::Minimum() || value > NumericLimits<int32_t>::Maximum()) {
				throw InvalidInputException(
				    "Column '%s' is %lld for read '%s', outside the 32-bit signed range of a SAM integer tag",
				    tag.column->column, static_cast<long long>(value), read_id);
			}
			auto value32 = static_cast<int32_t>(value);
			appended = bam_aux_append(record.get(), tag.column->tag, 'i', sizeof(int32_t),
			                          reinterpret_cast<uint8_t *>(&value32));
		} else {
			string value = UnifiedVectorFormat::GetData<string_t>(tag.data)[data_idx].GetString();
			appended = bam_aux_append(record.get(), tag.column->tag, 'Z', (int)value.size() + 1,
			                          reinterpret_cast<const uint8_t *>(value.c_str()));
		}
		if (appended < 0) {
			throw IOException("Failed to append %s tag to BAM record for read: %s", tag.column->tag, read_id);
		}
	}
}

// What bam_write1 refuses beyond what bam_set1 does, checked here so that a batch's data errors are raised while it is
// prepared and flush_batch writes every record it is given: positions and TLEN past 32 bits, and a CIGAR of more than
// 65535 operations that spans 2^28 reference bases or more.
void SAMRecordEncoder::CheckBAMLimits() const {
	const auto &c = record->core;
	// Positions are held 0-based; the columns, and these errors, are 1-based.
	const auto check_position = [&](const char *column, hts_pos_t pos) {
		if (pos > NumericLimits<int32_t>::Maximum()) {
			throw InvalidInputException(
			    "Column '%s' is %lld for read '%s', past %lld, the largest a BAM record holds (SAM text holds it)",
			    column, static_cast<long long>(pos + 1), read_id,
			    static_cast<long long>(NumericLimits<int32_t>::Maximum()) + 1);
		}
	};
	check_position("position", c.pos);
	check_position("mate_position", c.mpos);
	if (c.isize < NumericLimits<int32_t>::Minimum() || c.isize > NumericLimits<int32_t>::Maximum()) {
		throw InvalidInputException("Column 'template_length' is %lld for read '%s', outside the 32-bit signed range a "
		                            "BAM record holds (SAM text holds it)",
		                            static_cast<long long>(c.isize), read_id);
	}
	constexpr hts_pos_t LONG_CIGAR_SPAN = 1 << 28;
	if (c.n_cigar > 0xffff) {
		const auto span = bam_cigar2rlen(c.n_cigar, bam_get_cigar(record.get()));
		if (span >= LONG_CIGAR_SPAN) {
			throw InvalidInputException("Column 'cigar' of read '%s' has %lld operations spanning %lld reference "
			                            "bases: a BAM record holds more than 65535 only if they span fewer than %lld "
			                            "(SAM text holds it)",
			                            read_id, static_cast<long long>(c.n_cigar), static_cast<long long>(span),
			                            static_cast<long long>(LONG_CIGAR_SPAN));
		}
	}
}

void SAMRecordEncoder::Build(idx_t i) {
	auto flags_ptr = UnifiedVectorFormat::GetData<uint16_t>(flags_data);
	auto position_ptr = UnifiedVectorFormat::GetData<int64_t>(position_data);
	auto mapq_ptr = UnifiedVectorFormat::GetData<uint8_t>(mapq_data);
	auto cigar_ptr = UnifiedVectorFormat::GetData<string_t>(cigar_data);
	auto mate_position_ptr = UnifiedVectorFormat::GetData<int64_t>(mate_position_data);
	auto template_length_ptr = UnifiedVectorFormat::GetData<int64_t>(template_length_data);

	auto cigar_idx = cigar_data.sel->get_index(i);

	read_id = ExtractSamIdField(read_id_data, i, fdata.read_id_type);
	CheckBAMReadName(read_id.data(), read_id.size());

	// flags must be present: it is what tells us whether the record even claims to be aligned,
	// so every leniency below depends on it. A NULL here has no defensible reading.
	uint16_t flags = RequiredScalar(flags_ptr, flags_data, i, "flags", read_id);
	const bool row_unmapped = (flags & BAM_FUNMAP) != 0;

	string reference = ExtractSamIdField(reference_data, i, fdata.reference_type);

	// position and mapq are meaningless on an unmapped record, and SAM's own encoding for
	// "absent" is 0 -- unambiguous, with no alignment that could be silently discarded. So a
	// NULL is tolerated there and rejected on a mapped record. This is what lets
	// align_minimap2's include_unmapped output be written straight out; those rows carry NULL
	// position/mapq/cigar by design.
	int64_t position = RequiredScalar(position_ptr, position_data, i, "position", read_id, row_unmapped);
	uint8_t mapq = RequiredScalar(mapq_ptr, mapq_data, i, "mapq", read_id, row_unmapped);

	// Same rule for cigar, with '*' as the sentinel instead of 0. Note this must consult
	// validity rather than reading the slot: GetString() on an invalid slot is undefined and
	// observably varies -- "" from a materialized table, garbage from a constant-NULL vector.
	string cigar_str;
	if (cigar_data.validity.RowIsValid(cigar_idx)) {
		cigar_str = cigar_ptr[cigar_idx].GetString();
	} else if (row_unmapped) {
		cigar_str = "*";
	} else {
		throw InvalidInputException("Column 'cigar' is NULL for read: %s (only an unmapped record, FLAG 0x4, may "
		                            "omit it -- otherwise supply a CIGAR, or '*' plus the 0x4 flag)",
		                            read_id);
	}

	string mate_reference = ExtractSamIdField(mate_reference_data, i, fdata.mate_reference_type);
	int64_t mate_position = RequiredScalar(mate_position_ptr, mate_position_data, i, "mate_position", read_id);
	int64_t template_length = RequiredScalar(template_length_ptr, template_length_data, i, "template_length", read_id);

	// Parse CIGAR
	if (!ParseCIGAR(cigar_str, cigar_buffer)) {
		throw InvalidInputException("Failed to parse CIGAR string '%s' for read: %s", cigar_str, read_id);
	}

	// A record that claims to be mapped must carry a CIGAR. samtools accepts this on read but
	// silently rewrites the record as unmapped ("mapped query must have a CIGAR; treated as
	// unmapped"), so writing one degrades the alignment with nothing in the SQL layer to show
	// for it. Rewriting the caller's flags would be worse than refusing, so fail loud.
	// Unmapped records legitimately have no CIGAR -- that is the read_alignments round trip.
	if (cigar_buffer.empty() && (flags & BAM_FUNMAP) == 0) {
		throw InvalidInputException("Mapped read '%s' has no CIGAR (flags=%d): set the unmapped flag (0x4) or "
		                            "supply a CIGAR",
		                            read_id, flags);
	}

	// Validate position values before casting
	if (position < 0) {
		throw InvalidInputException("Invalid position value %lld for read: %s (must be non-negative)", position,
		                            read_id);
	}
	if (position > 0 && (position - 1) > HTS_POS_MAX) {
		throw InvalidInputException("Position value %lld exceeds maximum allowed position for read: %s", position,
		                            read_id);
	}
	if (mate_position < 0) {
		throw InvalidInputException("Invalid mate_position value %lld for read: %s (must be non-negative)",
		                            mate_position, read_id);
	}
	if (mate_position > 0 && (mate_position - 1) > HTS_POS_MAX) {
		throw InvalidInputException("Mate position value %lld exceeds maximum allowed position for read: %s",
		                            mate_position, read_id);
	}

	// Prepare SEQ and QUAL
	size_t l_seq = 0;
	const char *final_seq = nullptr;
	const char *final_qual = nullptr;

	if (relations && relations->sequence_data.has_value()) {
		PrepareSeqQual(relations->sequence_data.value(), read_id, flags, cigar_buffer, seq_buffer, qual_buffer);
		l_seq = seq_buffer.size();
		final_seq = seq_buffer.c_str();
		final_qual = qual_buffer.empty() ? nullptr : reinterpret_cast<const char *>(qual_buffer.data());
	}

	int32_t tid = ReferenceTID(reference, "reference");
	// "=" means the same as reference
	int32_t mtid = mate_reference == "=" ? tid : ReferenceTID(mate_reference, "mate_reference");

	// l_qname is the length WITHOUT null terminator (bam_set1 adds it internally)
	if (bam_set1(record.get(), read_id.length(), read_id.c_str(), flags, tid,
	             static_cast<hts_pos_t>(position - 1), // Convert to 0-based
	             mapq, cigar_buffer.size(), cigar_buffer.data(), mtid,
	             static_cast<hts_pos_t>(mate_position - 1), // Convert to 0-based
	             static_cast<hts_pos_t>(template_length), l_seq, final_seq, final_qual, 0) < 0) {
		throw IOException("Failed to build BAM record for read: " + read_id);
	}

	AppendTags(i);

	if (fdata.format == SAMOutputFormat::BAM) {
		CheckBAMLimits();
	}
}

void SAMRecordEncoder::Encode(idx_t row, RecordBuffer &out) {
	Build(row);
	if (fdata.format == SAMOutputFormat::BAM) {
		StoreBAMRecord(*record, out);
		return;
	}
	// Headered output's header, or this encoder's own references for headerless SAM
	const auto &header = shared_header ? *shared_header : references->Header();
	if (sam_format1(&header, record.get(), &line) < 0) {
		throw IOException("Failed to write SAM record for read: " + read_id);
	}
	auto data = out.Reserve(line.l + 1);
	memcpy(data, line.s, line.l);
	data[line.l] = '\n';
}

//===--------------------------------------------------------------------===//
// Prepare and Flush (DuckDB's batch API)
//===--------------------------------------------------------------------===//
static unique_ptr<PreparedBatchData> SAMCopyPrepareBatch(ClientContext &context, FunctionData &bind_data,
                                                         GlobalFunctionData &gstate_p,
                                                         unique_ptr<ColumnDataCollection> rows) {
	SAMRecordEncoder encoder(bind_data.Cast<SAMCopyBindData>(), gstate_p.Cast<SAMCopyGlobalState>());
	return RecordPrepareBatch(context, encoder, std::move(rows));
}

// Writes one record as the matching sam_write1 call would: BAM through WriteBAMRecord; SAM text as unthreaded
// sam_write1 writes a line, ending a gzip block before a line that does not fit. So the bytes are those of writing each
// record with sam_write1, whatever the batches. Caller holds gstate.lock.
static void WriteSAMRecord(const SAMCopyBindData &fdata, SAMCopyGlobalState &gstate, const_data_ptr_t data,
                           idx_t size) {
	auto fp = gstate.sam_file.get();
	if (fdata.format == SAMOutputFormat::BAM) {
		WriteBAMRecord(fp, gstate.header.get(), data, size, SAMFormatName(fdata.format), gstate.file_path);
		return;
	}
	errno = 0;
	if (fdata.use_gzip) {
		auto bgzf = fp->fp.bgzf;
		if (bgzf_flush_try(bgzf, NumericCast<ssize_t>(size)) < 0 ||
		    bgzf_write(bgzf, data, size) != NumericCast<ssize_t>(size)) {
			ThrowCopyWriteFailed(SAMFormatName(fdata.format), gstate.file_path, "bgzf_write");
		}
	} else if (hwrite(fp->fp.hfile, data, size) != NumericCast<ssize_t>(size)) {
		ThrowCopyWriteFailed(SAMFormatName(fdata.format), gstate.file_path, "hwrite");
	}
}

// DuckDB calls flush_batch for one batch at a time per file, in input order when order is kept.
static void SAMCopyFlushBatch(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                              PreparedBatchData &batch_p) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();
	auto &batch = batch_p.Cast<RecordBatchData>();
	lock_guard<mutex> glock(gstate.lock);
	RecordFlushBatch(
	    context, batch, [&]() { return SAMRecordEncoder(fdata, gstate); },
	    [&](const_data_ptr_t data, idx_t size) { WriteSAMRecord(fdata, gstate, data, size); });
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void SAMCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();
	lock_guard<mutex> glock(gstate.lock);
	CloseSAMCopyFile(gstate.sam_file.release(), SAMFormatName(fdata.format), gstate.file_path);
}

//===--------------------------------------------------------------------===//
// Function Definition
//===--------------------------------------------------------------------===//
// FORMAT SAM and BAM are two registrations each: the public one ("sam") declares the options, the extension and a plan
// callback, which rebinds the statement as the writer ("_miint_sam"). See SequenceCopyPlan.
static constexpr const char *SAM_WRITER_FORMAT = "_miint_sam";
static constexpr const char *BAM_WRITER_FORMAT = "_miint_bam";

static void SAMCopyOptions(ClientContext &context, CopyOptionsInput &input) {
	DeclareCopyOption(input, "include_header", LogicalType::BOOLEAN);
	DeclareCopyOption(input, "compression", LogicalType::VARCHAR);
	DeclareCompressionLevelCopyOption(input);
	DeclareCopyOption(input, "reference_lengths", LogicalType::VARCHAR);
	DeclareCopyOption(input, "sequence_data", LogicalType::VARCHAR);
}

static BoundStatement SAMCopyPlan(Binder &binder, CopyStatement &stmt) {
	return SequenceCopyPlan(binder, stmt, SAM_WRITER_FORMAT);
}

static BoundStatement BAMCopyPlan(Binder &binder, CopyStatement &stmt) {
	return SequenceCopyPlan(binder, stmt, BAM_WRITER_FORMAT);
}

CopyFunction CopySAMFunction::GetFunction() {
	CopyFunction function("sam");
	function.copy_options = SAMCopyOptions;
	function.plan = SAMCopyPlan;
	function.extension = "sam";
	return function;
}

CopyFunction CopySAMFunction::GetBAMFunction() {
	CopyFunction function("bam");
	function.copy_options = SAMCopyOptions;
	function.plan = BAMCopyPlan;
	function.extension = "bam";
	return function;
}

// FORMAT _miint_sam / _miint_bam: the writers. Named directly, they write without the default batch cap.
static CopyFunction GetWriterFunction(const char *name, copy_to_bind_t bind, const char *extension) {
	CopyFunction function(name);
	function.copy_options = SAMCopyOptions;
	function.copy_to_bind = bind;
	function.copy_to_initialize_global = SAMCopyInitializeGlobal;
	function.copy_to_finalize = SAMCopyFinalize;
	function.execution_mode = SequenceCopyExecutionMode;
	function.prepare_batch = SAMCopyPrepareBatch;
	function.flush_batch = SAMCopyFlushBatch;
	function.extension = extension;
	return function;
}

void CopySAMFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
	loader.RegisterFunction(GetBAMFunction());
	loader.RegisterFunction(GetWriterFunction(SAM_WRITER_FORMAT, SAMCopyBind, "sam"));
	loader.RegisterFunction(GetWriterFunction(BAM_WRITER_FORMAT, BAMCopyBind, "bam"));
}

} // namespace duckdb
