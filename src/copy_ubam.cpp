#include "copy_ubam.hpp"
#include "copy_format_common.hpp" // ResolveSequenceRecordId, SequenceCopyPlan, RecordPrepareBatch/RecordFlushBatch
#include "copy_option_parsing.hpp"
#include "htslib_raii.hpp"
#include "QualScore.hpp"       // MAX_PHRED_QUALITY
#include "id_column_utils.hpp" // IsAllowedIdType, AllowedIdTypeList
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/operator/cast_operators.hpp" // TryCast
#include "duckdb/common/types.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/planner/binder.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/list_vector.hpp"
#include <htslib-1.22.1/htslib/hts.h>
#include <htslib-1.22.1/htslib/sam.h>
#include <algorithm>
#include <cerrno>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct UBAMTagColumn {
	string tag;       // 2-char BAM aux tag name (e.g. "zm"), verbatim
	idx_t col_idx;    // position of the source column in the COPY projection
	LogicalType type; // integer type of the source column

	bool operator==(const UBAMTagColumn &o) const {
		return tag == o.tag && col_idx == o.col_idx && type == o.type;
	}
};

struct UBAMCopyBindData : public FunctionData {
	string file_path;
	vector<string> names;
	idx_t read_id_idx = DConstants::INVALID_INDEX;
	idx_t sequence1_idx = DConstants::INVALID_INDEX;
	idx_t qual1_idx = DConstants::INVALID_INDEX;
	LogicalType read_id_type;
	int compression_level = -1; // -1 => HTSlib default (6)

	// @RG (optional). rg_line is the full "@RG\t...\n" header text, built at bind;
	// rg_id is the ID value, emitted as a per-record RG:Z tag.
	bool has_read_group = false;
	string rg_line;
	string rg_id;

	// Per-read integer aux tags (optional).
	vector<UBAMTagColumn> tags;

	unique_ptr<FunctionData> Copy() const override {
		auto r = make_uniq<UBAMCopyBindData>();
		r->file_path = file_path;
		r->names = names;
		r->read_id_idx = read_id_idx;
		r->sequence1_idx = sequence1_idx;
		r->qual1_idx = qual1_idx;
		r->read_id_type = read_id_type;
		r->compression_level = compression_level;
		r->has_read_group = has_read_group;
		r->rg_line = rg_line;
		r->rg_id = rg_id;
		r->tags = tags;
		return std::move(r);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<UBAMCopyBindData>();
		return file_path == o.file_path && names == o.names && read_id_type == o.read_id_type &&
		       compression_level == o.compression_level && has_read_group == o.has_read_group && rg_line == o.rg_line &&
		       rg_id == o.rg_id && tags == o.tags;
	}
};

//===--------------------------------------------------------------------===//
// Bind helpers
//===--------------------------------------------------------------------===//
static bool IsIntegerColumnType(const LogicalType &t) {
	switch (t.id()) {
	case LogicalTypeId::TINYINT:
	case LogicalTypeId::SMALLINT:
	case LogicalTypeId::INTEGER:
	case LogicalTypeId::BIGINT:
	case LogicalTypeId::UTINYINT:
	case LogicalTypeId::USMALLINT:
	case LogicalTypeId::UINTEGER:
	case LogicalTypeId::UBIGINT:
		return true;
	default:
		return false;
	}
}

// Build the @RG header line (and capture the ID) from a READ_GROUP struct option.
// Field names are upper-cased (SAM @RG codes are uppercase: ID, PL, DS, SM, ...);
// ID is required and always emitted first per the SAM spec.
static void ParseReadGroup(const Value &rg, string &out_line, string &out_id) {
	if (rg.type().id() != LogicalTypeId::STRUCT) {
		throw BinderException(
		    "READ_GROUP must be a struct, e.g. READ_GROUP {ID: 'qiita', PL: 'PACBIO', DS: 'READTYPE=CCS'}");
	}
	auto &child_types = StructType::GetChildTypes(rg.type());
	auto children = StructValue::GetChildren(rg);

	string id;
	string rest; // non-ID fields, in user order
	for (idx_t i = 0; i < child_types.size(); i++) {
		string field = StringUtil::Upper(child_types[i].first.GetIdentifierName());
		// @RG fields are 2-character SAM codes (ID, PL, DS, SM, LB, ...). Reject
		// anything else at bind rather than emitting a malformed @RG line that
		// sam_hdr_add_lines fails on later with an opaque error.
		if (field.size() != 2) {
			throw BinderException("READ_GROUP field '%s' must be a 2-character SAM @RG code (e.g. ID, PL, DS, SM)",
			                      field);
		}
		if (children[i].IsNull()) {
			throw BinderException("READ_GROUP field '%s' must not be NULL", field);
		}
		string val = children[i].ToString();
		if (val.find('\t') != string::npos || val.find('\n') != string::npos) {
			throw BinderException("READ_GROUP field '%s' must not contain tab or newline", field);
		}
		if (field == "ID") {
			id = val;
		} else {
			rest += "\t" + field + ":" + val;
		}
	}
	if (id.empty()) {
		throw BinderException("READ_GROUP requires an 'ID' field");
	}
	out_id = id;
	out_line = "@RG\tID:" + id + rest + "\n";
}

// Resolve the TAGS struct ({tag: 'column'}) into (tag, col_idx, type) entries.
// Tag names are taken verbatim (2-char, case-sensitive -- unquoted struct keys are
// lower-cased by the parser, so uppercase tags need quoted keys). Values name an
// integer column in the projection.
static void ParseTags(const Value &tags, const vector<string> &names, vector<UBAMTagColumn> &out,
                      const vector<LogicalType> &sql_types) {
	if (tags.type().id() != LogicalTypeId::STRUCT) {
		throw BinderException("TAGS must be a struct, e.g. TAGS {zm: zmw}");
	}
	auto &child_types = StructType::GetChildTypes(tags.type());
	auto children = StructValue::GetChildren(tags);
	for (idx_t i = 0; i < child_types.size(); i++) {
		const string &tag = child_types[i].first.GetIdentifierName();
		if (tag.size() != 2) {
			throw BinderException("TAGS tag name '%s' must be exactly 2 characters", tag);
		}
		if (StringUtil::CIEquals(tag, "RG")) {
			throw BinderException("TAGS may not set the 'RG' tag; use READ_GROUP instead");
		}
		if (children[i].IsNull()) {
			throw BinderException("TAGS column name for tag '%s' must not be NULL", tag);
		}
		string col = children[i].ToString();

		idx_t col_idx = DConstants::INVALID_INDEX;
		for (idx_t j = 0; j < names.size(); j++) {
			if (names[j] == col) {
				col_idx = j;
				break;
			}
		}
		if (col_idx == DConstants::INVALID_INDEX) {
			throw BinderException("TAGS column '%s' (for tag '%s') not found in the input", col, tag);
		}
		if (!IsIntegerColumnType(sql_types[col_idx])) {
			throw BinderException("TAGS column '%s' (for tag '%s') must be an integer type of at most 64 bits "
			                      "(TINYINT, SMALLINT, INTEGER, BIGINT, UTINYINT, USMALLINT, UINTEGER or UBIGINT), "
			                      "not %s",
			                      col, tag, sql_types[col_idx].ToString());
		}
		out.push_back(UBAMTagColumn {tag, col_idx, sql_types[col_idx]});
	}
}

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> UBAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                             const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<UBAMCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = IdentifiersToStrings(names);
	RefuseRemoteCopyTarget("UBAM", result->file_path);

	// Locate required columns (case-insensitive, matching the SAM/BAM writer).
	for (idx_t i = 0; i < names.size(); i++) {
		if (names[i] == "read_id") {
			result->read_id_idx = i;
		} else if (names[i] == "sequence1") {
			result->sequence1_idx = i;
		} else if (names[i] == "qual1") {
			result->qual1_idx = i;
		}
	}
	if (result->read_id_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT UBAM requires 'read_id' column");
	}
	if (result->sequence1_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT UBAM requires 'sequence1' column");
	}
	if (result->qual1_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT UBAM requires 'qual1' column");
	}

	// Types.
	if (!IsAllowedIdType(sql_types[result->read_id_idx])) {
		throw BinderException("Column 'read_id' must be %s", AllowedIdTypeList());
	}
	result->read_id_type = sql_types[result->read_id_idx];
	if (sql_types[result->sequence1_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'sequence1' must be VARCHAR");
	}
	if (sql_types[result->qual1_idx].id() != LogicalTypeId::LIST ||
	    ListType::GetChildType(sql_types[result->qual1_idx]).id() != LogicalTypeId::UTINYINT) {
		throw BinderException("Column 'qual1' must be UTINYINT[]");
	}

	// Options.
	for (auto &option : input.info.options) {
		if (option.first == "read_group") {
			ParseReadGroup(CopyOptionValue(option.first, option.second), result->rg_line, result->rg_id);
			result->has_read_group = true;
		} else if (option.first == "tags") {
			ParseTags(CopyOptionValue(option.first, option.second), IdentifiersToStrings(names), result->tags,
			          sql_types);
		} else if (option.first == "compression_level") {
			result->compression_level = CopyOptionCompressionLevel(option.second);
		} else if (option.first == "reference_lengths") {
			throw BinderException("COPY FORMAT UBAM does not accept REFERENCE_LENGTHS: a uBAM's header has no @SQ "
			                      "lines, as its reads are unaligned. Use FORMAT BAM for aligned records.");
		} else {
			ThrowUnreadCopyOption("UBAM", option.first);
		}
	}

	return std::move(result);
}

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
// One output file. prepare_batch reads none of it: a record is built from bind data alone (architecture.md, "Record
// writers on DuckDB's batch API", Invariants).
struct UBAMCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	string file_path;
	SAMFilePtr sam_file;
	SAMHeaderPtr header;
};

static unique_ptr<GlobalFunctionData> UBAMCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                               const string &file_path) {
	auto &fdata = bind_data.Cast<UBAMCopyBindData>();
	auto gstate = make_uniq<UBAMCopyGlobalState>();
	gstate->file_path = file_path;

	string mode = fdata.compression_level >= 0 ? "wb" + std::to_string(fdata.compression_level) : "wb";
	gstate->sam_file = SAMFilePtr(sam_open(file_path.c_str(), mode.c_str()));
	if (!gstate->sam_file) {
		throw IOException("Failed to open file for writing: " + file_path);
	}

	// BGZF blocks are compressed on a worker pool, which does not change the bytes (embedded-tools.md, htslib).
	auto n_threads = NumericCast<int>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	if (n_threads > 1) {
		hts_set_threads(gstate->sam_file.get(), n_threads);
	}

	// @HD and the optional @RG, and no @SQ: the reads are unaligned.
	gstate->header = SAMHeaderPtr(sam_hdr_init());
	if (!gstate->header) {
		throw IOException("Failed to create BAM header");
	}
	if (sam_hdr_add_line(gstate->header.get(), "HD", "VN", "1.6", "SO", "unknown", NULL) < 0) {
		throw IOException("Failed to add @HD line to uBAM header");
	}
	if (fdata.has_read_group) {
		if (sam_hdr_add_lines(gstate->header.get(), fdata.rg_line.c_str(), fdata.rg_line.size()) < 0) {
			throw IOException("Failed to add @RG line to uBAM header");
		}
	}
	// Written now, so that an empty COPY still writes a valid (header-only) BAM.
	errno = 0;
	if (sam_hdr_write(gstate->sam_file.get(), gstate->header.get()) < 0) {
		ThrowCopyWriteFailed("UBAM", file_path, "sam_hdr_write");
	}

	return std::move(gstate);
}

//===--------------------------------------------------------------------===//
// Record Encoder
//===--------------------------------------------------------------------===//
// The ENCODER RecordPrepareBatch and RecordFlushBatch take (copy_format_common.hpp). It reads a row and raises its data
// errors, then builds it into a bam1_t and stores it (StoreBAMRecord), which the flush writes (WriteBAMRecord).
// Made for one prepare_batch or flush_batch call, so one thread uses it. It reads only bind data.
class UBAMRecordEncoder {
public:
	explicit UBAMRecordEncoder(const UBAMCopyBindData &fdata_p)
	    : fdata(fdata_p), record(bam_init1()), tag_data(fdata_p.tags.size()), tag_values(fdata_p.tags.size()),
	      tag_present(fdata_p.tags.size()) {
		if (!record) {
			throw IOException("Failed to allocate BAM record");
		}
	}

	void SetChunk(DataChunk &input) {
		input.data[fdata.read_id_idx].ToUnifiedFormat(read_id_data);
		input.data[fdata.sequence1_idx].ToUnifiedFormat(seq1_data);
		qual1.Init(input.data[fdata.qual1_idx]);
		for (idx_t t = 0; t < fdata.tags.size(); t++) {
			input.data[fdata.tags[t].col_idx].ToUnifiedFormat(tag_data[t]);
		}
	}

	// Checks row `row` of the chunk and appends its record to `out`
	void Encode(idx_t row, RecordBuffer &out);

	// Checks row `row` of the chunk, raising every data error Encode would
	void Check(idx_t row) {
		Read(row);
	}

private:
	// Reads row `row` of the chunk into the members below. Every data error a row can have is raised here, naming the
	// row's read_id as FORMAT FASTQ's errors do.
	void Read(idx_t row);
	int32_t TagValue(idx_t t, idx_t idx) const;
	template <class T>
	int32_t WideTagValue(idx_t t, idx_t idx) const;

	string ReadId() const {
		return string(id_ptr, id_size);
	}

	const UBAMCopyBindData &fdata;
	BAMRecordPtr record;
	UnifiedVectorFormat read_id_data, seq1_data;
	QualColumn qual1;
	vector<UnifiedVectorFormat> tag_data;

	// The row Read last read
	string id_buf; // a BIGINT/UUID read_id as text
	const char *id_ptr = nullptr;
	idx_t id_size = 0;
	const char *seq_ptr = nullptr;
	idx_t seq_size = 0;
	const uint8_t *qual_ptr = nullptr;
	vector<int32_t> tag_values;
	vector<bool> tag_present; // a NULL is no tag
};

// A cell of an integer type with values past int32's range.
template <class T>
int32_t UBAMRecordEncoder::WideTagValue(idx_t t, idx_t idx) const {
	auto v = UnifiedVectorFormat::GetData<T>(tag_data[t])[idx];
	int32_t result;
	if (TryCast::Operation(v, result, false)) {
		return result;
	}
	auto &tag = fdata.tags[t];
	throw InvalidInputException(
	    "Column '%s' (tag '%s') is %s for read '%s', outside the 32-bit signed range of a SAM integer tag",
	    fdata.names[tag.col_idx], tag.tag, std::to_string(v), ReadId());
}

// A tag cell as a SAM integer tag's 32-bit signed value; one outside that range is an error rather than a value cut to
// 32 bits.
int32_t UBAMRecordEncoder::TagValue(idx_t t, idx_t idx) const {
	auto &fmt = tag_data[t];
	switch (fdata.tags[t].type.id()) {
	case LogicalTypeId::TINYINT:
		return UnifiedVectorFormat::GetData<int8_t>(fmt)[idx];
	case LogicalTypeId::SMALLINT:
		return UnifiedVectorFormat::GetData<int16_t>(fmt)[idx];
	case LogicalTypeId::INTEGER:
		return UnifiedVectorFormat::GetData<int32_t>(fmt)[idx];
	case LogicalTypeId::UTINYINT:
		return UnifiedVectorFormat::GetData<uint8_t>(fmt)[idx];
	case LogicalTypeId::USMALLINT:
		return UnifiedVectorFormat::GetData<uint16_t>(fmt)[idx];
	case LogicalTypeId::BIGINT:
		return WideTagValue<int64_t>(t, idx);
	case LogicalTypeId::UINTEGER:
		return WideTagValue<uint32_t>(t, idx);
	case LogicalTypeId::UBIGINT:
		return WideTagValue<uint64_t>(t, idx);
	default:
		throw InternalException("UBAM tag: unsupported integer type '%s'", fdata.tags[t].type.ToString());
	}
}

void UBAMRecordEncoder::Read(idx_t row) {
	// read_id -> QNAME, unchanged. NULL is an error: a record needs a name.
	if (!read_id_data.validity.RowIsValid(read_id_data.sel->get_index(row))) {
		throw InvalidInputException("NULL value in read_id column");
	}
	ResolveSequenceRecordId(read_id_data, row, fdata.read_id_type, id_buf, id_ptr, id_size);
	CheckBAMReadName(id_ptr, id_size);

	auto seq_row = seq1_data.sel->get_index(row);
	if (!seq1_data.validity.RowIsValid(seq_row)) {
		ThrowNullColumn("sequence1", ReadId());
	}
	auto &seq = UnifiedVectorFormat::GetData<string_t>(seq1_data)[seq_row];
	seq_ptr = seq.GetData();
	seq_size = seq.GetSize();
	if (seq_size == 0) {
		throw InvalidInputException("Column 'sequence1' is empty for read '%s': a BAM record without bases reads as "
		                            "one whose sequence was left out",
		                            ReadId());
	}

	// qual1: raw Phred scores, one per base, none NULL.
	auto qual_row = qual1.data.sel->get_index(row);
	if (!qual1.data.validity.RowIsValid(qual_row)) {
		ThrowNullColumn("qual1", ReadId());
	}
	auto &entry = qual1.Check(qual_row, seq_size, "qual1", "", [&]() { return ReadId(); });
	qual_ptr = qual1.scores + entry.offset;
	uint8_t max_score = 0;
	for (idx_t k = 0; k < seq_size; k++) {
		max_score = MaxValue(max_score, qual_ptr[k]);
	}
	// Past what SAM text can show; a first score of 255 would also mark the BAM record as having no scores at all.
	if (max_score > miint::MAX_PHRED_QUALITY) {
		auto first =
		    *std::find_if(qual_ptr, qual_ptr + seq_size, [](uint8_t q) { return q > miint::MAX_PHRED_QUALITY; });
		throw InvalidInputException(
		    "Column 'qual1' has a quality score of %d for read '%s', past %d, the highest SAM text can show",
		    static_cast<int64_t>(first), ReadId(), static_cast<int64_t>(miint::MAX_PHRED_QUALITY));
	}

	for (idx_t t = 0; t < fdata.tags.size(); t++) {
		auto tag_idx = tag_data[t].sel->get_index(row);
		tag_present[t] = tag_data[t].validity.RowIsValid(tag_idx);
		if (tag_present[t]) {
			tag_values[t] = TagValue(t, tag_idx);
		}
	}
}

void UBAMRecordEncoder::Encode(idx_t row, RecordBuffer &out) {
	Read(row);
	// Unmapped (flag=4), no reference/cigar, mapq unavailable (255). QUAL is raw Phred, exactly what bam_set1 expects
	// (and what read_fastx/read_sequences_sam emit for qual1).
	if (bam_set1(record.get(), id_size, id_ptr, BAM_FUNMAP, /*tid*/ -1, /*pos*/ -1, /*mapq*/ 255, /*n_cigar*/ 0,
	             /*cigar*/ nullptr, /*mtid*/ -1, /*mpos*/ -1, /*isize*/ 0, seq_size, seq_ptr,
	             reinterpret_cast<const char *>(qual_ptr), /*l_aux*/ 0) < 0) {
		throw IOException("Failed to build uBAM record for read: " + ReadId());
	}

	// Per-record RG:Z (matches the @RG ID) so pbbam-based tools (lima) accept it.
	if (fdata.has_read_group) {
		if (bam_aux_append(record.get(), "RG", 'Z', NumericCast<int>(fdata.rg_id.size() + 1),
		                   reinterpret_cast<const uint8_t *>(fdata.rg_id.c_str())) < 0) {
			throw IOException("Failed to append RG tag to uBAM record");
		}
	}

	// Per-read integer aux tags.
	for (idx_t t = 0; t < fdata.tags.size(); t++) {
		if (!tag_present[t]) {
			continue;
		}
		if (bam_aux_append(record.get(), fdata.tags[t].tag.c_str(), 'i', sizeof(int32_t),
		                   reinterpret_cast<const uint8_t *>(&tag_values[t])) < 0) {
			throw IOException("Failed to append %s tag to uBAM record", fdata.tags[t].tag);
		}
	}

	StoreBAMRecord(*record, out);
}

//===--------------------------------------------------------------------===//
// Prepare and Flush (DuckDB's batch API)
//===--------------------------------------------------------------------===//
static unique_ptr<PreparedBatchData> UBAMCopyPrepareBatch(ClientContext &context, FunctionData &bind_data,
                                                          GlobalFunctionData &gstate_p,
                                                          unique_ptr<ColumnDataCollection> rows) {
	UBAMRecordEncoder encoder(bind_data.Cast<UBAMCopyBindData>());
	return RecordPrepareBatch(context, encoder, std::move(rows));
}

// DuckDB calls flush_batch for one batch at a time per file, in input order when order is kept.
static void UBAMCopyFlushBatch(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                               PreparedBatchData &batch_p) {
	auto &fdata = bind_data.Cast<UBAMCopyBindData>();
	auto &gstate = gstate_p.Cast<UBAMCopyGlobalState>();
	auto &batch = batch_p.Cast<RecordBatchData>();
	lock_guard<mutex> glock(gstate.lock);
	RecordFlushBatch(
	    context, batch, [&]() { return UBAMRecordEncoder(fdata); },
	    [&](const_data_ptr_t data, idx_t size) {
		    WriteBAMRecord(gstate.sam_file.get(), gstate.header.get(), data, size, "UBAM", gstate.file_path);
	    });
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void UBAMCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<UBAMCopyGlobalState>();
	lock_guard<mutex> glock(gstate.lock);
	CloseSAMCopyFile(gstate.sam_file.release(), "UBAM", gstate.file_path);
}

//===--------------------------------------------------------------------===//
// Registration
//===--------------------------------------------------------------------===//
// FORMAT UBAM is two registrations: the public one ("ubam") declares the options, the extension and a plan callback,
// which rebinds the statement as the writer ("_miint_ubam"). See SequenceCopyPlan.
static constexpr const char *UBAM_WRITER_FORMAT = "_miint_ubam";

// READ_GROUP and TAGS are structs, checked by ParseReadGroup and ParseTags. REFERENCE_LENGTHS is
// declared only to be refused with the reason, which DuckDB's unknown-option error would not give.
static void UBAMCopyOptions(ClientContext &context, CopyOptionsInput &input) {
	DeclareCopyOption(input, "read_group", LogicalType::ANY);
	DeclareCopyOption(input, "tags", LogicalType::ANY);
	DeclareCompressionLevelCopyOption(input);
	DeclareCopyOption(input, "reference_lengths", LogicalType::ANY);
}

static BoundStatement UBAMCopyPlan(Binder &binder, CopyStatement &stmt) {
	return SequenceCopyPlan(binder, stmt, UBAM_WRITER_FORMAT);
}

// A uBAM is a BAM file, so files DuckDB names end in .bam. A .bam path with no FORMAT still binds FORMAT BAM: DuckDB
// takes a missing format from the path's extension, by name.
CopyFunction CopyUBAMFunction::GetFunction() {
	CopyFunction function("ubam");
	function.copy_options = UBAMCopyOptions;
	function.plan = UBAMCopyPlan;
	function.extension = "bam";
	return function;
}

// FORMAT _miint_ubam: the writer. Named directly, it writes without the default batch cap.
static CopyFunction GetWriterFunction() {
	CopyFunction function(UBAM_WRITER_FORMAT);
	function.copy_options = UBAMCopyOptions;
	function.copy_to_bind = UBAMCopyBind;
	function.copy_to_initialize_global = UBAMCopyInitializeGlobal;
	function.copy_to_finalize = UBAMCopyFinalize;
	function.execution_mode = SequenceCopyExecutionMode;
	function.prepare_batch = UBAMCopyPrepareBatch;
	function.flush_batch = UBAMCopyFlushBatch;
	function.extension = "bam";
	return function;
}

void CopyUBAMFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
	loader.RegisterFunction(GetWriterFunction());
}

} // namespace duckdb
