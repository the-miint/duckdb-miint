#include "copy_ubam.hpp"
#include "copy_format_common.hpp" // ResolveSequenceRecordId (read_id VARCHAR/BIGINT/UUID dispatch)
#include "htslib_raii.hpp"
#include "id_column_utils.hpp" // IsAllowedIdType, AllowedIdTypeList
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include <htslib-1.22.1/htslib/hts.h>
#include <htslib-1.22.1/htslib/sam.h>

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
		string field = StringUtil::Upper(child_types[i].first);
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
		const string &tag = child_types[i].first;
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
			throw BinderException("TAGS column '%s' (for tag '%s') must be an integer type", col, tag);
		}
		out.push_back(UBAMTagColumn {tag, col_idx, sql_types[col_idx]});
	}
}

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> UBAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                             const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<UBAMCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;

	// Locate required columns (case-insensitive, matching the SAM/BAM writer).
	for (idx_t i = 0; i < names.size(); i++) {
		if (StringUtil::CIEquals(names[i], "read_id")) {
			result->read_id_idx = i;
		} else if (StringUtil::CIEquals(names[i], "sequence1")) {
			result->sequence1_idx = i;
		} else if (StringUtil::CIEquals(names[i], "qual1")) {
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
		if (StringUtil::CIEquals(option.first, "read_group")) {
			ParseReadGroup(option.second[0], result->rg_line, result->rg_id);
			result->has_read_group = true;
		} else if (StringUtil::CIEquals(option.first, "tags")) {
			ParseTags(option.second[0], names, result->tags, sql_types);
		} else if (StringUtil::CIEquals(option.first, "compression_level")) {
			result->compression_level = option.second[0].GetValue<int32_t>();
			if (result->compression_level < 0 || result->compression_level > 9) {
				throw BinderException("COMPRESSION_LEVEL must be between 0 and 9, got %d", result->compression_level);
			}
		} else if (StringUtil::CIEquals(option.first, "reference_lengths")) {
			throw BinderException("COPY FORMAT UBAM does not accept REFERENCE_LENGTHS: a uBAM is headerless "
			                      "(unaligned reads, no @SQ). Use FORMAT BAM for aligned records.");
		} else {
			throw BinderException("Unknown option for COPY FORMAT UBAM: %s", option.first);
		}
	}

	return std::move(result);
}

//===--------------------------------------------------------------------===//
// Global / Local State
//===--------------------------------------------------------------------===//
struct UBAMCopyGlobalState : public GlobalFunctionData {
	mutex write_lock;
	SAMFilePtr sam_file;
	SAMHeaderPtr header;
	std::atomic<bool> header_written {false};
};

static unique_ptr<GlobalFunctionData> UBAMCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                               const string &file_path) {
	auto &fdata = bind_data.Cast<UBAMCopyBindData>();
	auto gstate = make_uniq<UBAMCopyGlobalState>();

	string mode = fdata.compression_level >= 0 ? "wb" + std::to_string(fdata.compression_level) : "wb";
	gstate->sam_file = SAMFilePtr(sam_open(file_path.c_str(), mode.c_str()));
	if (!gstate->sam_file) {
		throw IOException("Failed to open file for writing: " + file_path);
	}

	// Offload bgzf compression to a worker pool (records still emitted in order).
	auto n_threads = NumericCast<int>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	if (n_threads > 1) {
		hts_set_threads(gstate->sam_file.get(), n_threads);
	}

	// Headerless-by-design: @HD + optional @RG, NO @SQ (these reads are unaligned).
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

	return std::move(gstate);
}

struct UBAMCopyLocalState : public LocalFunctionData {
	BAMRecordPtr record;
	string id_buf; // reused for BIGINT/UUID read_id stringification

	UBAMCopyLocalState() {
		record = BAMRecordPtr(bam_init1());
		if (!record) {
			throw IOException("Failed to allocate BAM record");
		}
	}
};

static unique_ptr<LocalFunctionData> UBAMCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	return make_uniq<UBAMCopyLocalState>();
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
// Read one integer cell as int64, dispatching on the column's physical type.
static int64_t ReadIntCell(const UnifiedVectorFormat &fmt, const LogicalType &type, idx_t idx) {
	switch (type.id()) {
	case LogicalTypeId::TINYINT:
		return UnifiedVectorFormat::GetData<int8_t>(fmt)[idx];
	case LogicalTypeId::SMALLINT:
		return UnifiedVectorFormat::GetData<int16_t>(fmt)[idx];
	case LogicalTypeId::INTEGER:
		return UnifiedVectorFormat::GetData<int32_t>(fmt)[idx];
	case LogicalTypeId::BIGINT:
		return UnifiedVectorFormat::GetData<int64_t>(fmt)[idx];
	case LogicalTypeId::UTINYINT:
		return UnifiedVectorFormat::GetData<uint8_t>(fmt)[idx];
	case LogicalTypeId::USMALLINT:
		return UnifiedVectorFormat::GetData<uint16_t>(fmt)[idx];
	case LogicalTypeId::UINTEGER:
		return UnifiedVectorFormat::GetData<uint32_t>(fmt)[idx];
	case LogicalTypeId::UBIGINT:
		return NumericCast<int64_t>(UnifiedVectorFormat::GetData<uint64_t>(fmt)[idx]);
	default:
		throw InternalException("UBAM tag: unsupported integer type '%s'", type.ToString());
	}
}

static void UBAMCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                         LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<UBAMCopyBindData>();
	auto &gstate = gstate_p.Cast<UBAMCopyGlobalState>();
	auto &lstate = lstate_p.Cast<UBAMCopyLocalState>();
	const idx_t n = input.size();

	UnifiedVectorFormat read_id_data, seq1_data, qual1_data;
	input.data[fdata.read_id_idx].ToUnifiedFormat(n, read_id_data);
	input.data[fdata.sequence1_idx].ToUnifiedFormat(n, seq1_data);
	input.data[fdata.qual1_idx].ToUnifiedFormat(n, qual1_data);

	auto seq1_strings = UnifiedVectorFormat::GetData<string_t>(seq1_data);
	auto qual1_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual1_data);
	auto qual1_child = FlatVector::GetData<uint8_t>(ListVector::GetEntry(input.data[fdata.qual1_idx]));

	// Tag columns: unify once per chunk.
	vector<UnifiedVectorFormat> tag_formats(fdata.tags.size());
	for (idx_t t = 0; t < fdata.tags.size(); t++) {
		input.data[fdata.tags[t].col_idx].ToUnifiedFormat(n, tag_formats[t]);
	}

	// Write the header once (first thread to arrive).
	if (!gstate.header_written.load()) {
		lock_guard<mutex> glock(gstate.write_lock);
		if (!gstate.header_written.load()) {
			if (sam_hdr_write(gstate.sam_file.get(), gstate.header.get()) < 0) {
				throw IOException("Failed to write uBAM header");
			}
			gstate.header_written.store(true);
		}
	}

	for (idx_t row = 0; row < n; row++) {
		// read_id -> QNAME (verbatim; opaque). NULL is an error: a record needs a name.
		auto id_row = read_id_data.sel->get_index(row);
		if (!read_id_data.validity.RowIsValid(id_row)) {
			throw InvalidInputException("NULL value in read_id column (row %llu)", row);
		}
		const char *id_ptr;
		idx_t id_size;
		ResolveSequenceRecordId(read_id_data, row, fdata.read_id_type, lstate.id_buf, id_ptr, id_size);

		// sequence1.
		auto seq_row = seq1_data.sel->get_index(row);
		if (!seq1_data.validity.RowIsValid(seq_row)) {
			throw InvalidInputException("NULL value in sequence1 column (row %llu)", row);
		}
		const char *seq_ptr = seq1_strings[seq_row].GetData();
		idx_t seq_size = seq1_strings[seq_row].GetSize();

		// qual1 (raw Phred bytes, must match sequence length).
		auto qual_row = qual1_data.sel->get_index(row);
		if (!qual1_data.validity.RowIsValid(qual_row)) {
			throw InvalidInputException("NULL value in qual1 column (row %llu)", row);
		}
		idx_t qual_len = qual1_entries[qual_row].length;
		const uint8_t *qual_ptr = qual1_child + qual1_entries[qual_row].offset;
		if (qual_len != seq_size) {
			throw InvalidInputException(
			    "Quality score length (%llu) does not match sequence length (%llu) for row %llu", qual_len, seq_size,
			    row);
		}

		// Build the record thread-locally (no lock): unmapped (flag=4), no
		// reference/cigar, mapq unavailable (255). QUAL is raw Phred, exactly what
		// bam_set1 expects (and what read_fastx/read_sequences_sam emit for qual1).
		if (bam_set1(lstate.record.get(), id_size, id_ptr, BAM_FUNMAP, /*tid*/ -1, /*pos*/ -1, /*mapq*/ 255,
		             /*n_cigar*/ 0, /*cigar*/ nullptr, /*mtid*/ -1, /*mpos*/ -1, /*isize*/ 0, seq_size, seq_ptr,
		             reinterpret_cast<const char *>(qual_ptr), /*l_aux*/ 0) < 0) {
			throw IOException("Failed to build uBAM record for read: " + string(id_ptr, id_size));
		}

		// Per-record RG:Z (matches the @RG ID) so pbbam-based tools (lima) accept it.
		if (fdata.has_read_group) {
			if (bam_aux_append(lstate.record.get(), "RG", 'Z', NumericCast<int>(fdata.rg_id.size() + 1),
			                   reinterpret_cast<const uint8_t *>(fdata.rg_id.c_str())) < 0) {
				throw IOException("Failed to append RG tag to uBAM record");
			}
		}

		// Per-read integer aux tags. A NULL cell omits the tag for that record.
		for (idx_t t = 0; t < fdata.tags.size(); t++) {
			auto &tag = fdata.tags[t];
			auto &fmt = tag_formats[t];
			auto tag_idx = fmt.sel->get_index(row);
			if (!fmt.validity.RowIsValid(tag_idx)) {
				continue;
			}
			// BAM 'i' aux tags are int32. A value outside that range fails loud here
			// (checked cast) rather than silently truncating as copy_sam.cpp does.
			int32_t v = NumericCast<int32_t>(ReadIntCell(fmt, tag.type, tag_idx));
			if (bam_aux_append(lstate.record.get(), tag.tag.c_str(), 'i', sizeof(int32_t),
			                   reinterpret_cast<uint8_t *>(&v)) < 0) {
				throw IOException("Failed to append %s tag to uBAM record", tag.tag);
			}
		}

		// Only the actual write touches shared state (htslib samFile is not
		// thread-safe for concurrent writes); the record was built lock-free above.
		{
			lock_guard<mutex> glock(gstate.write_lock);
			if (sam_write1(gstate.sam_file.get(), gstate.header.get(), lstate.record.get()) < 0) {
				throw IOException("Failed to write uBAM record");
			}
		}
	}
}

//===--------------------------------------------------------------------===//
// Combine / Finalize
//===--------------------------------------------------------------------===//
static void UBAMCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                            LocalFunctionData &lstate_p) {
	// HTSlib buffers internally; nothing to combine.
}

static void UBAMCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<UBAMCopyGlobalState>();
	// Empty input: still emit a valid (header-only) BAM.
	if (!gstate.header_written.load()) {
		if (sam_hdr_write(gstate.sam_file.get(), gstate.header.get()) < 0) {
			throw IOException("Failed to write uBAM header for empty file");
		}
		gstate.header_written.store(true);
	}
	// File closed by the smart pointer.
}

//===--------------------------------------------------------------------===//
// Registration
//===--------------------------------------------------------------------===//
CopyFunction CopyUBAMFunction::GetFunction() {
	CopyFunction function("ubam");
	function.copy_to_bind = UBAMCopyBind;
	function.copy_to_initialize_global = UBAMCopyInitializeGlobal;
	function.copy_to_initialize_local = UBAMCopyInitializeLocal;
	function.copy_to_sink = UBAMCopySink;
	function.copy_to_combine = UBAMCopyCombine;
	function.copy_to_finalize = UBAMCopyFinalize;
	// Deliberately no `function.extension`: a .bam path auto-detects as FORMAT BAM
	// (aligned). uBAM is a special case that must be requested via FORMAT UBAM.
	return function;
}

void CopyUBAMFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
