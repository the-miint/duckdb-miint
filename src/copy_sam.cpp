#include "copy_sam.hpp"
#include "copy_format_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct SAMCopyBindData : public FunctionData {
	bool include_header = true;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	idx_t flush_size = 1024 * 1024;
	string file_path;
	vector<string> names;
	std::unordered_map<string, int64_t> reference_lengths;

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<SAMCopyBindData>();
		result->include_header = include_header;
		result->compression = compression;
		result->flush_size = flush_size;
		result->file_path = file_path;
		result->names = names;
		result->reference_lengths = reference_lengths;
		return std::move(result);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<SAMCopyBindData>();
		return include_header == other.include_header && compression == other.compression &&
		       file_path == other.file_path && reference_lengths == other.reference_lengths;
	}
};

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
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> SAMCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                            const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<SAMCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;

	// Detect columns
	SAMColumnIndices indices;
	indices.FindIndices(names);

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
	for (auto &option : input.info.options) {
		if (StringUtil::CIEquals(option.first, "include_header")) {
			result->include_header = option.second[0].GetValue<bool>();
		} else if (StringUtil::CIEquals(option.first, "compression")) {
			result->compression = DetectCompressionType(result->file_path, option.second[0]);
		} else if (StringUtil::CIEquals(option.first, "reference_lengths")) {
			const auto &map_value = option.second[0];
			if (map_value.type().id() != LogicalTypeId::MAP) {
				throw BinderException("reference_lengths must be a MAP type (e.g., MAP{'ref1': 100, 'ref2': 200})");
			}
			const auto &entries = MapValue::GetChildren(map_value);
			for (const auto &entry : entries) {
				const auto &kv = StructValue::GetChildren(entry);
				const auto key = StringValue::Get(kv[0]);
				const auto value = kv[1].GetValue<int64_t>();
				if (value < 0) {
					throw BinderException("reference_lengths values must be non-negative");
				}
				result->reference_lengths.emplace(key, value);
			}
		} else {
			throw BinderException("Unknown option for COPY FORMAT SAM: %s", option.first);
		}
	}
	
	// Validate reference_lengths requirement
	if (result->include_header && result->reference_lengths.empty()) {
		throw BinderException("COPY FORMAT SAM with INCLUDE_HEADER=true requires REFERENCE_LENGTHS parameter "
		                      "(e.g., REFERENCE_LENGTHS=MAP{'chr1': 248956422, 'chr2': 242193529})");
	}

	// Auto-detect compression from file extension if not specified
	if (result->compression == FileCompressionType::UNCOMPRESSED) {
		bool compression_specified = false;
		for (auto &option : input.info.options) {
			if (StringUtil::CIEquals(option.first, "compression")) {
				compression_specified = true;
				break;
			}
		}
		if (!compression_specified) {
			result->compression = DetectCompressionType(result->file_path, Value());
		}
	}

	return std::move(result);
}

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct SAMCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	unique_ptr<CopyFileHandle> file;
	bool header_written = false;
};

static unique_ptr<GlobalFunctionData> SAMCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                               const string &file_path) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &fs = FileSystem::GetFileSystem(context);

	auto gstate = make_uniq<SAMCopyGlobalState>();
	gstate->file = make_uniq<CopyFileHandle>(fs, file_path, fdata.compression);

	return std::move(gstate);
}

//===--------------------------------------------------------------------===//
// Local State
//===--------------------------------------------------------------------===//
struct SAMCopyLocalState : public LocalFunctionData {
	unique_ptr<FormatWriterState> writer_state;
};

static unique_ptr<LocalFunctionData> SAMCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto lstate = make_uniq<SAMCopyLocalState>();
	
	lstate->writer_state = make_uniq<FormatWriterState>(context.client, fdata.flush_size);
	
	return std::move(lstate);
}

//===--------------------------------------------------------------------===//
// Helper Functions
//===--------------------------------------------------------------------===//
static void AppendTag(stringstream &ss, const char *tag_name, int64_t value, bool &first_tag) {
	if (!first_tag) {
		ss << "\t";
	}
	ss << tag_name << ":i:" << value;
	first_tag = false;
}

static void AppendTag(stringstream &ss, const char *tag_name, const string &value, bool &first_tag) {
	if (!first_tag) {
		ss << "\t";
	}
	ss << tag_name << ":Z:" << value;
	first_tag = false;
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void SAMCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                        LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();
	auto &lstate = lstate_p.Cast<SAMCopyLocalState>();

	// Get column indices
	SAMColumnIndices indices;
	indices.FindIndices(fdata.names);

	// Process each row - update local reference lengths
	UnifiedVectorFormat read_id_data, flags_data, reference_data, position_data, stop_position_data, mapq_data;
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

	if (indices.stop_position_idx != DConstants::INVALID_INDEX) {
		input.data[indices.stop_position_idx].ToUnifiedFormat(input.size(), stop_position_data);
	}
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

	auto stop_position_ptr =
	    indices.stop_position_idx != DConstants::INVALID_INDEX
	        ? UnifiedVectorFormat::GetData<int64_t>(stop_position_data)
	        : nullptr;
	auto tag_as_ptr = indices.tag_as_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_as_data)
	                      : nullptr;
	auto tag_xs_ptr = indices.tag_xs_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_xs_data)
	                      : nullptr;
	auto tag_ys_ptr = indices.tag_ys_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_ys_data)
	                      : nullptr;
	auto tag_xn_ptr = indices.tag_xn_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_xn_data)
	                      : nullptr;
	auto tag_xm_ptr = indices.tag_xm_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_xm_data)
	                      : nullptr;
	auto tag_xo_ptr = indices.tag_xo_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_xo_data)
	                      : nullptr;
	auto tag_xg_ptr = indices.tag_xg_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_xg_data)
	                      : nullptr;
	auto tag_nm_ptr = indices.tag_nm_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<int64_t>(tag_nm_data)
	                      : nullptr;
	auto tag_yt_ptr = indices.tag_yt_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<string_t>(tag_yt_data)
	                      : nullptr;
	auto tag_md_ptr = indices.tag_md_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<string_t>(tag_md_data)
	                      : nullptr;
	auto tag_sa_ptr = indices.tag_sa_idx != DConstants::INVALID_INDEX
	                      ? UnifiedVectorFormat::GetData<string_t>(tag_sa_data)
	                      : nullptr;

	// Get reference to local buffer
	auto &stream = *lstate.writer_state->stream;
	
	// Build records in local buffer - NO LOCK
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
		string cigar = cigar_ptr[cigar_idx].GetString();
		string mate_reference = mate_reference_ptr[mate_reference_idx].GetString();
		int64_t mate_position = mate_position_ptr[mate_position_idx];
		int64_t template_length = template_length_ptr[template_length_idx];

		// Build SAM record: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
		stringstream record;
		record << read_id << "\t" << flags << "\t" << reference << "\t" << position << "\t"
		       << static_cast<int>(mapq) << "\t" << cigar << "\t" << mate_reference << "\t" << mate_position
		       << "\t" << template_length << "\t*\t*";
		// Add optional tags
		if (tag_as_ptr && indices.tag_as_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_as_data.sel->get_index(i);
			if (tag_as_data.validity.RowIsValid(tag_idx)) {
				record << "\tAS:i:" << tag_as_ptr[tag_idx];
			}
		}
		if (tag_xs_ptr && indices.tag_xs_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_xs_data.sel->get_index(i);
			if (tag_xs_data.validity.RowIsValid(tag_idx)) {
				record << "\tXS:i:" << tag_xs_ptr[tag_idx];
			}
		}
		if (tag_ys_ptr && indices.tag_ys_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_ys_data.sel->get_index(i);
			if (tag_ys_data.validity.RowIsValid(tag_idx)) {
				record << "\tYS:i:" << tag_ys_ptr[tag_idx];
			}
		}
		if (tag_xn_ptr && indices.tag_xn_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_xn_data.sel->get_index(i);
			if (tag_xn_data.validity.RowIsValid(tag_idx)) {
				record << "\tXN:i:" << tag_xn_ptr[tag_idx];
			}
		}
		if (tag_xm_ptr && indices.tag_xm_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_xm_data.sel->get_index(i);
			if (tag_xm_data.validity.RowIsValid(tag_idx)) {
				record << "\tXM:i:" << tag_xm_ptr[tag_idx];
			}
		}
		if (tag_xo_ptr && indices.tag_xo_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_xo_data.sel->get_index(i);
			if (tag_xo_data.validity.RowIsValid(tag_idx)) {
				record << "\tXO:i:" << tag_xo_ptr[tag_idx];
			}
		}
		if (tag_xg_ptr && indices.tag_xg_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_xg_data.sel->get_index(i);
			if (tag_xg_data.validity.RowIsValid(tag_idx)) {
				record << "\tXG:i:" << tag_xg_ptr[tag_idx];
			}
		}
		if (tag_nm_ptr && indices.tag_nm_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_nm_data.sel->get_index(i);
			if (tag_nm_data.validity.RowIsValid(tag_idx)) {
				record << "\tNM:i:" << tag_nm_ptr[tag_idx];
			}
		}
		if (tag_yt_ptr && indices.tag_yt_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_yt_data.sel->get_index(i);
			if (tag_yt_data.validity.RowIsValid(tag_idx)) {
				record << "\tYT:Z:" << tag_yt_ptr[tag_idx].GetString();
			}
		}
		if (tag_md_ptr && indices.tag_md_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_md_data.sel->get_index(i);
			if (tag_md_data.validity.RowIsValid(tag_idx)) {
				record << "\tMD:Z:" << tag_md_ptr[tag_idx].GetString();
			}
		}
		if (tag_sa_ptr && indices.tag_sa_idx != DConstants::INVALID_INDEX) {
			auto tag_idx = tag_sa_data.sel->get_index(i);
			if (tag_sa_data.validity.RowIsValid(tag_idx)) {
				record << "\tSA:Z:" << tag_sa_ptr[tag_idx].GetString();
			}
		}

		record << "\n";
		
		// Write record to local buffer
		string record_str = record.str();
		stream.WriteData(const_data_ptr_cast(record_str.c_str()), record_str.size());
		lstate.writer_state->written_anything = true;
	}

	// Check if we need to flush (buffer exceeded threshold)
	if (lstate.writer_state->stream->GetPosition() >= lstate.writer_state->flush_size) {
		// Acquire lock only for writing
		lock_guard<mutex> glock(gstate.lock);
		
		// Write header if needed (only first flush)
		if (!gstate.header_written && fdata.include_header) {
			stringstream header;
			for (const auto &ref : fdata.reference_lengths) {
				header << "@SQ\tSN:" << ref.first << "\tLN:" << ref.second << "\n";
			}
			gstate.file->WriteString(header.str());
			gstate.header_written = true;
		} else if (!gstate.header_written && !fdata.include_header) {
			gstate.header_written = true;
		}
		
		// Flush buffer (already holding lock, so write directly)
		if (lstate.writer_state->written_anything) {
			gstate.file->Write(lstate.writer_state->stream->GetData(), lstate.writer_state->stream->GetPosition());
			lstate.writer_state->Reset();
		}
	}
}

//===--------------------------------------------------------------------===//
// Combine
//===--------------------------------------------------------------------===//
static void SAMCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                           LocalFunctionData &lstate_p) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();
	auto &lstate = lstate_p.Cast<SAMCopyLocalState>();
	
	// Acquire lock for header check and flush
	lock_guard<mutex> glock(gstate.lock);
	
	// Write header if needed (only first combine)
	if (!gstate.header_written && fdata.include_header) {
		stringstream header;
		for (const auto &ref : fdata.reference_lengths) {
			header << "@SQ\tSN:" << ref.first << "\tLN:" << ref.second << "\n";
		}
		gstate.file->WriteString(header.str());
		gstate.header_written = true;
	} else if (!gstate.header_written && !fdata.include_header) {
		gstate.header_written = true;
	}
	
	// Flush any remaining data in local buffer (already holding lock)
	if (lstate.writer_state->written_anything) {
		gstate.file->Write(lstate.writer_state->stream->GetData(), lstate.writer_state->stream->GetPosition());
		lstate.writer_state->Reset();
	}
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void SAMCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &fdata = bind_data.Cast<SAMCopyBindData>();
	auto &gstate = gstate_p.Cast<SAMCopyGlobalState>();

	// Write header if we haven't written anything yet (empty file with header)
	if (!gstate.header_written && fdata.include_header) {
		stringstream header;
		for (const auto &ref : fdata.reference_lengths) {
			header << "@SQ\tSN:" << ref.first << "\tLN:" << ref.second << "\n";
		}
		gstate.file->WriteString(header.str());
	}

	gstate.file->Close();
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

} // namespace duckdb
