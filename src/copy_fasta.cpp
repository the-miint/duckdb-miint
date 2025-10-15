#include "copy_fasta.hpp"
#include "copy_format_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include <sstream>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct FastaCopyBindData : public FunctionData {
	bool interleave = false;
	bool id_as_sequence_index = false;
	bool include_comment = false;
	bool use_compression = false;
	string file_path;
	bool is_paired = false;
	vector<string> names;
	
	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<FastaCopyBindData>();
		result->interleave = interleave;
		result->id_as_sequence_index = id_as_sequence_index;
		result->include_comment = include_comment;
		result->use_compression = use_compression;
		result->file_path = file_path;
		result->is_paired = is_paired;
		result->names = names;
		return std::move(result);
	}
	
	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<FastaCopyBindData>();
		return interleave == other.interleave &&
		       id_as_sequence_index == other.id_as_sequence_index &&
		       include_comment == other.include_comment &&
		       use_compression == other.use_compression &&
		       file_path == other.file_path &&
		       is_paired == other.is_paired;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> FastaCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                              const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<FastaCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;
	
	// Detect columns
	ColumnIndices indices;
	indices.FindIndices(names);
	
	bool has_sequence1 = indices.sequence1_idx != DConstants::INVALID_INDEX;
	bool has_sequence2 = indices.sequence2_idx != DConstants::INVALID_INDEX;
	bool has_read_id = indices.read_id_idx != DConstants::INVALID_INDEX;
	bool has_sequence_index = indices.sequence_index_idx != DConstants::INVALID_INDEX;
	
	// Validate required columns
	ValidateRequiredColumns(has_read_id, has_sequence1, "FASTA");
	
	result->is_paired = has_sequence2;
	
	// Parse common parameters
	CommonCopyParameters common_params;
	bool has_interleave_param = false;
	
	for (auto &option : input.info.options) {
		if (StringUtil::CIEquals(option.first, "interleave")) {
			has_interleave_param = true;
		} else if (!StringUtil::CIEquals(option.first, "id_as_sequence_index") &&
		           !StringUtil::CIEquals(option.first, "include_comment") &&
		           !StringUtil::CIEquals(option.first, "compression")) {
			throw BinderException("Unknown option for COPY FORMAT FASTA: %s", option.first);
		}
	}
	
	common_params.ParseFromOptions(input.info.options, result->file_path);
	
	result->interleave = common_params.interleave;
	result->id_as_sequence_index = common_params.id_as_sequence_index;
	result->include_comment = common_params.include_comment;
	result->use_compression = common_params.use_compression;
	
	// Validate paired-end parameters
	ValidatePairedEndParameters(result->is_paired, has_interleave_param, result->interleave, result->file_path);
	
	// Validate sequence_index parameter
	ValidateSequenceIndexParameter(result->id_as_sequence_index, has_sequence_index);
	
	return std::move(result);
}

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct FastaCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	unique_ptr<CopyFileHandle> file_r1;
	unique_ptr<CopyFileHandle> file_r2;
	bool is_paired = false;
	bool interleave = false;
};

static unique_ptr<GlobalFunctionData> FastaCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                                 const string &file_path) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	auto &fs = FileSystem::GetFileSystem(context);
	
	auto gstate = make_uniq<FastaCopyGlobalState>();
	gstate->is_paired = fdata.is_paired;
	gstate->interleave = fdata.interleave;
	
	if (fdata.is_paired && !fdata.interleave) {
		// Split mode: open two files
		string path_r1 = SubstituteOrientation(file_path, "R1");
		string path_r2 = SubstituteOrientation(file_path, "R2");
		
		gstate->file_r1 = make_uniq<CopyFileHandle>(fs, path_r1, fdata.use_compression);
		gstate->file_r2 = make_uniq<CopyFileHandle>(fs, path_r2, fdata.use_compression);
	} else {
		// Interleaved or single-end: open one file
		gstate->file_r1 = make_uniq<CopyFileHandle>(fs, file_path, fdata.use_compression);
	}
	
	return std::move(gstate);
}

//===--------------------------------------------------------------------===//
// Local State
//===--------------------------------------------------------------------===//
struct FastaCopyLocalState : public LocalFunctionData {};

static unique_ptr<LocalFunctionData> FastaCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	return make_uniq<FastaCopyLocalState>();
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void WriteFastaRecord(CopyFileHandle &file, const string &id, const string &seq, const string &comment) {
	stringstream record;
	record << ">" << id;
	if (!comment.empty()) {
		record << " " << comment;
	}
	record << "\n" << seq << "\n";
	
	file.Write(record.str());
}

static void FastaCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                          LocalFunctionData &lstate, DataChunk &input) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	auto &gstate = gstate_p.Cast<FastaCopyGlobalState>();
	
	lock_guard<mutex> glock(gstate.lock);
	
	// Get column indices
	ColumnIndices indices;
	indices.FindIndices(fdata.names);
	
	// Process each row
	UnifiedVectorFormat read_id_data, sequence_index_data, comment_data;
	UnifiedVectorFormat seq1_data, seq2_data;
	
	input.data[indices.read_id_idx].ToUnifiedFormat(input.size(), read_id_data);
	auto read_ids = UnifiedVectorFormat::GetData<string_t>(read_id_data);
	
	if (fdata.id_as_sequence_index) {
		input.data[indices.sequence_index_idx].ToUnifiedFormat(input.size(), sequence_index_data);
	}
	
	if (fdata.include_comment && indices.comment_idx != DConstants::INVALID_INDEX) {
		input.data[indices.comment_idx].ToUnifiedFormat(input.size(), comment_data);
	}
	
	input.data[indices.sequence1_idx].ToUnifiedFormat(input.size(), seq1_data);
	auto seq1_strings = UnifiedVectorFormat::GetData<string_t>(seq1_data);
	
	bool is_paired = fdata.is_paired;
	if (is_paired) {
		input.data[indices.sequence2_idx].ToUnifiedFormat(input.size(), seq2_data);
	}
	
	for (idx_t row = 0; row < input.size(); row++) {
		auto row_idx = read_id_data.sel->get_index(row);
		
		// Get ID
		string id;
		if (fdata.id_as_sequence_index) {
			auto seq_idx_data = UnifiedVectorFormat::GetData<int64_t>(sequence_index_data);
			auto seq_idx_row = sequence_index_data.sel->get_index(row);
			id = to_string(seq_idx_data[seq_idx_row]);
		} else {
			id = read_ids[row_idx].GetString();
		}
		
		// Get comment
		string comment;
		if (fdata.include_comment && indices.comment_idx != DConstants::INVALID_INDEX) {
			auto comment_strings = UnifiedVectorFormat::GetData<string_t>(comment_data);
			auto comment_row = comment_data.sel->get_index(row);
			if (comment_data.validity.RowIsValid(comment_row)) {
				comment = comment_strings[comment_row].GetString();
			}
		}
		
		// Get R1 data
		auto seq1_row = seq1_data.sel->get_index(row);
		string seq1 = seq1_strings[seq1_row].GetString();
		
		// Write R1 record
		WriteFastaRecord(*gstate.file_r1, id, seq1, comment);
		
		// Handle R2 for paired-end
		if (is_paired) {
			auto seq2_strings = UnifiedVectorFormat::GetData<string_t>(seq2_data);
			auto seq2_row = seq2_data.sel->get_index(row);
			string seq2 = seq2_strings[seq2_row].GetString();
			
			if (fdata.interleave) {
				// Write R2 to same file
				WriteFastaRecord(*gstate.file_r1, id, seq2, comment);
			} else {
				// Write R2 to separate file
				WriteFastaRecord(*gstate.file_r2, id, seq2, comment);
			}
		}
	}
}

//===--------------------------------------------------------------------===//
// Combine
//===--------------------------------------------------------------------===//
static void FastaCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate,
                             LocalFunctionData &lstate) {
	// No-op for FASTA
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void FastaCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<FastaCopyGlobalState>();
	lock_guard<mutex> glock(gstate.lock);
	
	if (gstate.file_r1) {
		gstate.file_r1->Close();
	}
	if (gstate.file_r2) {
		gstate.file_r2->Close();
	}
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
CopyFunction CopyFastaFunction::GetFunction() {
	CopyFunction func("fasta");
	func.copy_to_bind = FastaCopyBind;
	func.copy_to_initialize_local = FastaCopyInitializeLocal;
	func.copy_to_initialize_global = FastaCopyInitializeGlobal;
	func.copy_to_sink = FastaCopySink;
	func.copy_to_combine = FastaCopyCombine;
	func.copy_to_finalize = FastaCopyFinalize;
	return func;
}

} // namespace duckdb
