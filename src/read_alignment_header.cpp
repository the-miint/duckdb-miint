#include "read_alignment_header.hpp"
#include "SAMReader.hpp"
#include "remote_file_helper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"

// htslib comes in via SAMReader.hpp, which owns the SAMFilePtr/SAMHeaderPtr RAII
// wrappers used below (note the version-prefixed include path this tree uses).

namespace duckdb {

unique_ptr<FunctionData> ReadAlignmentHeaderTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                vector<LogicalType> &return_types,
                                                                vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_alignment_header: path cannot be NULL");
	}
	data->path = input.inputs[0].ToString();
	if (data->path.empty()) {
		throw InvalidInputException("read_alignment_header: path cannot be empty");
	}

	// Local paths only, by design. read_alignments streams remote files through
	// hfile_duckdb_open; that plumbing is deliberately not duplicated here, so say
	// so plainly instead of failing deeper with an htslib open error.
	if (miint::RemoteFileHelper::IsRemotePath(data->path)) {
		throw InvalidInputException("read_alignment_header: remote paths are not supported ('%s'). Download the file "
		                            "first, or read its records with read_alignments, which does stream remotely.",
		                            data->path);
	}

	auto &fs = FileSystem::GetFileSystem(context);
	if (!fs.FileExists(data->path)) {
		throw IOException("File not found: " + data->path);
	}

	names = {"tid", "reference", "length"};
	// length is BIGINT, from sam_hdr_tid2len (hts_pos_t). A 32-bit read is wrong:
	// data/sam/foo_large_positions.bam carries LN=3000000000, which is negative as
	// an int32.
	return_types = {LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::BIGINT};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadAlignmentHeaderTableFunction::InitGlobal(ClientContext &context,
                                                                                  TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Raw htslib rather than SAMReader: SAMReader keeps its sam_hdr_t private with
	// no accessor, and its require_references modes would reject a zero-@SQ header.
	// This mirrors the header probe already in read_alignments.cpp.
	miint::SAMFilePtr fp(sam_open(data.path.c_str(), "r"));
	if (!fp) {
		throw IOException("Failed to open alignment file: " + data.path);
	}
	miint::SAMHeaderPtr hdr(sam_hdr_read(fp.get()));
	if (!hdr) {
		throw IOException("Failed to read alignment header from: " + data.path);
	}

	// A header with no @SQ lines is legitimate -- an unaligned BAM/SAM has no
	// references, including the uBAM files miint itself writes -- so this yields
	// zero rows rather than raising.
	int n_ref = sam_hdr_nref(hdr.get());
	gstate->entries.reserve(n_ref > 0 ? static_cast<size_t>(n_ref) : 0);
	for (int tid = 0; tid < n_ref; tid++) {
		const char *name = sam_hdr_tid2name(hdr.get(), tid);
		if (!name) {
			throw IOException("Alignment header in '%s' has no name for reference index %d", data.path, tid);
		}
		hts_pos_t len = sam_hdr_tid2len(hdr.get(), tid);
		if (len <= 0) {
			// htslib documents tid2len as strictly positive on success.
			throw IOException("Alignment header in '%s' reports a non-positive length for reference '%s'", data.path,
			                  name);
		}
		gstate->entries.push_back(SQEntry {tid, std::string(name), static_cast<int64_t>(len)});
	}

	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState>
ReadAlignmentHeaderTableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void ReadAlignmentHeaderTableFunction::Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, gstate.entries.size() - gstate.cursor);

	auto tids = FlatVector::GetData<int32_t>(output.data[0]);
	auto lengths = FlatVector::GetData<int64_t>(output.data[2]);
	for (idx_t i = 0; i < count; i++) {
		const auto &e = gstate.entries[gstate.cursor + i];
		tids[i] = e.tid;
		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], e.reference);
		lengths[i] = e.length;
	}

	gstate.cursor += count;
	output.SetCardinality(count);
}

TableFunction ReadAlignmentHeaderTableFunction::GetFunction() {
	return TableFunction("read_alignment_header", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
}

void ReadAlignmentHeaderTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
