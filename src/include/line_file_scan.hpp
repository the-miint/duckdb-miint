#pragma once

#include "duckdb_seq_stream.hpp"
#include "line_splitter.hpp"
#include "sam_style_tags.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace duckdb {

// Shared machinery for table functions that turn each line of a text file into at most
// one row (read_paf, read_gfa). Files go through DuckDB's file system, so gzip,
// multi-member gzip and remote paths share one code path. Each file is read start to end
// by one thread; several files are read in parallel.

struct LineScanData : public TableFunctionData {
	std::string function_name;
	std::vector<std::string> file_paths;
	bool include_filepath = false;
};

// A chunk stops filling once its lines hold this many bytes, so a file of megabase records
// (GFA segments that carry their sequence) is not held in memory thousands of records at a time.
constexpr idx_t LINE_SCAN_MAX_CHUNK_BYTES = 16 * 1024 * 1024;

// Resolve a reader's arguments: a VARCHAR path or glob, or a VARCHAR[] of literal paths,
// plus include_filepath. Rejects NULL, an empty list and stdin with messages prefixed by the
// function's name, and a local path that does not exist with "File not found: <path>", as
// read_fastx does. `names` and `return_types` are the reader's own columns; a filepath column
// is appended when requested.
unique_ptr<FunctionData> BindLineScan(ClientContext &context, TableFunctionBindInput &input, vector<std::string> &names,
                                      vector<LogicalType> &return_types);

// A line-scan reader's table function: one path argument, include_filepath, the shared
// global state, and a row estimate from the files' sizes.
TableFunction MakeLineScanFunction(const std::string &name, table_function_t execute, table_function_bind_t bind,
                                   table_function_init_local_t init_local);

struct LineScanGlobalState : public GlobalTableFunctionState {
	std::mutex lock;
	size_t next_file_idx = 0;
	idx_t max_threads;

	explicit LineScanGlobalState(idx_t max_threads) : max_threads(max_threads) {
	}
	idx_t MaxThreads() const override {
		return max_threads;
	}
};

// The file a thread is reading.
struct LineFile {
	std::string path;
	miint::DuckDBSeqStreamHandle stream;
	// Holds `stream`'s raw pointer, taken when it is built, so it must be declared after
	// `stream` (members are built in order), and `stream` must not be reset or replaced.
	miint::LineSplitter lines;

	LineFile(FileSystem &fs, const std::string &path, const std::string &function_name);
};

// A reader's per-thread state derives from this. The file is destroyed at its end, so open
// descriptors scale with threads rather than with files.
struct LineScanLocalState : public LocalTableFunctionState {
	std::unique_ptr<LineFile> file;
};

// Write `tags` as row `row` of a MAP(VARCHAR, VARCHAR) vector.
void SetTagMap(Vector &map_vec, idx_t row, const std::vector<miint::SamStyleTag> &tags);

// Fill `output` with rows from the calling thread's files. `parse(line)` turns one line into
// the reader's current record and returns false for a line that holds none (a comment);
// a std::runtime_error from it is reported as "<function>: <file>: line N: <message>".
// `write(row)` then stores that record at `row`. A chunk holds one file's rows only, so the
// filepath column, the last one when requested, is constant.
template <typename Parse, typename Write>
void ScanLines(ClientContext &context, TableFunctionInput &data_p, DataChunk &output, Parse &&parse, Write &&write) {
	auto &data = data_p.bind_data->Cast<LineScanData>();
	auto &global = data_p.global_state->Cast<LineScanGlobalState>();
	auto &local = data_p.local_state->Cast<LineScanLocalState>();

	idx_t row = 0;
	idx_t chunk_bytes = 0;
	std::string chunk_path = local.file ? local.file->path : std::string();
	while (row < STANDARD_VECTOR_SIZE && chunk_bytes < LINE_SCAN_MAX_CHUNK_BYTES) {
		if (!local.file) {
			if (row > 0) {
				break;
			}
			size_t idx;
			{
				std::lock_guard<std::mutex> guard(global.lock);
				if (global.next_file_idx >= data.file_paths.size()) {
					break;
				}
				idx = global.next_file_idx++;
			}
			local.file =
			    make_uniq<LineFile>(FileSystem::GetFileSystem(context), data.file_paths[idx], data.function_name);
			chunk_path = local.file->path;
		}

		std::string_view line;
		if (!local.file->lines.Next(line)) {
			local.file.reset();
			continue;
		}
		bool is_record;
		try {
			is_record = parse(line);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException(data.function_name + ": " + local.file->path + ": line " +
			                            std::to_string(local.file->lines.LineNumber()) + ": " + e.what());
		}
		if (is_record) {
			write(row);
			row++;
			chunk_bytes += line.size();
		}
	}
	if (data.include_filepath && row > 0) {
		SetResultVectorFilepath(output.data[output.ColumnCount() - 1], chunk_path);
	}
	output.SetCardinality(row);
}

} // namespace duckdb
