// read_sequences_sam's batch index as DuckDB reads it: what get_partition_data returns after each chunk, where
// max_batch_bytes closes a chunk, and how many threads a scan asks for. The SQL tests check the order an order-keeping
// sink restores, which one batch per file would restore as well. That batching would make a sink hold each file whole,
// which giving each chunk its own index avoids: so this pins one batch index per chunk, counted from 0 in each file,
// within that file's own range, a file's chunks past the stride kept in that range, and chunks no bigger than
// max_batch_bytes unless they hold one read.

#include <catch2/catch_test_macros.hpp>

#include "read_sequences_sam.hpp"
#include "scan_batch_index.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <string>
#include <vector>

namespace {

using namespace duckdb;

// `count` uBAM files of 300 reads of 200 bases each. A read counts its name, 200 bases and 200 qualities, 402 to 404
// bytes, so 16 KiB chunks hold about 40 reads and each file spans several.
struct UbamFiles {
	std::filesystem::path dir;
	std::vector<std::string> paths;

	UbamFiles(Connection &con, int count) {
		dir = std::filesystem::temp_directory_path() /
		      ("miint_read_sequences_sam_batch_index_" +
		       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directories(dir);
		for (int file = 0; file < count; file++) {
			auto path = (dir / ("f" + std::to_string(file) + ".bam")).string();
			auto written =
			    con.Query("COPY (SELECT 'r' || i AS read_id, repeat('ACGT', 50) AS sequence1, "
			              "list_transform(range(200), lambda j: 30::UTINYINT) AS qual1 FROM range(300) t(i)) "
			              "TO '" +
			              path + "' (FORMAT UBAM)");
			INFO(written->GetError());
			REQUIRE_FALSE(written->HasError());
			paths.push_back(path);
		}
	}

	~UbamFiles() {
		std::error_code ignored;
		std::filesystem::remove_all(dir, ignored);
	}
};

// One emitted chunk: its batch index, and the bytes max_batch_bytes counts (each read's name, bases and qualities).
struct Chunk {
	idx_t batch_index;
	idx_t rows;
	idx_t bytes;
	idx_t first_read_bytes;
};

idx_t ReadBytes(DataChunk &output, idx_t row) {
	auto read_id = output.GetValue(1, row);
	auto sequence1 = output.GetValue(3, row);
	auto qual1 = output.GetValue(5, row);
	return StringValue::Get(read_id).size() + StringValue::Get(sequence1).size() +
	       (qual1.IsNull() ? 0 : ListValue::GetChildren(qual1).size());
}

// Every chunk a scan emits, by file. One local state reads the files in turn, as one thread does when it claims another
// file.
std::vector<std::vector<Chunk>> ReadChunks(ClientContext &context, ReadSequencesSamTableFunction::Data &data,
                                           ReadSequencesSamTableFunction::GlobalState &global) {
	ReadSequencesSamTableFunction::LocalState local;
	vector<LogicalType> types(data.types.begin(), data.types.end());
	DataChunk output;
	output.Initialize(Allocator::DefaultAllocator(), types);

	std::vector<std::vector<Chunk>> chunks(data.file_paths.size());
	while (true) {
		output.Reset();
		TableFunctionInput input(&data, &local, &global);
		ReadSequencesSamTableFunction::Execute(context, input, output);
		if (output.size() == 0) {
			break;
		}
		auto partition_info = OperatorPartitionInfo::BatchIndex();
		TableFunctionGetPartitionInput partition_input(&data, &local, &global, partition_info);
		Chunk chunk {ReadSequencesSamTableFunction::GetPartitionData(context, partition_input).batch_index,
		             output.size(), 0, 0};
		for (idx_t row = 0; row < output.size(); row++) {
			auto read_bytes = ReadBytes(output, row);
			if (row == 0) {
				chunk.first_read_bytes = read_bytes;
			}
			chunk.bytes += read_bytes;
		}
		chunks[local.current_file_idx].push_back(chunk);
	}
	return chunks;
}

} // namespace

TEST_CASE("read_sequences_sam gives each chunk its own batch index, counted from 0 in its file's range, and closes a "
          "chunk before the read that would pass max_batch_bytes",
          "[read_sequences_sam][batch_index]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &context = *con.context;
	UbamFiles files(con, 2);
	const idx_t max_batch_bytes = 16 * 1024;

	ReadSequencesSamTableFunction::Data data(files.paths, false, false, max_batch_bytes);
	ReadSequencesSamTableFunction::GlobalState global(files.paths, FileSystem::GetFileSystem(context), false, 1);
	auto chunks = ReadChunks(context, data, global);

	const auto stride = ScanBatchIndex::Stride(files.paths.size());
	for (idx_t file = 0; file < files.paths.size(); file++) {
		INFO("file " << file);
		auto &file_chunks = chunks[file];
		// At the default of 64 MiB each file would be one chunk of 300 reads.
		REQUIRE(file_chunks.size() > 1);
		idx_t rows = 0;
		for (idx_t c = 0; c < file_chunks.size(); c++) {
			INFO("chunk " << c);
			REQUIRE(file_chunks[c].batch_index == file * stride + c);
			// Reads of 402 to 404 bytes put many in each chunk, and the chunk stays at most the limit: a read that
			// would pass it waits for the next chunk. (Only a chunk of one read may pass it.)
			REQUIRE(file_chunks[c].rows > 1);
			REQUIRE(file_chunks[c].bytes <= max_batch_bytes);
			// Every chunk but the file's last closed because its next read would have passed the limit: none closed
			// early.
			if (c + 1 < file_chunks.size()) {
				REQUIRE(file_chunks[c].bytes + file_chunks[c + 1].first_read_bytes > max_batch_bytes);
			}
			rows += file_chunks[c].rows;
		}
		REQUIRE(rows == 300);
	}
}

// The clamp keeps a file with more chunks than the stride inside its own range: its later chunks join one batch
// instead of taking the next file's indexes. No test file comes near a real stride (about 5 * 10^12 chunks for two
// files), so this one is set to 3.
TEST_CASE("read_sequences_sam gives a file's chunks past the stride its last index",
          "[read_sequences_sam][batch_index]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &context = *con.context;
	UbamFiles files(con, 2);

	ReadSequencesSamTableFunction::Data data(files.paths, false, false, 16 * 1024);
	ReadSequencesSamTableFunction::GlobalState global(files.paths, FileSystem::GetFileSystem(context), false, 1);
	global.batch_stride = 3;
	auto chunks = ReadChunks(context, data, global);

	for (idx_t file = 0; file < files.paths.size(); file++) {
		INFO("file " << file);
		REQUIRE(chunks[file].size() > 3);
		for (idx_t c = 0; c < chunks[file].size(); c++) {
			INFO("chunk " << c);
			REQUIRE(chunks[file][c].batch_index == file * 3 + std::min<idx_t>(c, 2));
		}
	}
}

TEST_CASE("read_sequences_sam asks for one thread per file", "[read_sequences_sam][batch_index]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &fs = FileSystem::GetFileSystem(*con.context);
	UbamFiles files(con, 10);

	// Ten files, so ten threads: the count follows the files.
	ReadSequencesSamTableFunction::GlobalState ten(files.paths, fs, false, 1);
	REQUIRE(ten.MaxThreads() == 10);

	// stdin is always a scan's only path (Bind refuses it in a list), so it is this case.
	ReadSequencesSamTableFunction::GlobalState one({files.paths[0]}, fs, false, 1);
	REQUIRE(one.MaxThreads() == 1);
}
