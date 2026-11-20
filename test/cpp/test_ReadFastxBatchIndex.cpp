// read_fastx's batch index as DuckDB reads it: what get_partition_data returns after each chunk. The SQL tests check
// the order an order-keeping sink restores, which one batch per file would restore as well. That batching would make
// a sink hold each file whole, which the design avoids: so this pins one batch index per chunk, counted from 0 in
// each file, within that file's own range.

#include <catch2/catch_test_macros.hpp>

#include "read_fastx.hpp"
#include "scan_batch_index.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <chrono>
#include <filesystem>
#include <string>
#include <vector>

namespace {

using namespace duckdb;

struct TwoFastqFiles {
	std::filesystem::path dir;
	std::vector<std::string> paths;

	explicit TwoFastqFiles(Connection &con) {
		dir = std::filesystem::temp_directory_path() /
		      ("miint_read_fastx_batch_index_" +
		       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directories(dir);
		// 300 records of 200 bases each, so 16 KB chunks hold about 40 records, and each file spans several.
		for (int file = 0; file < 2; file++) {
			auto path = (dir / ("f" + std::to_string(file) + ".fastq")).string();
			auto written =
			    con.Query("COPY (SELECT 'r' || i AS read_id, repeat('ACGT', 50) AS sequence1, "
			              "list_transform(range(200), lambda j: 30::UTINYINT) AS qual1 FROM range(300) t(i)) "
			              "TO '" +
			              path + "' (FORMAT FASTQ)");
			INFO(written->GetError());
			REQUIRE_FALSE(written->HasError());
			paths.push_back(path);
		}
	}

	~TwoFastqFiles() {
		std::error_code ignored;
		std::filesystem::remove_all(dir, ignored);
	}
};

} // namespace

TEST_CASE("read_fastx gives each chunk its own batch index, counted from 0 in its file's range",
          "[read_fastx][batch_index]") {
	DuckDB db(nullptr);
	Connection con(db);
	auto &context = *con.context;
	TwoFastqFiles files(con);

	// One local state reads both files in turn, as one thread does when it claims a second file.
	ReadFastxTableFunction::Data data(files.paths, std::nullopt, false, false, 33, 16 * 1024, false);
	ReadFastxTableFunction::GlobalState global(FileSystem::GetFileSystem(context), files.paths, std::nullopt, false,
	                                           false, false);
	ReadFastxTableFunction::LocalState local;
	vector<LogicalType> types(data.types.begin(), data.types.end());
	DataChunk output;
	output.Initialize(Allocator::DefaultAllocator(), types);

	std::vector<std::vector<idx_t>> indexes(files.paths.size());
	while (true) {
		output.Reset();
		TableFunctionInput input(&data, &local, &global);
		ReadFastxTableFunction::Execute(context, input, output);
		if (output.size() == 0) {
			break;
		}
		auto partition_info = OperatorPartitionInfo::BatchIndex();
		TableFunctionGetPartitionInput partition_input(&data, &local, &global, partition_info);
		indexes[local.current_file_idx].push_back(
		    ReadFastxTableFunction::GetPartitionData(context, partition_input).batch_index);
	}

	const auto stride = ScanBatchIndex::Stride(files.paths.size());
	for (idx_t file = 0; file < files.paths.size(); file++) {
		INFO("file " << file);
		REQUIRE(indexes[file].size() > 1);
		for (idx_t chunk = 0; chunk < indexes[file].size(); chunk++) {
			INFO("chunk " << chunk);
			REQUIRE(indexes[file][chunk] == file * stride + chunk);
		}
	}
}
