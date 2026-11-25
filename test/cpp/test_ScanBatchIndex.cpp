// The batch index arithmetic of multi-file readers (scan_batch_index.hpp), checked against the rules DuckDB enforces
// (docs/internals/duckdb-engine-notes.md, "Output order costs threads and memory"): indexes that never decrease along
// files in claim order and their chunks, stay below BATCH_INCREMENT - 2, and give each file a range of its own.

#include <catch2/catch_test_macros.hpp>

#include "scan_batch_index.hpp"

#include <algorithm>
#include <vector>

using duckdb::idx_t;
using duckdb::PipelineBuildState;
using duckdb::ScanBatchIndex;

namespace {

// The count of indexes DuckDB accepts: 0 to BATCH_INCREMENT - 3.
constexpr idx_t INDEX_COUNT = PipelineBuildState::BATCH_INCREMENT - 2;

} // namespace

TEST_CASE("ScanBatchIndex: the stride splits the index range between the files", "[batch_index]") {
	REQUIRE(ScanBatchIndex::Stride(1) == INDEX_COUNT);
	REQUIRE(ScanBatchIndex::Stride(4) == INDEX_COUNT / 4);
	// No files is treated as one: Stride never divides by zero.
	REQUIRE(ScanBatchIndex::Stride(0) == INDEX_COUNT);
	for (idx_t files : {idx_t(1), idx_t(3), idx_t(4000), idx_t(596046), idx_t(10000000)}) {
		const auto stride = ScanBatchIndex::Stride(files);
		REQUIRE(stride >= 1);
		REQUIRE(stride * files <= INDEX_COUNT);
	}
}

TEST_CASE("ScanBatchIndex: an index is the file's block plus the chunk", "[batch_index]") {
	const auto stride = ScanBatchIndex::Stride(10);
	REQUIRE(ScanBatchIndex::Index(0, 0, stride) == 0);
	REQUIRE(ScanBatchIndex::Index(0, 5, stride) == 5);
	REQUIRE(ScanBatchIndex::Index(3, 0, stride) == 3 * stride);
	REQUIRE(ScanBatchIndex::Index(3, 7, stride) == 3 * stride + 7);
}

TEST_CASE("ScanBatchIndex: a file with more chunks than the stride repeats its last index", "[batch_index]") {
	const idx_t stride = 4;
	REQUIRE(ScanBatchIndex::Index(2, 3, stride) == 2 * stride + 3);
	REQUIRE(ScanBatchIndex::Index(2, 4, stride) == 2 * stride + 3);
	REQUIRE(ScanBatchIndex::Index(2, 1000, stride) == 2 * stride + 3);
	// Still below the next file's first index.
	REQUIRE(ScanBatchIndex::Index(2, 1000, stride) < ScanBatchIndex::Index(3, 0, stride));
}

TEST_CASE("ScanBatchIndex: indexes never decrease along files and chunks, and stay below DuckDB's limit",
          "[batch_index]") {
	for (idx_t files : {idx_t(1), idx_t(2), idx_t(7), idx_t(596046)}) {
		const auto stride = ScanBatchIndex::Stride(files);
		// Files in claim order, each once.
		std::vector<idx_t> claimed {0, files / 2, files - 1};
		claimed.erase(std::unique(claimed.begin(), claimed.end()), claimed.end());
		idx_t previous = 0;
		for (idx_t file : claimed) {
			for (idx_t chunk : {idx_t(0), idx_t(1), stride - 1, stride, idx_t(-1) / 2}) {
				const auto index = ScanBatchIndex::Index(file, chunk, stride);
				REQUIRE(index >= previous);
				REQUIRE(index < INDEX_COUNT);
				previous = index;
			}
		}
		// The last file's first index is above every index of the file before it.
		if (files > 1) {
			REQUIRE(ScanBatchIndex::Index(files - 1, 0, stride) > ScanBatchIndex::Index(files - 2, idx_t(-1), stride));
		}
	}
}
