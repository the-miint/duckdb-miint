#pragma once

#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/execution/partition_info.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/parallel/pipeline.hpp"

namespace duckdb {

// Batch indexes for a multi-file reader that claims whole files in increasing order and never puts two files in one
// chunk. A chunk's index is its file's ordinal times a stride, plus the chunk's ordinal within the file. DuckDB's own
// readers take one index per claimed unit from a counter; for a whole file that would make the file one batch, which
// an order-keeping sink holds whole.
//
// DuckDB's rules for the indexes, and which consumers use them: docs/internals/duckdb-engine-notes.md, "Output order
// costs threads and memory". The stride keeps every index below the bound; files claimed in increasing order keep
// each thread's indexes from decreasing or falling below the lowest one held; and each file's own range keeps two
// threads from sharing one.
struct ScanBatchIndex {
	// The chunks per file that get an index of their own, for `file_count` files, so that every index stays below
	// DuckDB's limit.
	static idx_t Stride(idx_t file_count) {
		return (PipelineBuildState::BATCH_INCREMENT - 2) / MaxValue<idx_t>(file_count, 1);
	}

	// The index of the `chunk`-th chunk of the `file`-th file. Chunks past the stride repeat the file's last index,
	// so they join one batch instead of reaching into the next file's.
	static idx_t Index(idx_t file, idx_t chunk, idx_t stride) {
		return file * stride + MinValue<idx_t>(chunk, stride - 1);
	}

	// Called from Bind: keeps `function`'s batch index only for a scan of several files. A one-file scan reads on one
	// thread, where a batch index gains nothing and costs an order-kept CREATE TABLE AS (the engine notes above).
	// Not pursued: grouping a file's chunks into about 64 MiB per index, so that short-read batches fill a row group.
	static void KeepForSeveralFiles(BoundTableFunction &function, idx_t file_count) {
		if (file_count < 2) {
			function.get_partition_data = nullptr;
		}
	}

	// A get_partition_data result: `index`, the batch index of the chunk just emitted. Partition columns are another
	// kind of partition data, which these readers do not provide.
	static OperatorPartitionData PartitionData(const TableFunctionGetPartitionInput &input, idx_t index) {
		if (input.partition_info.RequiresPartitionColumns()) {
			throw NotImplementedException("This table function does not provide partition columns");
		}
		return OperatorPartitionData(index);
	}
};

} // namespace duckdb
