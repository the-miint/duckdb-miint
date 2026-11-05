#pragma once

namespace duckdb {
class ExtensionLoader;
//! Registers the replacement scan that reads a file named in FROM with the miint reader its
//! suffix names, as DuckDB reads FROM 'x.parquet': FROM 'reads.fq.gz' is read_fastx('reads.fq.gz').
void RegisterFileSuffixScan(ExtensionLoader &loader);
} // namespace duckdb
