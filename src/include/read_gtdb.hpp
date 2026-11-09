#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// read_gtdb_metadata(source := NULL, release := NULL, refresh := false) reads GTDB's
// genome metadata: a domain column ('Bacteria' or 'Archaea', from the file's name), then
// every column of the file, VARCHAR, with GTDB's 'none' read as NULL.
// read_gtdb_qc_failed(...) reads GTDB's list of genomes that failed its quality checks,
// every column VARCHAR, values as written.
//
// Pass exactly one of source (a file, or a directory holding one release's files) and
// release. Binding reads only the files' header lines; the rows are read by DuckDB's CSV
// reader when the query runs.
class ReadGtdbFunctions {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
