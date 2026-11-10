#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// read_gtdb_metadata(source := NULL, release := NULL, refresh := false) reads GTDB's
// genome metadata: a domain column ('Bacteria' or 'Archaea', from the file's name), then
// every column of the file, VARCHAR, with GTDB's 'none' read as NULL.
// read_gtdb_qc_failed(...) reads GTDB's list of genomes that failed its quality checks,
// every column VARCHAR, values as written.
//
// Pass exactly one of source (a file, a URL, or a directory holding one release's files)
// and release ('232', '214.1'). A release's files are downloaded into a cache when the
// query runs, and read from there after. Binding reads only the files' header lines,
// from GTDB for a release not yet cached; the rows are read by DuckDB's CSV reader.
class ReadGtdbFunctions {
public:
	// Where GTDB serves its releases: the default of the miint_gtdb_url setting.
	static constexpr const char *DEFAULT_URL = "https://data.gtdb.ecogenomic.org/releases";

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
