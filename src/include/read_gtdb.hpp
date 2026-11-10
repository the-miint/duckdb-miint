#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// read_gtdb_metadata(source := NULL, release := NULL, refresh := false) reads GTDB's
// genome metadata: a domain column ('Bacteria' or 'Archaea', from the file's name), then
// every column of the file, VARCHAR, with GTDB's 'none' read as NULL.
// read_gtdb_qc_failed(...) reads GTDB's list of genomes that failed its quality checks,
// every column VARCHAR, values as written.
//
// The user-facing contract is docs/insdc_ncbi.md, "Read GTDB's genome metadata".
class ReadGtdbFunctions {
public:
	// Where GTDB serves its releases: the default of the miint_gtdb_url setting.
	static constexpr const char *DEFAULT_URL = "https://data.gtdb.ecogenomic.org/releases";

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
