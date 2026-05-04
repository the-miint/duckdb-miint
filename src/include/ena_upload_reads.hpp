// SPDX-License-Identifier: MIT
//
// `ena_upload_reads` table function. Streams a multi-sample reads relation
// to ENA's Webin upload area, returning per-file metadata (sample, filename,
// filetype, md5, bytes, layout). Transport is selected by `target_url`:
// `aspera://...` uses `ascp --mode=send`; `file://...` writes locally
// (used by the in-tree SQL test bypass).
//
// The encode → gzip → MD5 pipeline is shared across transports; only the
// final write target differs. Helpers (layout detection, ascp argv,
// URL parsing) are in `ena_upload_helpers.hpp` and unit-tested separately.

#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class ENAUploadReadsTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
