#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// COPY ... TO 'x.bam' (FORMAT UBAM, ...) -- writes an unaligned reads BAM (uBAM):
// SEQ/QUAL from read_id/sequence1/qual1, no @SQ header, an optional settable @RG,
// and optional per-read integer aux tags. The reads-side counterpart of FORMAT
// FASTQ. See docs/writing.md and issue #156.
class CopyUBAMFunction {
public:
	static CopyFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
