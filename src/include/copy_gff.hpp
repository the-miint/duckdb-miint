#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// COPY ... TO 'x.gff' (FORMAT GFF). Two dialects, chosen by DIALECT:
//   - 'gff3' (default), from read_gff's columns (seqid, source, type, position, stop_position, score, strand, phase,
//     attributes): a line per row;
//   - 'prodigal', from predict_genes' columns: Prodigal's own GFF, with each contig's "# Sequence Data" and
//     "# Model Data" lines, as `prodigal -f gff` writes it.
// The rows become lines in SQL (the _miint_gff3_lines / _miint_prodigal_gff_lines macros), and DuckDB's CSV writer
// writes them, as DuckDB's own FORMAT JSON does: compression, PARTITION_BY and the rest are the CSV writer's.
class CopyGffFunction {
public:
	static CopyFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
