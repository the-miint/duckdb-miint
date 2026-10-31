#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// COPY ... TO 'x.gff' (FORMAT GFF). Two inputs:
//   - read_gff's columns (seqid, source, type, position, stop_position, score, strand, phase, attributes): GFF3, a
//     line per row;
//   - predict_genes' columns (detected by gene_index): Prodigal's own GFF, with each contig's "# Sequence Data" and
//     "# Model Data" lines, as `prodigal -f gff` writes it.
class CopyGffFunction {
public:
	static CopyFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
