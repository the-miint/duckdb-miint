#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// predict_genes(relation, genome_id := 'col', method := 'prodigal', ...): prokaryotic gene calls, one row per gene,
// for each genome (the distinct non-NULL values of the genome_id column) in a relation of contigs (read_id,
// sequence1). Prodigal runs in the gpl-boundary daemon (>= 0.5.0). See docs/gene_prediction.md.
class PredictGenesTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
