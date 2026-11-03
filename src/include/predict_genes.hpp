#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// The Prodigal predict_genes is built for. predict_genes refuses a gpl-boundary whose --version lists another. A
// literal too, for prodigal_gff_lines' SQL (miint_macros.hpp).
#define PRODIGAL_VERSION_LITERAL "2.6.5"
constexpr const char *kProdigalVersion = PRODIGAL_VERSION_LITERAL;

// predict_genes(relation, genome_id := 'col', method := 'prodigal', ...): prokaryotic gene calls, one row per gene,
// for each genome (the distinct non-NULL values of the genome_id column) in a relation of contigs (read_id,
// sequence1). Prodigal runs in the gpl-boundary daemon (>= 0.5.0). See docs/gene_prediction.md.
class PredictGenesTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
