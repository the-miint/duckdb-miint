#pragma once
#ifdef MIINT_HAS_SYLPH

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

// Pull sylph.h directly so SylphGenomeSketchParams is available as a field type
// on Data below — Bind() seeds it via sylph_genome_sketch_params_default() and
// the named-parameter dispatch only overrides individual fields the caller set,
// keeping sylph as the single source of truth for the sketch defaults.
#include "sylph.h"

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// Build a sylph `.syldb` reference database from a table/view of reference
// sequences (read_fastx output or equivalent), grouping contigs into genomes by
// a genome-key column. Wraps the sylph index-builder FFI
// (sylph_index_builder_*), which sketches each genome from in-memory sequence
// bytes and serializes the Vec<GenomeSketch> to disk. Returns a single status
// row.
//
// Input contract (validated at bind time via ValidateSequenceTableSchema):
//   read_id   VARCHAR|BIGINT (required) — contig name; the first contig of each
//                                         genome supplies its first_contig_name
//   sequence1 VARCHAR        (required) — the contig sequence
//   <genome_id> (required named param)  — grouping column; one GenomeSketch per
//                                         distinct value, all its contigs merged
//
// Contigs within a genome are fed in `order_by` order (default: read_id) so the
// resulting sketch is reproducible; pass order_by := 'sequence_index' to match
// original FASTA order. Given the same genomes + params, the resulting sketches
// are content-identical to `sylph sketch` (verified in test/sql/sylph_index_create.test
// and ext/sylph's builder tests) — but not byte-identical: sylph orders genomes
// in the file by parallel-completion, so its output isn't byte-deterministic.
//
// Like rype_index_create, the build is a synchronous side effect performed in
// InitGlobal; Execute only emits the status row. Single-threaded.
class SylphIndexCreateTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string source_table;
		std::string output_path;
		std::string genome_id_col;            // grouping column (required)
		std::string order_by_col = "read_id"; // within-genome contig order

		// Source has a `comment` column (read_fastx does). If so, the contig name
		// is rebuilt as the full header `read_id || ' ' || comment` to match
		// sylph's first_contig_name; otherwise it's just read_id.
		bool has_comment = false;

		// Seeded from sylph defaults in Bind(); named params override individual
		// fields. 0 fields mean "use sylph's default" (k=31, c=200, spacing=30).
		SylphGenomeSketchParams sketch_params;

		vector<std::string> names;
		vector<LogicalType> types;

		Data()
		    : names({"output_path", "k", "c", "num_genomes", "status"}),
		      types({LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::UBIGINT,
		             LogicalType::VARCHAR}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		// The build happens in InitGlobal; Execute emits one status row and
		// `done` guards against a second emit.
		bool done = false;
		size_t num_genomes = 0;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	struct LocalState : public LocalTableFunctionState {};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb

#endif // MIINT_HAS_SYLPH
