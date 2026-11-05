#pragma once

#include "sra_run.hpp"

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <memory>
#include <string>
#include <vector>

namespace duckdb {

// read_sra(path): the reads of local .sra files as fasterq-dump writes them, in read_fastx's
// schema. The scan is parallel within a run: a run is split into work units of consecutive spots,
// decoded by any thread, and emitted in spot order through their batch indexes.
class ReadSraTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> paths;
		miint::sra::SplitMode split = miint::sra::SplitMode::Split3;
		bool include_technical = false;
		bool include_filepath = false;
		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);

	// The accession fasterq-dump prints for a local file: its name without a .sra or .sralite
	// extension. Exposed for the NCBI fetcher, which names runs by accession directly.
	static std::string AccessionFromPath(const std::string &path);
};

} // namespace duckdb
