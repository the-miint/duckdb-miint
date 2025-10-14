#define DUCKDB_EXTENSION_MAIN

#include "miint_extension.hpp"
#include <alignment_flag_functions.hpp>
#include <alignment_functions.hpp>
#include <kseq++/seqio.hpp>
#include <read_fastx.hpp>
#include <read_sam.hpp>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

namespace fs = std::filesystem;

namespace duckdb {

static void LoadInternal(ExtensionLoader &loader) {
	// TODO: use [[nodiscard]] throughout in headers
	// TODO: //! comment on headers
	ReadFastxTableFunction read_fastx;
	loader.RegisterFunction(read_fastx.GetFunction());

	ReadSAMTableFunction read_sam;
	loader.RegisterFunction(read_sam.GetFunction());

	AlignmentFlagFunctions::Register(loader);
	AlignmentSeqIdentityFunction::Register(loader);

	// QualFilterScalarFunction find_low_quality_window;
	// loader.RegisterFunction(find_low_quality_window.GetFunction());
}

void MiintExtension::Load(ExtensionLoader &loader) {
	hts_set_log_level(HTS_LOG_ERROR);

	LoadInternal(loader);
}
std::string MiintExtension::Name() {
	return "miint";
}

std::string MiintExtension::Version() const {
#ifdef EXT_VERSION_MIINT
	return EXT_VERSION_MIINT;
#else
	return "unversioned";
#endif
}

} // namespace duckdb

extern "C" {

DUCKDB_CPP_EXTENSION_ENTRY(miint, loader) {
	duckdb::LoadInternal(loader);
}
}
