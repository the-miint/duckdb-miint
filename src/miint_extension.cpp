#define DUCKDB_EXTENSION_MAIN

#include "miint_extension.hpp"
#include <alignment_flag_functions.hpp>
#include <alignment_functions.hpp>
#include <compress_intervals.hpp>
#include <copy_fasta.hpp>
#include <copy_fastq.hpp>
#include <copy_sam.hpp>
#include <kseq++/seqio.hpp>
#include <read_fastx.hpp>
#include <read_sam.hpp>
#include <miint_macros.hpp>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include <hdf5.h>

namespace fs = std::filesystem;

namespace duckdb {

static void poke() {
	// Open HDF5 file
	auto file_id = H5Fopen("foo.biom", H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0) {
		throw IOException("Failed to open HDF5 file");
	}
}

static void LoadInternal(ExtensionLoader &loader) {
	// TODO: use [[nodiscard]] throughout in headers
	// TODO: //! comment on headers
	ReadFastxTableFunction::Register(loader);
	ReadSAMTableFunction::Register(loader);
	AlignmentFlagFunctions::Register(loader);
	AlignmentSeqIdentityFunction::Register(loader);
	CompressIntervalsFunction::Register(loader);
	CopyFastqFunction::Register(loader);
	CopyFastaFunction::Register(loader);
	CopySAMFunction::Register(loader);
	MIINTMacros::Register(loader);
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
