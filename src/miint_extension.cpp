#define DUCKDB_EXTENSION_MAIN

#include "miint_extension.hpp"
#include <alignment_flag_functions.hpp>
#include <alignment_functions.hpp>
#include <compress_intervals.hpp>
#include <copy_biom.hpp>
#include <copy_fasta.hpp>
#include <copy_fastq.hpp>
#include <copy_sam.hpp>
#include <kseq++/seqio.hpp>
#include <read_fastx.hpp>
#include <read_alignments.hpp>
#include <read_biom.hpp>
#include <miint_macros.hpp>
#include <sequence_functions.hpp>
#include <woltka_batched.hpp>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include <hdf5.h>

namespace fs = std::filesystem;

namespace duckdb {

static void SetDependencyLogging() {
	// HTSlib and HDF5 had runtime logging behaviors. Disable globally
	// for now. It's unclear whether these should be exposed or the
	// exact benefit, we will defer that decision for the future.
	hts_set_log_level(HTS_LOG_ERROR);
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
}

static void LoadInternal(ExtensionLoader &loader) {
	// TODO: use [[nodiscard]] throughout in headers
	// TODO: //! comment on headers
	ReadFastxTableFunction::Register(loader);
	ReadAlignmentsTableFunction::Register(loader);
	ReadBIOMTableFunction::Register(loader);

	AlignmentFlagFunctions::Register(loader);
	AlignmentSeqIdentityFunction::Register(loader);
	AlignmentQueryLengthFunction::Register(loader);
	AlignmentQueryCoverageFunction::Register(loader);
	CompressIntervalsFunction::Register(loader);
	SequenceFunctions::Register(loader);

	CopyBiomFunction::Register(loader);
	CopyFastqFunction::Register(loader);
	CopyFastaFunction::Register(loader);
	CopySAMFunction::Register(loader);

	MIINTMacros::Register(loader);
	WoltkaOguPerSampleBatchedTableFunction::Register(loader);
}

void MiintExtension::Load(ExtensionLoader &loader) {
	SetDependencyLogging();
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
