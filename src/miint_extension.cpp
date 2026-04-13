#define DUCKDB_EXTENSION_MAIN

#include "miint_extension.hpp"
#include "aspera_utils.hpp"
#include <csignal>
#include <alignment_flag_functions.hpp>
#include <alignment_slice.hpp>
#include <alignment_functions.hpp>
#include <compress_intervals.hpp>
#include <compute_coverage_depth.hpp>
#include <copy_fasta.hpp>
#include <copy_fastq.hpp>
#include <copy_newick.hpp>
#include <copy_sam.hpp>
#include <read_jplace_newick.hpp>
#include <read_newick.hpp>
#include <tree_resolve_placement.hpp>
#include <kseq++/seqio.hpp>
#include <read_fastx.hpp>
#include <read_alignments.hpp>
#include <read_sequences_sam.hpp>
#include <read_sequences_sff.hpp>
#include <align_minimap2.hpp>
#include <align_minimap2_sharded.hpp>
#include <save_minimap2_index.hpp>
#include <read_ncbi_fasta.hpp>
#include <read_ncbi.hpp>
#include <read_ncbi_annotation.hpp>
#include <read_ena.hpp>
#include <read_ena_attributes.hpp>
#include <read_ena_sequences.hpp>
#include <miint_macros.hpp>
#include "duckdb/main/extension_helper.hpp"
#include "duckdb/main/database.hpp"
#include <sequence_functions.hpp>
#include <formula_function.hpp>
#include <massql_function.hpp>
#include <mzml_peak_pair_function.hpp>
#include <align_mafft.hpp>
#include <deblur_table_function.hpp>
#include <align_pairwise_functions.hpp>
#include <read_mzml.hpp>
#include <read_mzml_chromatograms.hpp>
#include <read_mzxml.hpp>
#include <rype_classify.hpp>
#include <rype_extract.hpp>
#include <rype_log_ratio.hpp>
#ifdef MIINT_HAS_VSEARCH
#include <uchime_ref.hpp>
#include <uchime_denovo.hpp>
#include <mask_dust_function.hpp>
#include <merge_pairs_function.hpp>
#include <search_sequences.hpp>
#include <cluster_sequences.hpp>
#endif
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include <htslib-1.22.1/htslib/hts.h>
#include <kseq++/config.hpp>
#include <zlib.h>

#ifdef MIINT_HAS_HDF5
#include <copy_biom.hpp>
#include <read_biom.hpp>
#include <hdf5.h>
#endif

#ifdef MIINT_HAS_BOWTIE2
#include <align_bowtie2.hpp>
#include <align_bowtie2_sharded.hpp>
#endif

namespace fs = std::filesystem;

namespace duckdb {

void SetDependencyLogging() {
	// HTSlib and HDF5 had runtime logging behaviors. Disable globally
	// for now. It's unclear whether these should be exposed or the
	// exact benefit, we will defer that decision for the future.
	hts_set_log_level(HTS_LOG_ERROR);
#ifdef MIINT_HAS_HDF5
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
#endif
}

void SetupSignalHandling() {
#if defined(MIINT_HAS_BOWTIE2) || MIINT_ASPERA_SUPPORTED
	// Ignore SIGPIPE globally so that writes to closed pipes return EPIPE instead of
	// killing the process. Needed for Bowtie2Aligner, AsperaProcess, and other subprocess
	// management where pipes may close unexpectedly.
	// Note: This is a PROCESS-WIDE setting that persists for the lifetime of the process.
	// Setting it once at extension load is thread-safe (vs calling signal() from multiple threads).
	std::signal(SIGPIPE, SIG_IGN);
#endif
}

static void MiintVersionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
#ifdef EXT_VERSION_MIINT
	result.Reference(Value(EXT_VERSION_MIINT));
#else
	result.Reference(Value("unversioned"));
#endif
}

struct MiintVersionsData : public TableFunctionData {
	vector<pair<string, string>> versions;
	bool done = false;
};

static unique_ptr<FunctionData> MiintVersionsBind(ClientContext &context, TableFunctionBindInput &input,
                                                  vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MiintVersionsData>();
	names = {"library", "version"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR};

#ifdef EXT_VERSION_MIINT
	data->versions.emplace_back("miint", EXT_VERSION_MIINT);
#else
	data->versions.emplace_back("miint", "unversioned");
#endif
	data->versions.emplace_back("htslib", hts_version());
	data->versions.emplace_back("minimap2", MINIMAP2_GIT_VERSION);
	data->versions.emplace_back("kseq++", KSEQPP_PROJECT_VERSION);
	data->versions.emplace_back("WFA2-lib", WFA2_GIT_VERSION);
#ifdef MIINT_HAS_HDF5
#ifdef H5_VERS_STR
	data->versions.emplace_back("HDF5", H5_VERS_STR);
#elif defined(H5_VERSION)
	data->versions.emplace_back("HDF5", H5_VERSION);
#else
	data->versions.emplace_back("HDF5", std::to_string(H5_VERS_MAJOR) + "." + std::to_string(H5_VERS_MINOR) + "." +
	                                        std::to_string(H5_VERS_RELEASE));
#endif
#endif
	data->versions.emplace_back("zlib", zlibVersion());
	data->versions.emplace_back("rype", RYPE_GIT_VERSION);
	return data;
}

static void MiintVersionsExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->CastNoConst<MiintVersionsData>();
	if (data.done) {
		output.SetCardinality(0);
		return;
	}
	idx_t count = data.versions.size();
	for (idx_t i = 0; i < count; i++) {
		FlatVector::GetData<string_t>(output.data[0])[i] =
		    StringVector::AddString(output.data[0], data.versions[i].first);
		FlatVector::GetData<string_t>(output.data[1])[i] =
		    StringVector::AddString(output.data[1], data.versions[i].second);
	}
	output.SetCardinality(count);
	data.done = true;
}

static void LoadInternal(ExtensionLoader &loader) {
	// TODO: use [[nodiscard]] throughout in headers
	// TODO: //! comment on headers

	ScalarFunction version_func("miint_version", {}, LogicalType::VARCHAR, MiintVersionFunction);
	loader.RegisterFunction(version_func);

	TableFunction versions_table_func("miint_versions", {}, MiintVersionsExecute, MiintVersionsBind);
	loader.RegisterFunction(versions_table_func);

	ReadFastxTableFunction::Register(loader);
	ReadAlignmentsTableFunction::Register(loader);
	ReadSequencesSamTableFunction::Register(loader);
	ReadSequencesSFFTableFunction::Register(loader);
#ifdef MIINT_HAS_HDF5
	ReadBIOMTableFunction::Register(loader);
#endif
	ReadNewickTableFunction::Register(loader);
	ReadJplaceNewickTableFunction::Register(loader);
	TreeResolvePlacementTableFunction::Register(loader);
	AlignMinimap2TableFunction::Register(loader);
	AlignMinimap2ShardedTableFunction::Register(loader);
	SaveMinimap2IndexTableFunction::Register(loader);
#ifdef MIINT_HAS_BOWTIE2
	AlignBowtie2TableFunction::Register(loader);
	AlignBowtie2ShardedTableFunction::Register(loader);
	RegisterBowtie2AvailableFunction(loader);
#else
	// Stub: bowtie2_available() always returns false when Bowtie2 support is compiled out
	ScalarFunction bowtie2_stub(
	    "bowtie2_available", {}, LogicalType::BOOLEAN,
	    [](DataChunk &args, ExpressionState &state, Vector &result) { result.Reference(Value::BOOLEAN(false)); });
	loader.RegisterFunction(bowtie2_stub);
#endif
	ReadNCBIFastaTableFunction::Register(loader);
	ReadNCBITableFunction::Register(loader);
	ReadNCBIAnnotationTableFunction::Register(loader);
	ReadENATableFunction::Register(loader);
	ReadENAAttributesTableFunction::Register(loader);
	ReadENASequencesTableFunction::Register(loader);

	AlignmentFlagFunctions::Register(loader);
	AlignmentSeqIdentityFunction::Register(loader);
	AlignmentQueryLengthFunction::Register(loader);
	AlignmentQueryCoverageFunction::Register(loader);
	CompressIntervalsFunction::Register(loader);
	ComputeCoverageDepthFunction::Register(loader);
	AlignmentSliceTableFunction::Register(loader);
	SequenceFunctions::Register(loader);
	FormulaFunction::Register(loader);
	MassQLFunction::Register(loader);
	MzmlPeakPairFunction::Register(loader);

	AlignPairwiseScoreFunction::Register(loader);
	AlignPairwiseCigarFunction::Register(loader);
	AlignPairwiseFullFunction::Register(loader);

	AlignMafftTableFunction::Register(loader);
	DeblurTableFunction::Register(loader);

#ifdef MIINT_HAS_HDF5
	CopyBiomFunction::Register(loader);
#endif
	CopyFastqFunction::Register(loader);
	CopyFastaFunction::Register(loader);
	CopyNewickFunction::Register(loader);
	CopySAMFunction::Register(loader);

	ReadMzMLTableFunction::Register(loader);
	ReadMzMLChromatogramsTableFunction::Register(loader);
	ReadMzXMLTableFunction::Register(loader);

	// Ensure dependency extensions are loaded before registering macros that
	// reference their functions (e.g., read_jplace needs json's read_json,
	// parse_gff_attributes needs core_functions' map_from_entries).
	auto &instance = loader.GetDatabaseInstance();
	for (const auto *dep : {"json", "core_functions"}) {
		if (!instance.ExtensionIsLoaded(dep)) {
			ExtensionHelper::TryAutoLoadExtension(instance, dep);
		}
#ifndef DUCKDB_BUILD_LOADABLE_EXTENSION
		// LoadExtension links against symbols not available in loadable extension builds
		if (!instance.ExtensionIsLoaded(dep)) {
			DuckDB db_wrapper(instance);
			ExtensionHelper::LoadExtension(db_wrapper, dep);
		}
#endif
	}

	MIINTMacros::Register(loader);

	RypeClassifyTableFunction::Register(loader);
	RypeExtractMinimizerSetTableFunction::Register(loader);
	RypeExtractStrandMinimizersTableFunction::Register(loader);
	RypeLogRatioTableFunction::Register(loader);
#ifdef MIINT_HAS_VSEARCH
	UchimeRefTableFunction::Register(loader);
	UchimeDenovoTableFunction::Register(loader);
	MaskDustFunction::Register(loader);
	MergePairsFunction::Register(loader);
	SearchSequencesTableFunction::Register(loader);
	ClusterSequencesTableFunction::Register(loader);
#endif
}

void MiintExtension::Load(ExtensionLoader &loader) {
	SetDependencyLogging();
	SetupSignalHandling();
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
	duckdb::SetDependencyLogging();
	duckdb::SetupSignalHandling();
	duckdb::LoadInternal(loader);
}
}
