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
#include <copy_ubam.hpp>
#include <read_jplace_newick.hpp>
#include <read_newick.hpp>
#include <shear_tree.hpp>
#include <tree_resolve_placement.hpp>
#include <kseq++/seqio.hpp>
#include <read_fastx.hpp>
#include <read_alignments.hpp>
#include <read_sequences_sam.hpp>
#include <read_sequences_sff.hpp>
#include <align_minimap2.hpp>
#include <align_minimap2_sharded.hpp>
#include <save_minimap2_index.hpp>
#include <save_bowtie2_index.hpp>
#include <read_ncbi_fasta.hpp>
#include <read_ncbi.hpp>
#include <read_ncbi_annotation.hpp>
#include <read_ncbi_taxdump.hpp>
#include <read_ncbi_lineage.hpp>
#include <blast_search.hpp>
#include <read_ena.hpp>
#include <read_ena_attributes.hpp>
#include <read_ena_searchable_fields.hpp>
#include <read_ena_sequences.hpp>
#ifdef MIINT_HAS_CURL
#include <curl_send.hpp>
#endif
#include <ena_lifecycle_functions.hpp>
#include <ena_modify_functions.hpp>
#include <ena_secret.hpp>
#include <ena_storage.hpp>
#include <ena_upload_reads.hpp>
#include <miint_log.hpp>
#include <miint_macros.hpp>
#include "duckdb/main/extension_helper.hpp"
#include "duckdb/main/database.hpp"
#include <sequence_functions.hpp>
#include <qc_functions.hpp>
#include <formula_function.hpp>
#include <massql_function.hpp>
#include <woltka_ogu_function.hpp>
#include <mzml_peak_pair_function.hpp>
#ifdef MIINT_HAS_MAFFT
#include <align_mafft.hpp>
#endif
#ifdef MIINT_HAS_ABPOA
#include <align_abpoa.hpp>
#include <consensus_abpoa.hpp>
#endif
#ifdef MIINT_HAS_SORTMERNA
#include <align_sortmerna.hpp>
#include <align_sortmerna_rrna.hpp>
#endif
#include <deblur_table_function.hpp>
#include <align_pairwise_wfa2_functions.hpp>
#include <align_pairwise_ksw2_functions.hpp>
#include <align_pairwise_ksw2_dual_affine_functions.hpp>
#include <align_pairwise_ksw2_splice_functions.hpp>
#include <extract_linked_amplicon.hpp>
#include <match_short_barcodes.hpp>
#include <compute_pileup.hpp>
#include <compute_msa_consensus.hpp>
#include <read_mzml.hpp>
#include <read_mzml_chromatograms.hpp>
#include <read_mzxml.hpp>
#include <rype_classify.hpp>
#include <rype_extract.hpp>
#include <rype_index_create.hpp>
#include <rype_log_ratio.hpp>
#ifdef MIINT_HAS_VSEARCH
#include <uchime_ref.hpp>
#include <uchime_denovo.hpp>
#include <mask_dust_function.hpp>
#include <merge_pairs_function.hpp>
#include <search_sequences.hpp>
#include <cluster_sequences.hpp>
#endif
#ifdef MIINT_HAS_UNIFRAC
#include <unifrac_table_functions.hpp>
#endif
#ifdef MIINT_HAS_SYLPH
#include <sylph_index_create.hpp>
#include <sylph_profile.hpp>
#endif
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>
#include <duckdb/main/config.hpp>
#include <duckdb/storage/storage_extension.hpp>
#include <htslib-1.22.1/htslib/hts.h>
#include <kseq++/config.hpp>
#include <zlib.h>
#ifdef MIINT_HAS_LIBDEFLATE
#include <libdeflate.h>
#endif

#ifdef MIINT_HAS_HDF5
#include <copy_biom.hpp>
#include <read_biom.hpp>
#include <hdf5.h>
#endif

#ifdef MIINT_HAS_GPL_BOUNDARY
#include <align_bowtie2.hpp>
#include <align_bowtie2_sharded.hpp>
#include <install_gpl_boundary.hpp>
#include <phylogeny_fasttree.hpp>
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
#if MIINT_ASPERA_SUPPORTED || defined(MIINT_HAS_GPL_BOUNDARY)
	// Ignore SIGPIPE globally so that writes to closed pipes return EPIPE instead of
	// killing the process. Needed for AsperaProcess and gpl-boundary's ChildProcess +
	// Session, where pipes may close unexpectedly.
	// Note: This is a PROCESS-WIDE setting that persists for the lifetime of the process.
	// Setting it once at extension load is thread-safe (vs calling signal() from multiple threads).
	// gpl_boundary::Session ALSO uses pthread_sigmask + sigtimedwait per-thread inside its
	// hot loop — that's per-thread defense; this is the global no-kill backstop.
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
#ifdef MIINT_HAS_LIBDEFLATE
	data->versions.emplace_back("libdeflate", LIBDEFLATE_VERSION_STRING);
#endif
#ifdef MIINT_HAS_CURL
	data->versions.emplace_back("libcurl", miint::GetCurlVersion());
#endif
	data->versions.emplace_back("rype", RYPE_GIT_VERSION);
#ifdef MIINT_HAS_VSEARCH
	data->versions.emplace_back("vsearch", VSEARCH_GIT_VERSION);
#endif
	data->versions.emplace_back("mafft", MAFFT_GIT_VERSION);
#ifdef MIINT_HAS_ABPOA
	data->versions.emplace_back("abpoa", ABPOA_GIT_VERSION);
#endif
#ifdef MIINT_HAS_SYLPH
	data->versions.emplace_back("sylph", SYLPH_GIT_VERSION);
#endif
#ifdef MIINT_HAS_UNIFRAC
	data->versions.emplace_back("unifrac", UNIFRAC_GIT_VERSION);
	data->versions.emplace_back("scikit-bio-binaries", SKBB_GIT_VERSION);
#endif
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

	miint::RegisterMiintLogType(loader.GetDatabaseInstance());
	miint::RegisterENASecretType(loader);

	auto &ena_db_config = DBConfig::GetConfig(loader.GetDatabaseInstance());
	StorageExtension::Register(ena_db_config, "ena", make_shared_ptr<ENAStorageExtension>());

	// Session-scoped flag that flips ena.* INSERT from ADD to VALIDATE
	// (server-side dry-run). When true, the next INSERT INTO ena.X builds a
	// VALIDATE envelope instead of ADD: no accessions are assigned, the
	// alias-collision check is skipped (the user may legitimately validate
	// against an alias they've already registered), and the
	// ena.submission_log row records action='VALIDATE' with empty
	// object_aliases / object_accessions. Default false → existing ADD path.
	ena_db_config.AddExtensionOption("miint_ena_validate_only",
	                                 "When true, the next INSERT INTO ena.* is sent as a Webin V2 VALIDATE "
	                                 "(server-side dry-run) instead of ADD. No accessions are returned. "
	                                 "Default false.",
	                                 LogicalType::BOOLEAN, Value::BOOLEAN(false));

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
	ShearTreeTableFunction::Register(loader);
	AlignMinimap2TableFunction::Register(loader);
	AlignMinimap2ShardedTableFunction::Register(loader);
	SaveMinimap2IndexTableFunction::Register(loader);
#ifdef MIINT_HAS_GPL_BOUNDARY
	AlignBowtie2TableFunction::Register(loader);
	AlignBowtie2ShardedTableFunction::Register(loader);
	SaveBowtie2IndexTableFunction::Register(loader);
	RegisterBowtie2AvailableFunction(loader);
#else
	// Stub: bowtie2_available() always returns false when the gpl-boundary
	// subsystem is compiled out (e.g., on WASM/Windows). The table functions
	// (align_bowtie2 / align_bowtie2_sharded) are also unavailable in that
	// case — the daemon is the only path.
	ScalarFunction bowtie2_stub(
	    "bowtie2_available", {}, LogicalType::BOOLEAN,
	    [](DataChunk &args, ExpressionState &state, Vector &result) { result.Reference(Value::BOOLEAN(false)); });
	loader.RegisterFunction(bowtie2_stub);
#endif
	ReadNCBIFastaTableFunction::Register(loader);
	ReadNCBITableFunction::Register(loader);
	ReadNCBIAnnotationTableFunction::Register(loader);
	ReadNCBITaxdumpTableFunction::Register(loader);
	ReadNCBITaxdumpMergedTableFunction::Register(loader);
	ReadNCBILineageTableFunction::Register(loader);
	BlastSearchTableFunction::Register(loader);
	ReadENATableFunction::Register(loader);
	ReadENAAttributesTableFunction::Register(loader);
	ReadENASearchableFieldsTableFunction::Register(loader);
	ReadENASequencesTableFunction::Register(loader);
	ENAUploadReadsTableFunction::Register(loader);
	miint::RegisterENALifecycleTableFunctions(loader);
	miint::RegisterENAModifyTableFunctions(loader);

	AlignmentFlagFunctions::Register(loader);
	AlignmentSeqIdentityFunction::Register(loader);
	CigarSequenceIdentityFunction::Register(loader);
	CigarQueryLengthFunction::Register(loader);
	CigarQueryCoverageFunction::Register(loader);
	CompressIntervalsFunction::Register(loader);
	ComputeCoverageDepthFunction::Register(loader);
	AlignmentSliceTableFunction::Register(loader);
	SequenceFunctions::Register(loader);
	QcFunctions::Register(loader);
	FormulaFunction::Register(loader);
	MassQLFunction::Register(loader);
	WoltkaOguFunction::Register(loader);
	MzmlPeakPairFunction::Register(loader);
#ifdef MIINT_HAS_UNIFRAC
	RegisterUnifracPcoa(loader);
	RegisterUnifracPermanova(loader);
	RegisterUnifracFaithPD(loader);
	RegisterUnifracDistances(loader);
#endif

	AlignPairwiseWfa2ScoreFunction::Register(loader);
	AlignPairwiseWfa2CigarFunction::Register(loader);
	AlignPairwiseWfa2FullFunction::Register(loader);

	AlignPairwiseKsw2ScoreFunction::Register(loader);
	AlignPairwiseKsw2CigarFunction::Register(loader);
	AlignPairwiseKsw2FullFunction::Register(loader);

	AlignPairwiseKsw2DualAffineScoreFunction::Register(loader);
	AlignPairwiseKsw2DualAffineCigarFunction::Register(loader);
	AlignPairwiseKsw2DualAffineFullFunction::Register(loader);

	AlignPairwiseKsw2SpliceScoreFunction::Register(loader);
	AlignPairwiseKsw2SpliceCigarFunction::Register(loader);
	AlignPairwiseKsw2SpliceFullFunction::Register(loader);

	ExtractLinkedAmpliconFunction::Register(loader);
	MatchShortBarcodesTableFunction::Register(loader);
	ComputePileupTableFunction::Register(loader);
	ComputeMsaConsensusFunction::Register(loader);

#ifdef MIINT_HAS_MAFFT
	AlignMafftTableFunction::Register(loader);
#endif
#ifdef MIINT_HAS_ABPOA
	AlignAbpoaTableFunction::Register(loader);
	ConsensusAbpoaTableFunction::Register(loader);
#endif
#ifdef MIINT_HAS_SORTMERNA
	AlignSortMeRNATableFunction::Register(loader);
	AlignSortMeRNARRNATableFunction::Register(loader);
#endif
#ifdef MIINT_HAS_SYLPH
	SylphProfileTableFunction::Register(loader);
	SylphIndexCreateTableFunction::Register(loader);
#endif
#ifdef MIINT_HAS_GPL_BOUNDARY
	PhylogenyFastTreeTableFunction::Register(loader);
	PhylogenyFastTreeAvailableScalar::Register(loader);
	InstallGplBoundaryScalar::Register(loader);
#else
	// Stub: phylogeny_fasttree_available() always returns false when gpl-boundary
	// support is compiled out (e.g., on WASM/Windows).
	ScalarFunction phylogeny_fasttree_stub(
	    "phylogeny_fasttree_available", {}, LogicalType::BOOLEAN,
	    [](DataChunk &args, ExpressionState &state, Vector &result) { result.Reference(Value::BOOLEAN(false)); });
	loader.RegisterFunction(phylogeny_fasttree_stub);

	// Stub: install_gpl_boundary() reports the platform doesn't support
	// gpl-boundary at all, so installing wouldn't help.
	const auto install_struct = LogicalType::STRUCT({{"installed", LogicalType::BOOLEAN},
	                                                 {"path", LogicalType::VARCHAR},
	                                                 {"version", LogicalType::VARCHAR},
	                                                 {"message", LogicalType::VARCHAR}});
	const auto install_stub_exec = [](DataChunk &args, ExpressionState &state, Vector &result) {
		auto &entries = StructVector::GetEntries(result);
		auto installed_data = FlatVector::GetData<bool>(*entries[0]);
		auto &path_vec = *entries[1];
		auto &version_vec = *entries[2];
		auto &message_vec = *entries[3];
		const idx_t n = args.size();
		const string msg = "install_gpl_boundary: this miint build was compiled without "
		                   "MIINT_ENABLE_GPL_BOUNDARY (typically WASM or Windows). gpl-boundary "
		                   "is not supported on this platform.";
		for (idx_t i = 0; i < n; i++) {
			installed_data[i] = false;
			FlatVector::GetData<string_t>(path_vec)[i] = StringVector::AddString(path_vec, "");
			FlatVector::GetData<string_t>(version_vec)[i] = StringVector::AddString(version_vec, "");
			FlatVector::GetData<string_t>(message_vec)[i] = StringVector::AddString(message_vec, msg);
		}
		result.SetVectorType(VectorType::CONSTANT_VECTOR);
	};
	// Mirror the real path's 0-arg + 1-arg(BOOLEAN) overload set so SQL that
	// passes `force` still type-checks on stub builds (returns the same
	// "unsupported platform" payload regardless of the flag).
	ScalarFunctionSet install_gpl_boundary_stub_set("install_gpl_boundary");
	install_gpl_boundary_stub_set.AddFunction(ScalarFunction({}, install_struct, install_stub_exec));
	install_gpl_boundary_stub_set.AddFunction(
	    ScalarFunction({LogicalType::BOOLEAN}, install_struct, install_stub_exec));
	loader.RegisterFunction(install_gpl_boundary_stub_set);
#endif
	DeblurTableFunction::Register(loader);

#ifdef MIINT_HAS_HDF5
	CopyBiomFunction::Register(loader);
#endif
	CopyFastqFunction::Register(loader);
	CopyFastaFunction::Register(loader);
	CopyNewickFunction::Register(loader);
	CopySAMFunction::Register(loader);
	CopyUBAMFunction::Register(loader);

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
	RypeIndexCreateTableFunction::Register(loader);
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
