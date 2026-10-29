cmake_minimum_required(VERSION 3.15)
include(FetchContent)

set(TARGET_NAME miint)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# from stobealign, seem needed for kseq++ note: bumped from c++17 to 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# How this build is laid out, and how to add a source or a dependency:
# docs/internals/build-system.md. Each embedded dependency lives in its own
# cmake/deps/<name>.cmake (sylph and st3 are part of rust_glue.cmake).

# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------
option(MIINT_ENABLE_CURL "Build with libcurl streaming-upload support" ON)
option(MIINT_ENABLE_HDF5 "Build with HDF5/BIOM support" ON)
option(MIINT_ENABLE_MAFFT
       "Build with MAFFT multiple sequence alignment support" ON)
option(MIINT_ENABLE_ABPOA "Build with abPOA partial-order alignment support" ON)
option(MIINT_ENABLE_SORTMERNA "Build with SortMeRNA rRNA alignment support" ON)
option(MIINT_ENABLE_SYLPH
       "Build with sylph FracMinHash relative-abundance profiling" ON)
option(MIINT_ENABLE_ST3
       "Build with SourceTracker3 (st3) microbial source attribution" ON)
option(
  MIINT_ENABLE_GPL_BOUNDARY
  "Build with gpl-boundary subsystem (process pipes + Arrow IPC over POSIX shm)"
  ON)
option(MIINT_ENABLE_UNIFRAC
       "Build with UniFrac (PCoA / PERMANOVA / Faith PD) support" ON)
option(MIINT_ENABLE_KREPP "Build with krepp phylogenetic placement support" ON)
option(MIINT_ENABLE_VSEARCH
       "Build with vsearch search, clustering and chimera detection" ON)
option(MIINT_ENABLE_LIBDEFLATE "Build htslib bgzf with libdeflate (Linux only)"
       ON)

# Platform matrix: the one place a platform turns a feature off. Each line says
# why in a few words; the full reason is in that feature's cmake/deps file.
if(EMSCRIPTEN)
  # curl: no network stack, and unresolved side-module imports trap at call
  # time. HDF5: its C++ static members become unresolvable GOT.mem imports.
  # MAFFT: 2-4 MB stack frames overrun duckdb-wasm's 1 MB stack. The rest rely
  # on POSIX processes, files or threads the browser lacks.
  set(MIINT_ENABLE_CURL OFF)
  set(MIINT_ENABLE_HDF5 OFF)
  set(MIINT_ENABLE_SORTMERNA OFF)
  set(MIINT_ENABLE_SYLPH OFF)
  set(MIINT_ENABLE_GPL_BOUNDARY OFF)
  set(MIINT_ENABLE_KREPP OFF)
  set(MIINT_ENABLE_MAFFT OFF)
  message(
    STATUS
      "Emscripten detected: disabling libcurl, HDF5, SortMeRNA, sylph, gpl-boundary, krepp, and MAFFT (WASM-incompatible)"
  )
elseif(WIN32)
  # All of these rely on POSIX APIs (mkdtemp, process pipes, shm, autotools);
  # MAFFT segfaults on MinGW. SortMeRNA's cmph needs <sys/time.h>, which MinGW
  # lacks (upstream CI tests MSVC only).
  set(MIINT_ENABLE_VSEARCH OFF)
  set(MIINT_ENABLE_MAFFT OFF)
  set(MIINT_ENABLE_ABPOA OFF)
  set(MIINT_ENABLE_SYLPH OFF)
  set(MIINT_ENABLE_GPL_BOUNDARY OFF)
  set(MIINT_ENABLE_UNIFRAC OFF)
  set(MIINT_ENABLE_KREPP OFF)
  if(MINGW)
    set(MIINT_ENABLE_SORTMERNA OFF)
    message(
      STATUS
        "MinGW detected: disabling vsearch, MAFFT, SortMeRNA, sylph, gpl-boundary, UniFrac, and krepp (cmph + gpl-boundary rely on POSIX APIs; libssu inmem assumes POSIX)"
    )
  else()
    message(
      STATUS
        "Windows detected: disabling vsearch, MAFFT, sylph, gpl-boundary, UniFrac, and krepp (require POSIX APIs)"
    )
  endif()
elseif(APPLE)
  # vcpkg's libcurl brings OpenSSL's MD5_*/SHA1_*, which clash with vsearch's,
  # and Apple's ld has no --allow-multiple-definition.
  set(MIINT_ENABLE_CURL OFF)
  message(
    STATUS
      "Apple detected: disabling libcurl streaming transport (vsearch+OpenSSL symbol clash)"
  )
endif()
# libdeflate: Windows' htslib is a hand-written config (no autotools probe) and
# wasm_eh is single-threaded in-browser, so neither benefits.
if(MIINT_ENABLE_LIBDEFLATE AND NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(
    STATUS "MIINT_ENABLE_LIBDEFLATE only supported on Linux; disabling here")
  set(MIINT_ENABLE_LIBDEFLATE OFF)
endif()

# ---------------------------------------------------------------------------
# Toolchain, system libraries and embedded dependencies
# ---------------------------------------------------------------------------
# Each file below appends to the MIINT_* registration lists that
# miint_configure_target() applies. Include order is link order: a library must
# come after every library that uses it. Start every list empty: this file runs
# inside DuckDB's tree, so a same-named variable from there or the cache would
# otherwise reach all three targets.
foreach(
  list
  EXTENSION_SOURCES
  TEST_SOURCES
  DEFINITIONS
  INCLUDE_DIRS
  LINK_LIBRARIES
  OBJECT_SOURCES
  TARGET_DEPENDENCIES
  WASM_LIBS
  SUBMODULE_VERSIONS)
  set(MIINT_${list} "")
endforeach()
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/miint_helpers.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/system_deps.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/openmp.cmake)

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})

include(ExternalProject)

foreach(dep htslib minimap2 wfa2 rust_glue)
  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/deps/${dep}.cmake)
endforeach()
foreach(
  feature
  MAFFT
  ABPOA
  KREPP
  VSEARCH
  SORTMERNA
  CURL
  UNIFRAC
  GPL_BOUNDARY
  HDF5)
  if(MIINT_ENABLE_${feature})
    string(TOLOWER ${feature} dep)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/deps/${dep}.cmake)
  endif()
endforeach()

# Libraries that the archives above depend on go last.
list(APPEND MIINT_LINK_LIBRARIES ZLIB::ZLIB expat::expat Threads::Threads)
list(APPEND MIINT_WASM_LIBS "$<TARGET_FILE:ZLIB::ZLIB>"
     "$<TARGET_FILE:expat::expat>")
if(HAVE_LIBZSTD)
  list(APPEND MIINT_LINK_LIBRARIES ${ZSTD_TARGET})
  list(APPEND MIINT_WASM_LIBS "$<TARGET_FILE:${ZSTD_TARGET}>")
endif()
if(MIINT_ENABLE_LIBDEFLATE)
  list(APPEND MIINT_LINK_LIBRARIES ${LIBDEFLATE_TARGET})
endif()
if(APPLE)
  # rype's chrono (via iana_time_zone) calls CoreFoundation.
  list(APPEND MIINT_LINK_LIBRARIES "-framework CoreFoundation")
elseif(WIN32)
  # HTSlib: -lregex for POSIX regex (regcomp/regexec/regfree). Rust
  # x86_64-pc-windows-gnu stdlib system deps (from rustc_target spec and std
  # #[link] attributes): ntdll, userenv, bcrypt, advapi32, ws2_32. MinGW driver
  # auto-adds: kernel32, user32, msvcrt, mingwex, mingw32, gcc.
  list(
    APPEND
    MIINT_LINK_LIBRARIES
    regex
    ws2_32
    ntdll
    userenv
    bcrypt
    advapi32)
elseif(UNIX)
  list(APPEND MIINT_LINK_LIBRARIES dl m)
endif()

list(JOIN MIINT_SUBMODULE_VERSIONS " " MIINT_SUBMODULE_VERSIONS_STR)
message(STATUS "Submodule versions: ${MIINT_SUBMODULE_VERSIONS_STR}")

# Headers of every dependency above, for every target in this directory
# (include_directories also reaches targets already created here, such as
# krepp). Subdirectories added above (abpoa, nanoarrow_ipc) do not inherit them.
include_directories(src/include third_party/microtar ext ${MIINT_INCLUDE_DIRS})

# ---------------------------------------------------------------------------
# Sources always compiled, sorted. A feature's own sources are in its cmake/deps
# file. unifrac_function_common.cpp and the beta-diversity / ordination wrappers
# that use it are UniFrac-independent despite the name: nothing here links
# libssu or libskbb.
# ---------------------------------------------------------------------------
list(
  APPEND
  MIINT_EXTENSION_SOURCES
  src/absquant.cpp
  src/absquant_cell_counts_function.cpp
  src/absquant_function.cpp
  src/absquant_orf_copies_function.cpp
  src/absquant_readers.cpp
  src/align_minimap2.cpp
  src/align_minimap2_sharded.cpp
  src/align_pairwise_ksw2_dual_affine_functions.cpp
  src/align_pairwise_ksw2_functions.cpp
  src/align_pairwise_ksw2_splice_functions.cpp
  src/align_pairwise_wfa2_functions.cpp
  src/alignment_flag_functions.cpp
  src/alignment_functions.cpp
  src/alignment_slice.cpp
  src/AlignmentSlicer.cpp
  src/aspera_send.cpp
  src/aspera_stream.cpp
  src/aspera_utils.cpp
  src/blast_client.cpp
  src/blast_parser.cpp
  src/blast_search.cpp
  src/catalog_utils.cpp
  src/cluster_kmeans.cpp
  src/cluster_kmeans_function.cpp
  src/cluster_upgma.cpp
  src/cluster_upgma_function.cpp
  src/community_distances.cpp
  src/community_distances_function.cpp
  src/compress_intervals.cpp
  src/compute_coverage_depth.cpp
  src/compute_msa_consensus.cpp
  src/compute_pileup.cpp
  src/copy_fasta.cpp
  src/copy_fastq.cpp
  src/copy_format_common.cpp
  src/copy_newick.cpp
  src/copy_sam.cpp
  src/copy_ubam.cpp
  src/CoverageDepthCalculator.cpp
  src/cumulative_coverage.cpp
  src/CumulativeCoverageAccumulator.cpp
  src/deblur.cpp
  src/deblur_table_function.cpp
  src/duckdb_seq_stream.cpp
  src/ena_alias_check.cpp
  src/ena_checklist.cpp
  src/ena_client.cpp
  src/ena_envelope_builder.cpp
  src/ena_experiments_insert.cpp
  src/ena_experiments_insert_op.cpp
  src/ena_insert_common.cpp
  src/ena_lifecycle.cpp
  src/ena_lifecycle_delete.cpp
  src/ena_lifecycle_functions.cpp
  src/ena_lifecycle_submit.cpp
  src/ena_modify_functions.cpp
  src/ena_parser.cpp
  src/ena_projects_insert.cpp
  src/ena_projects_insert_op.cpp
  src/ena_receipt_parser.cpp
  src/ena_reports_client.cpp
  src/ena_resolver_cache.cpp
  src/ena_run_info_extractor.cpp
  src/ena_runs_insert.cpp
  src/ena_runs_insert_op.cpp
  src/ena_samples_insert.cpp
  src/ena_samples_insert_op.cpp
  src/ena_secret.cpp
  src/ena_storage.cpp
  src/ena_submission_log.cpp
  src/ena_upload_helpers.cpp
  src/ena_upload_reads.cpp
  src/ensure_httpfs.cpp
  src/extract_linked_amplicon.cpp
  src/fastq_encoder.cpp
  src/formula_function.cpp
  src/formula_parser.cpp
  src/gfa_parser.cpp
  src/hfile_duckdb.cpp
  src/id_column_codec.cpp
  src/IntervalCompressor.cpp
  src/ks_2samp.cpp
  src/KsTwoSample.cpp
  src/KSW2Aligner.cpp
  src/line_file_scan.cpp
  src/line_splitter.cpp
  src/mass_tables.cpp
  src/massql_function.cpp
  src/massql_parser.cpp
  src/massql_transpiler.cpp
  src/match_short_barcodes.cpp
  src/miint_cache.cpp
  src/miint_extension.cpp
  src/miint_isa.cpp
  src/miint_log.cpp
  src/miint_parallel.cpp
  src/minimap2_part_cursor.cpp
  src/Minimap2Aligner.cpp
  src/mmvec.cpp
  src/mmvec_derived_functions.cpp
  src/mmvec_fit_function.cpp
  src/mmvec_relation.cpp
  src/MsaConsensus.cpp
  src/MsBinaryDecoder.cpp
  src/mzml_peak_pair_function.cpp
  src/MzMLReader.cpp
  src/MzXMLReader.cpp
  src/ncbi_client.cpp
  src/ncbi_parser.cpp
  src/NewickTree.cpp
  src/paf_parser.cpp
  src/per_run_reader.cpp
  src/per_sample_table_function.cpp
  src/phylo_ancestral_ml.cpp
  src/phylo_ancestral_parsimony.cpp
  src/phylo_ancestral_states.cpp
  src/phylo_independent_contrasts.cpp
  src/phylo_traits_reader.cpp
  src/pick_anchors.cpp
  src/pick_anchors_function.cpp
  src/PileupWalker.cpp
  src/placement_table_reader.cpp
  src/procrustes_core.cpp
  src/procrustes_function.cpp
  src/qc_algorithms.cpp
  src/qc_functions.cpp
  src/QualScore.cpp
  src/read_alignment_header.cpp
  src/read_alignments.cpp
  src/read_ena.cpp
  src/read_ena_attributes.cpp
  src/read_ena_searchable_fields.cpp
  src/read_ena_sequences.cpp
  src/read_fastx.cpp
  src/read_gfa.cpp
  src/read_jplace_newick.cpp
  src/read_mzml.cpp
  src/read_mzml_chromatograms.cpp
  src/read_mzxml.cpp
  src/read_ncbi.cpp
  src/read_ncbi_annotation.cpp
  src/read_ncbi_fasta.cpp
  src/read_ncbi_lineage.cpp
  src/read_ncbi_taxdump.cpp
  src/read_newick.cpp
  src/read_paf.cpp
  src/read_sequences_sam.cpp
  src/read_sequences_sff.cpp
  src/reference_table_reader.cpp
  src/remote_file_helper.cpp
  src/rype_classify.cpp
  src/rype_extract.cpp
  src/rype_filter.cpp
  src/rype_index_create.cpp
  src/rype_input_stream.cpp
  src/rype_log_ratio.cpp
  src/sam_style_tags.cpp
  src/SAMReader.cpp
  src/save_minimap2_index.cpp
  src/sc_common.cpp
  src/sc_coo_builder.cpp
  src/sc_cross_validate_function.cpp
  src/sc_feature_importances_function.cpp
  src/sc_fit_function.cpp
  src/sc_model_features_function.cpp
  src/sc_predict_function.cpp
  src/sc_predict_proba_function.cpp
  src/sc_rf_common.cpp
  src/sc_shap_function.cpp
  src/sequence_data_reader.cpp
  src/sequence_functions.cpp
  src/sequence_table_reader.cpp
  src/SequenceReader.cpp
  src/SFFReader.cpp
  src/shear_tree.cpp
  src/simulate_resemblance.cpp
  src/simulate_resemblance_function.cpp
  src/table_function_common.cpp
  src/taxdump_archive.cpp
  src/taxdump_parser.cpp
  src/taxonomy_lineage.cpp
  src/tree_resolve_multifurcations.cpp
  src/tree_resolve_placement.cpp
  src/tree_table_reader.cpp
  src/unifrac_function_common.cpp
  src/WFA2Aligner.cpp
  src/woltka_ogu_function.cpp
  src/zip_utils.cpp
  third_party/microtar/microtar.c)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/isa_dispatch.cmake)
miint_add_isa_variants(mmvec_kernel src/mmvec_kernel.cpp)

# Compiled into `tests` without MIINT_STATIC_BUILD (the tests binary has no
# libduckdb), so a source listed here is compiled again rather than shared with
# the extension. Only DuckDB-free sources belong here. miniz.cpp is for the
# ExtractFromZip tests.
list(
  APPEND
  MIINT_TEST_SOURCES
  duckdb/third_party/miniz/miniz.cpp
  src/absquant.cpp
  src/AlignmentSlicer.cpp
  src/aspera_utils.cpp
  src/blast_parser.cpp
  src/cluster_kmeans.cpp
  src/cluster_upgma.cpp
  src/community_distances.cpp
  src/CoverageDepthCalculator.cpp
  src/CumulativeCoverageAccumulator.cpp
  src/deblur.cpp
  src/ena_alias_check.cpp
  src/ena_checklist.cpp
  src/ena_envelope_builder.cpp
  src/ena_experiments_insert.cpp
  src/ena_lifecycle_submit.cpp
  src/ena_parser.cpp
  src/ena_projects_insert.cpp
  src/ena_receipt_parser.cpp
  src/ena_reports_client.cpp
  src/ena_resolver_cache.cpp
  src/ena_run_info_extractor.cpp
  src/ena_runs_insert.cpp
  src/ena_samples_insert.cpp
  src/ena_upload_helpers.cpp
  src/fastq_encoder.cpp
  src/formula_parser.cpp
  src/gfa_parser.cpp
  src/id_column_codec.cpp
  src/IntervalCompressor.cpp
  src/KsTwoSample.cpp
  src/KSW2Aligner.cpp
  src/line_splitter.cpp
  src/mass_tables.cpp
  src/massql_parser.cpp
  src/massql_transpiler.cpp
  src/miint_isa.cpp
  src/miint_parallel.cpp
  src/minimap2_part_cursor.cpp
  src/Minimap2Aligner.cpp
  src/mmvec.cpp
  src/mmvec_relation.cpp
  src/MsaConsensus.cpp
  src/MsBinaryDecoder.cpp
  src/MzMLReader.cpp
  src/MzXMLReader.cpp
  src/ncbi_parser.cpp
  src/NewickTree.cpp
  src/paf_parser.cpp
  src/pick_anchors.cpp
  src/PileupWalker.cpp
  src/procrustes_core.cpp
  src/qc_algorithms.cpp
  src/QualScore.cpp
  src/remote_file_helper.cpp
  src/sam_style_tags.cpp
  src/SAMReader.cpp
  src/sc_coo_builder.cpp
  src/SequenceReader.cpp
  src/SFFReader.cpp
  src/simulate_resemblance.cpp
  src/taxdump_archive.cpp
  src/taxdump_parser.cpp
  src/taxonomy_lineage.cpp
  src/WFA2Aligner.cpp
  src/zip_utils.cpp
  test/cpp/test_AbsQuant.cpp
  test/cpp/test_AlignBowtie2Sharded.cpp
  test/cpp/test_AlignmentFunctions.cpp
  test/cpp/test_AlignmentSlicer.cpp
  test/cpp/test_AncestralStatesContinuous.cpp
  test/cpp/test_AncestralStatesML.cpp
  test/cpp/test_AncestralStatesParsimony.cpp
  test/cpp/test_AsperaUtils.cpp
  test/cpp/test_blast_parser.cpp
  test/cpp/test_Bowtie2Nthreads.cpp
  test/cpp/test_Bowtie2SubjectsBatching.cpp
  test/cpp/test_ClusterKmeans.cpp
  test/cpp/test_ClusterUpgma.cpp
  test/cpp/test_CommunityDistances.cpp
  test/cpp/test_compute_pileup.cpp
  test/cpp/test_consensus_msa.cpp
  test/cpp/test_copy_ubam.cpp
  test/cpp/test_CopyFileHandlePipe.cpp
  test/cpp/test_CoverageDepthCalculator.cpp
  test/cpp/test_CumulativeCoverageAccumulator.cpp
  test/cpp/test_deblur.cpp
  test/cpp/test_DuckDBSeqStream.cpp
  test/cpp/test_ena_alias_check.cpp
  test/cpp/test_ena_checklist.cpp
  test/cpp/test_ena_envelope.cpp
  test/cpp/test_ena_experiments_insert.cpp
  test/cpp/test_ena_lifecycle_submit.cpp
  test/cpp/test_ena_md5.cpp
  test/cpp/test_ena_projects_insert.cpp
  test/cpp/test_ena_receipt.cpp
  test/cpp/test_ena_reports_client.cpp
  test/cpp/test_ena_runs_insert.cpp
  test/cpp/test_ena_samples_insert.cpp
  test/cpp/test_ena_secret.cpp
  test/cpp/test_ena_upload_reads.cpp
  test/cpp/test_ENAClient.cpp
  test/cpp/test_ENAResolverCache.cpp
  test/cpp/test_fastq_encoder.cpp
  test/cpp/test_gfa_parser.cpp
  test/cpp/test_hfile_duckdb.cpp
  test/cpp/test_id_column_codec.cpp
  test/cpp/test_IndependentContrasts.cpp
  test/cpp/test_InsertFullyResolved.cpp
  test/cpp/test_IntervalCompressor.cpp
  test/cpp/test_KsTwoSample.cpp
  test/cpp/test_KSW2Aligner.cpp
  test/cpp/test_line_splitter.cpp
  test/cpp/test_MassQLParser.cpp
  test/cpp/test_MiintIsa.cpp
  test/cpp/test_MiintParallel.cpp
  test/cpp/test_Minimap2Aligner.cpp
  test/cpp/test_MMvec.cpp
  test/cpp/test_MMvecRelation.cpp
  test/cpp/test_MsBinaryDecoder.cpp
  test/cpp/test_MzMLReader.cpp
  test/cpp/test_MzXMLReader.cpp
  test/cpp/test_ncbi_client.cpp
  test/cpp/test_NewickParser.cpp
  test/cpp/test_paf_parser.cpp
  test/cpp/test_PickAnchors.cpp
  test/cpp/test_Procrustes.cpp
  test/cpp/test_qc_algorithms.cpp
  test/cpp/test_QualScore.cpp
  test/cpp/test_RemoteFileHelper.cpp
  test/cpp/test_ResolveMultifurcations.cpp
  test/cpp/test_sam_style_tags.cpp
  test/cpp/test_SAMReader.cpp
  test/cpp/test_ScCooBuilder.cpp
  test/cpp/test_SequenceReader.cpp
  test/cpp/test_SFFReader.cpp
  test/cpp/test_ShardProgress.cpp
  test/cpp/test_ShearTree.cpp
  test/cpp/test_sigpipe_bootstrap.cpp
  test/cpp/test_SimulateCluster.cpp
  test/cpp/test_SimulateGradient.cpp
  test/cpp/test_taxdump_archive.cpp
  test/cpp/test_taxdump_parser.cpp
  test/cpp/test_taxonomy_lineage.cpp
  test/cpp/test_WFA2Aligner.cpp
  third_party/microtar/microtar.c)

miint_check_unique(MIINT_EXTENSION_SOURCES)
miint_check_unique(MIINT_TEST_SOURCES)

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
build_static_extension(${TARGET_NAME} ${MIINT_EXTENSION_SOURCES})

# WASM builds: DuckDB's build_loadable_extension runs a POST_BUILD emcc command
# to create the .wasm side module. It reads additional libraries to link from
# DUCKDB_EXTENSION_<NAME>_LINKED_LIBS — NOT from target_link_libraries (which
# only records metadata for STATIC library targets on Emscripten). Without this,
# symbols like gzclose become unresolved imports and fail at extension load time
# with "bad export type: undefined". extension_build_tools.cmake runs
# separate_arguments() on this space-separated string, then passes the result to
# add_custom_command COMMAND, where generator expressions are evaluated at build
# time.
if(EMSCRIPTEN)
  list(JOIN MIINT_WASM_LIBS " " DUCKDB_EXTENSION_MIINT_LINKED_LIBS)
  message(STATUS "WASM linked libs: ${DUCKDB_EXTENSION_MIINT_LINKED_LIBS}")
endif()

build_loadable_extension(${TARGET_NAME} " " ${MIINT_EXTENSION_SOURCES})

foreach(target ${EXTENSION_NAME} ${LOADABLE_EXTENSION_NAME})
  miint_configure_target(${target})
  target_compile_definitions(${target} PRIVATE RYPE_ARROW MIINT_STATIC_BUILD)
endforeach()
# The loadable extension is dlopen'd at runtime. When the host DuckDB is loaded
# via Python (RTLD_LOCAL), duckdb_je_* symbols from the host are not visible.
# Link stubs that forward to system malloc so the .so is self-contained.
target_sources(
  ${LOADABLE_EXTENSION_NAME}
  PRIVATE $<${MIINT_JEMALLOC}:${CMAKE_CURRENT_SOURCE_DIR}/src/mm_alloc_stubs.c>)

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")

FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0)
FetchContent_MakeAvailable(catch2)

add_executable(tests ${MIINT_TEST_SOURCES})
target_link_libraries(tests Catch2::Catch2WithMain)
miint_configure_target(tests)
target_include_directories(
  tests
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/cpp/include # Catch2 helper headers
          ${CMAKE_CURRENT_BINARY_DIR}/htslib/include
          # For the ExtractFromZip tests
          ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/miniz)
target_sources(
  tests
  PRIVATE $<${MIINT_JEMALLOC}:${CMAKE_CURRENT_SOURCE_DIR}/src/mm_alloc_stubs.c>)

# vsearch ships its own MD5/SHA1 symbols (libvsearch_core_a-md5.o,
# libvsearch_core_a-sha1.o), and OpenSSL's libcrypto — pulled in transitively
# when we link CURL::libcurl — exports the same names. GNU ld's
# --allow-multiple-definition picks the first-defined symbol; since each archive
# uses its own copy at compile time (intra-.o calls stay within the archive),
# runtime behaviour is preserved. libcurl itself uses `EVP_Digest*` rather than
# the colliding `MD5_*` C API, so the choice is benign.
#
# This workaround is only needed when both vsearch AND libcurl are linked. Apple
# disables MIINT_ENABLE_CURL (see the platform matrix) precisely because Apple's
# ld lacks an equivalent flag and vsearch must remain first-class on macOS;
# Windows auto-disables vsearch. So we only need this on Linux.
#
# Windows MinGW needs it for a different reason: duckdb's LogicalType static
# members (VARCHAR, DOUBLE, BIGINT, ...) are emitted into COMDAT sections in
# every translation unit that references them. GNU ld on Linux folds these by
# default; rtools42's MinGW ld on Windows does not, so the duckdb.dll link sees
# duplicate `duckdb::LogicalType::VARCHAR` symbols between
# libmiint_extension.a's TUs and duckdb_common's TUs.
#
# The static `${EXTENSION_NAME}` archive is consumed by DuckDB-built binaries
# (duckdb shell, unittest, plan_serializer) that we don't control. INTERFACE
# makes the link option propagate to whichever executable links against this
# archive, so the duplicate-symbol tolerance follows the static lib instead of
# being a per-target chore.
if((MIINT_ENABLE_VSEARCH
    AND MIINT_ENABLE_CURL
    AND UNIX
    AND NOT APPLE
   )
   OR (WIN32 AND MINGW))
  target_link_options(tests PRIVATE "-Wl,--allow-multiple-definition")
  target_link_options(${EXTENSION_NAME} INTERFACE
                      "-Wl,--allow-multiple-definition")
  target_link_options(${LOADABLE_EXTENSION_NAME} PRIVATE
                      "-Wl,--allow-multiple-definition")
endif()

# C++ tests that need a running DuckDB: they open an in-process database and
# drive miint code against it. `tests` above links no DuckDB, so these get a
# binary of their own, linking DuckDB and miint statically, built only where
# DuckDB builds its unittest runner (not in the WASM builds, for one).
if(BUILD_UNITTESTS)
  add_executable(
    tests_duckdb
    test/cpp/test_RypeInputStream.cpp test/cpp/test_RypeFilterParallel.cpp
    test/cpp/test_RypeFilterSchedule.cpp)
  target_link_libraries(tests_duckdb Catch2::Catch2WithMain duckdb_static)
  link_extension_libraries(tests_duckdb "")
endif()
