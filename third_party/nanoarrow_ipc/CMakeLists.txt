# Vendored nanoarrow + nanoarrow_ipc + flatcc runtime.
# See README.md for upstream source and updating instructions.

set(NANOARROW_IPC_SOURCES
    src/common/array.c
    src/common/array_stream.c
    src/common/schema.c
    src/common/utils.c
    src/ipc/codecs.c
    src/ipc/decoder.c
    src/ipc/encoder.c
    src/ipc/reader.c
    src/ipc/writer.c
    src/flatcc_runtime/builder.c
    src/flatcc_runtime/emitter.c
    src/flatcc_runtime/refmap.c
    src/flatcc_runtime/verifier.c)

# OBJECT library: avoids the export-set headache when miint_extension is
# install()'d into DuckDB's export set. Consumers link the object files
# directly without producing a separate .a that would need to be exported.
add_library(miint_nanoarrow_ipc OBJECT ${NANOARROW_IPC_SOURCES})

target_include_directories(miint_nanoarrow_ipc
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Hide all our symbols by default; consumers reach in through the headers.
# Combined with NANOARROW_NAMESPACE=miint (set in nanoarrow_config.h), this
# prevents collision with DuckDB's bundled nanoarrow (duckdb_nanoarrow C++
# namespace) and any other Arrow library that may end up in the same binary.
set_target_properties(miint_nanoarrow_ipc PROPERTIES
    C_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON)

# Upstream nanoarrow + flatcc accept warnings we don't want to inherit
# in miint's stricter warning-as-error mode.
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(miint_nanoarrow_ipc PRIVATE
        -w  # vendored upstream code; warnings owned by them
    )
endif()
