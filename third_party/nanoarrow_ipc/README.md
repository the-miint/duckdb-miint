# Vendored: nanoarrow_ipc + nanoarrow + flatcc runtime

Vendored to give miint Arrow IPC encode/decode without adding a heavy vcpkg
dependency. Used by the `gpl-boundary` subsystem (`src/gpl_boundary/arrow_ipc.cpp`)
to exchange RecordBatches with the long-lived Rust daemon over POSIX shared
memory.

## Source

- Upstream: https://github.com/apache/arrow-nanoarrow
- Tag: see `VERSION`
- License: Apache 2.0 (see `LICENSE`, `NOTICE`)

## flatcc runtime

`nanoarrow_ipc` depends on flatcc (FlatBuffers C compiler runtime) for the
FlatBuffers schema/message envelopes that Arrow IPC uses. The flatcc runtime
is vendored under `src/flatcc_runtime/` and `include/flatcc/`. See
`LICENSE.flatcc` (Apache 2.0).

## Symbol namespacing

Compiled with `-DNANOARROW_NAMESPACE=miint`. All `Arrow*` and `ArrowIpc*`
symbols are prefixed with `miint_` to avoid collision with DuckDB's bundled
nanoarrow (which lives in the C++ namespace `duckdb_nanoarrow`).

## Updating

To pull a newer upstream tag:

1. Update `VERSION` with the new tag + commit hash.
2. Re-run `scripts/vendor-nanoarrow-ipc.sh` (TODO if needed) or manually
   refresh:
   ```
   git clone --depth 1 -b <new-tag> https://github.com/apache/arrow-nanoarrow.git /tmp/nanoarrow
   # copy: src/nanoarrow/{nanoarrow.h,nanoarrow_ipc.h,nanoarrow_config.h.in}
   #       src/nanoarrow/common/{inline_array.h,inline_buffer.h,inline_types.h,*.c}
   #       src/nanoarrow/ipc/{*.c,flatcc_generated.h}
   #       thirdparty/flatcc/include/flatcc/{*.h,portable/*.h}
   #       thirdparty/flatcc/src/runtime/*.c
   ```
3. Re-run all C++ unit tests; pay special attention to round-trip tests
   in `test/cpp/test_gpl_boundary_arrow_ipc.cpp`.
