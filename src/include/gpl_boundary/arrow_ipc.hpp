#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

// Forward-declare from nanoarrow's C ABI. We re-export the standard structs.
struct ArrowSchema;
struct ArrowArray;

namespace duckdb {
namespace miint {
namespace gpl_boundary {

/// Encode an Arrow IPC stream from one schema and one or more record batches.
///
/// Inputs:
/// - `schema` — populated `ArrowSchema` (struct of named columns). The encoder
///   READS the schema — it does NOT call `release()` on it. Caller retains
///   full ownership and must release it themselves on every code path.
/// - `arrays` — array of `n_arrays` populated `ArrowArray`s, each conforming
///   to `schema`. Same ownership rules: the encoder reads the data, never
///   calls `release()`. Caller must release each array.
/// - `n_arrays` — count of batches in `arrays`. Must be ≥ 1.
///
/// Returns a buffer of IPC stream bytes — Schema message followed by one
/// RecordBatch message per array, each encapsulated with the IPC continuation
/// marker, then a zero-length End-Of-Stream marker. This is the format
/// gpl-boundary's `arrow-ipc` Rust crate consumes via `StreamDecoder`.
///
/// Throws `std::runtime_error` on internal allocation / FlatBuffers / I/O
/// errors. On exception, the schema and arrays are still owned by the caller
/// and must be released as usual.
std::vector<uint8_t> EncodeIpcStream(ArrowSchema *schema, ArrowArray *arrays, std::size_t n_arrays);

/// Decode an Arrow IPC stream from a contiguous byte buffer.
///
/// The buffer must remain valid for the lifetime of the decoder (the decoder
/// reads bytes lazily via an in-memory cursor). For the gpl-boundary use case
/// this is the lifetime of the `OutputShmRegion` that owns the mmap'd shm.
///
/// Usage:
///   IpcStreamDecoder decoder(bytes, len);
///   ArrowSchema schema {};
///   decoder.GetSchema(&schema);
///   ArrowArray batch {};
///   while (decoder.NextBatch(&batch)) { ... process ... batch.release(...); }
///   schema.release(&schema);
class IpcStreamDecoder {
public:
	/// `bytes` and `len` must remain valid for the decoder's lifetime.
	IpcStreamDecoder(const void *bytes, std::size_t len);
	~IpcStreamDecoder();

	IpcStreamDecoder(const IpcStreamDecoder &) = delete;
	IpcStreamDecoder &operator=(const IpcStreamDecoder &) = delete;
	IpcStreamDecoder(IpcStreamDecoder &&) noexcept;
	IpcStreamDecoder &operator=(IpcStreamDecoder &&) noexcept;

	/// Read the schema from the stream. Must be the first decode call.
	/// Throws `std::runtime_error` on malformed input.
	///
	/// **Lifetime contract:** `*out` must be a zero-initialized or
	/// already-released `ArrowSchema`. We zero `*out` before reading, which
	/// would silently leak any contents the caller had populated earlier.
	void GetSchema(ArrowSchema *out);

	/// Read the next record batch from the stream. Returns false on
	/// end-of-stream (caller's `out` is not modified). Throws on malformed
	/// input.
	///
	/// **Lifetime contract:** `*out` must be a zero-initialized or
	/// already-released `ArrowArray`. We zero `*out` before reading, which
	/// would silently leak any contents the caller had populated earlier.
	bool NextBatch(ArrowArray *out);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
