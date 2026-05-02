#include "gpl_boundary/arrow_ipc.hpp"

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "nanoarrow/nanoarrow.h"
#include "nanoarrow/nanoarrow_ipc.h"

namespace duckdb {
namespace miint {
namespace gpl_boundary {

namespace {

[[noreturn]] void throw_with_error(const char *prefix, const ArrowError &err, ArrowErrorCode code) {
	std::string msg = "gpl_boundary: ";
	msg += prefix;
	msg += " failed (code=" + std::to_string(code);
	if (err.message[0] != '\0') {
		msg += ", error=";
		msg += err.message;
	}
	msg += ")";
	throw std::runtime_error(msg);
}

} // namespace

// =============================================================================
// EncodeIpcStream
// =============================================================================

std::vector<uint8_t> EncodeIpcStream(ArrowSchema *schema, ArrowArray *arrays, std::size_t n_arrays) {
	if (!schema || !arrays || n_arrays == 0) {
		throw std::runtime_error("gpl_boundary: EncodeIpcStream requires schema + at least one array");
	}

	// Output buffer: this is what we'll return.
	ArrowBuffer output;
	ArrowBufferInit(&output);

	// Wire `output` into an ArrowIpcOutputStream.
	ArrowIpcOutputStream stream;
	std::memset(&stream, 0, sizeof(stream));
	ArrowErrorCode rc = ArrowIpcOutputStreamInitBuffer(&stream, &output);
	if (rc != NANOARROW_OK) {
		ArrowBufferReset(&output);
		throw std::runtime_error("gpl_boundary: ArrowIpcOutputStreamInitBuffer failed");
	}

	// Init the writer.
	ArrowIpcWriter writer;
	std::memset(&writer, 0, sizeof(writer));
	rc = ArrowIpcWriterInit(&writer, &stream);
	if (rc != NANOARROW_OK) {
		// `stream` does NOT own `output` (ArrowIpcOutputStreamInitBuffer holds
		// a pointer to `output`, never copies or owns the bytes). On WriterInit
		// failure the stream may have been moved into the writer's private
		// data, so the writer's reset would handle it; on the path before that
		// move, the stream's release callback is responsible for its private
		// wrapper. Cover both by releasing the stream and explicitly resetting
		// the output buffer.
		if (stream.release) {
			stream.release(&stream);
		}
		ArrowBufferReset(&output);
		throw std::runtime_error("gpl_boundary: ArrowIpcWriterInit failed");
	}

	// At this point the writer holds the stream; the stream points to `output`
	// (no ownership transfer of bytes — `output.data` is still allocated by
	// the buffer's allocator). We touch `output` only after the writer is reset.

	ArrowError err;
	ArrowErrorInit(&err);

	// Schema first.
	rc = ArrowIpcWriterWriteSchema(&writer, schema, &err);
	if (rc != NANOARROW_OK) {
		ArrowIpcWriterReset(&writer);
		// writer.Reset releases the underlying stream which releases `output`.
		throw_with_error("ArrowIpcWriterWriteSchema", err, rc);
	}

	// Each batch.
	for (std::size_t i = 0; i < n_arrays; ++i) {
		ArrowArrayView view;
		std::memset(&view, 0, sizeof(view));
		rc = ArrowArrayViewInitFromSchema(&view, schema, &err);
		if (rc != NANOARROW_OK) {
			// ArrowArrayViewInitFromSchema may have allocated child views
			// before failing; reset to free them.
			ArrowArrayViewReset(&view);
			ArrowIpcWriterReset(&writer);
			throw_with_error("ArrowArrayViewInitFromSchema", err, rc);
		}
		rc = ArrowArrayViewSetArray(&view, &arrays[i], &err);
		if (rc != NANOARROW_OK) {
			ArrowArrayViewReset(&view);
			ArrowIpcWriterReset(&writer);
			throw_with_error("ArrowArrayViewSetArray", err, rc);
		}
		rc = ArrowIpcWriterWriteArrayView(&writer, &view, &err);
		ArrowArrayViewReset(&view);
		if (rc != NANOARROW_OK) {
			ArrowIpcWriterReset(&writer);
			throw_with_error("ArrowIpcWriterWriteArrayView", err, rc);
		}
	}

	// End-of-stream sentinel: write a NULL array view, which the writer
	// translates into the canonical zero-length EOS marker that gpl-boundary's
	// arrow-ipc Rust crate requires when `StreamDecoder::with_require_alignment`
	// is true.
	rc = ArrowIpcWriterWriteArrayView(&writer, nullptr, &err);
	if (rc != NANOARROW_OK) {
		ArrowIpcWriterReset(&writer);
		throw_with_error("ArrowIpcWriterWriteArrayView(EOS)", err, rc);
	}

	// Flush. Resetting the writer flushes pending writes and tears down the
	// nested stream, leaving `output` populated with the encoded bytes.
	ArrowIpcWriterReset(&writer);

	// Move bytes into a std::vector. The output buffer's allocator is
	// nanoarrow's default; we can't transfer ownership, so we copy. The
	// expected sizes here are small (KB range for an MSA, MB range for
	// typical alignments) so a single copy is fine.
	std::vector<uint8_t> result(output.data, output.data + output.size_bytes);
	ArrowBufferReset(&output);
	return result;
}

// =============================================================================
// IpcStreamDecoder
// =============================================================================

struct IpcStreamDecoder::Impl {
	// `ArrowIpcInputStreamInitBuffer` does an `ArrowBufferMove` of the struct
	// fields (data pointer, size, allocator) — it does NOT memcpy the
	// underlying bytes. The decoder dereferences `bytes` lazily during
	// get_schema/get_next, so the caller-owned buffer MUST outlive us.
	// (Public API contract; do not relax without a deep audit of the
	// nanoarrow_ipc decoder's read path.)
	const void *bytes = nullptr;
	std::size_t len = 0;

	// nanoarrow's ArrowArrayStream API; we drive it manually after init.
	ArrowArrayStream stream;
	bool stream_inited = false;

	// Cached schema (released by the decoder on destruction unless ownership
	// has already been transferred to the caller via GetSchema).
	bool schema_taken = false;

	Impl() {
		std::memset(&stream, 0, sizeof(stream));
	}
	~Impl() {
		if (stream_inited && stream.release) {
			stream.release(&stream);
		}
	}
};

IpcStreamDecoder::IpcStreamDecoder(const void *bytes, std::size_t len) : impl_(std::make_unique<Impl>()) {
	if (!bytes || len == 0) {
		throw std::runtime_error("gpl_boundary: IpcStreamDecoder needs non-empty input");
	}
	impl_->bytes = bytes;
	impl_->len = len;

	// Build an ArrowBuffer that aliases (does NOT own) the caller's bytes.
	// nanoarrow's init-from-buffer takes ownership of the buffer's
	// allocator-managed bytes; using a no-op deallocator keeps us safe.
	ArrowBuffer aliased_buffer;
	ArrowBufferInit(&aliased_buffer);
	// `ArrowBufferInit` sets a default allocator. Override to a no-op so the
	// "free" path doesn't call free() on caller-owned bytes.
	aliased_buffer.allocator =
	    ArrowBufferDeallocator([](ArrowBufferAllocator *, uint8_t *, int64_t) { /* caller owns */ }, nullptr);
	aliased_buffer.data = const_cast<uint8_t *>(static_cast<const uint8_t *>(bytes));
	aliased_buffer.size_bytes = static_cast<int64_t>(len);
	aliased_buffer.capacity_bytes = static_cast<int64_t>(len);

	ArrowIpcInputStream input_stream;
	std::memset(&input_stream, 0, sizeof(input_stream));
	ArrowErrorCode rc = ArrowIpcInputStreamInitBuffer(&input_stream, &aliased_buffer);
	if (rc != NANOARROW_OK) {
		// init failed: we still own aliased_buffer. Releasing is a no-op
		// because of the deallocator override.
		throw std::runtime_error("gpl_boundary: ArrowIpcInputStreamInitBuffer failed");
	}

	ArrowIpcArrayStreamReaderOptions opts {};
	opts.field_index = -1;
	opts.use_shared_buffers = 0;

	rc = ArrowIpcArrayStreamReaderInit(&impl_->stream, &input_stream, &opts);
	if (rc != NANOARROW_OK) {
		// On reader-init failure, ownership of `input_stream` was NOT
		// transferred to the reader (the move only happens after a successful
		// private-data allocation). Release explicitly to avoid leaking the
		// input-stream wrapper struct.
		if (input_stream.release) {
			input_stream.release(&input_stream);
		}
		throw std::runtime_error("gpl_boundary: ArrowIpcArrayStreamReaderInit failed");
	}
	impl_->stream_inited = true;
}

IpcStreamDecoder::~IpcStreamDecoder() = default;
IpcStreamDecoder::IpcStreamDecoder(IpcStreamDecoder &&) noexcept = default;
IpcStreamDecoder &IpcStreamDecoder::operator=(IpcStreamDecoder &&) noexcept = default;

void IpcStreamDecoder::GetSchema(ArrowSchema *out) {
	if (!impl_ || !impl_->stream_inited) {
		throw std::runtime_error("gpl_boundary: IpcStreamDecoder not initialized");
	}
	std::memset(out, 0, sizeof(ArrowSchema));
	const int rc = impl_->stream.get_schema(&impl_->stream, out);
	if (rc != 0) {
		const char *err = impl_->stream.get_last_error(&impl_->stream);
		throw std::runtime_error(std::string("gpl_boundary: get_schema failed: ") + (err ? err : "unknown"));
	}
}

bool IpcStreamDecoder::NextBatch(ArrowArray *out) {
	if (!impl_ || !impl_->stream_inited) {
		throw std::runtime_error("gpl_boundary: IpcStreamDecoder not initialized");
	}
	std::memset(out, 0, sizeof(ArrowArray));
	const int rc = impl_->stream.get_next(&impl_->stream, out);
	if (rc != 0) {
		const char *err = impl_->stream.get_last_error(&impl_->stream);
		throw std::runtime_error(std::string("gpl_boundary: get_next failed: ") + (err ? err : "unknown"));
	}
	if (!out->release) {
		// nanoarrow's convention: a released-array (release=null) signals EOS.
		return false;
	}
	return true;
}

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
