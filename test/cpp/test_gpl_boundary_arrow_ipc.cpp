#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gpl_boundary/arrow_ipc.hpp"

// We pull in the vendored nanoarrow directly here ONLY in the test file to
// hand-build small ArrowSchema/ArrowArray fixtures. Production code
// (`src/gpl_boundary/arrow_ipc.cpp`) uses the same headers; this isn't a
// public dependency for callers.
#include "nanoarrow/nanoarrow.h"

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

using namespace duckdb::miint::gpl_boundary;

// =============================================================================
// Cycle 3.1 — Encode (Arrow C Data Interface → Arrow IPC stream bytes)
// =============================================================================
//
// We build a small schema (struct: name VARCHAR, count INT64), populate two
// rows, encode to IPC stream bytes, decode back, and assert per-cell equality.
// Self-roundtrip catches any framing / endianness / alignment bug without
// needing python or a pre-baked golden buffer.

namespace {

// Minimal helper: build an ArrowSchema for `struct<name: utf8, count: int64>`.
// Caller owns the schema and must call schema->release(schema) when done.
void make_schema(ArrowSchema *schema) {
	if (ArrowSchemaInitFromType(schema, NANOARROW_TYPE_STRUCT) != NANOARROW_OK) {
		throw std::runtime_error("schema init failed");
	}
	if (ArrowSchemaAllocateChildren(schema, 2) != NANOARROW_OK) {
		throw std::runtime_error("alloc children failed");
	}

	ArrowSchemaInitFromType(schema->children[0], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema->children[0], "name");
	ArrowSchemaInitFromType(schema->children[1], NANOARROW_TYPE_INT64);
	ArrowSchemaSetName(schema->children[1], "count");
}

// Build an ArrowArray with two rows: (name="alpha", count=10), (name="beta", count=42).
// Caller owns array and must call array->release(array) when done.
void make_array_two_rows(const ArrowSchema *schema, ArrowArray *array) {
	ArrowError err {};
	if (ArrowArrayInitFromSchema(array, schema, &err) != NANOARROW_OK) {
		throw std::runtime_error("array init failed");
	}
	if (ArrowArrayStartAppending(array) != NANOARROW_OK) {
		throw std::runtime_error("start appending failed");
	}
	// Row 0
	ArrowStringView alpha {"alpha", 5};
	if (ArrowArrayAppendString(array->children[0], alpha) != NANOARROW_OK) {
		throw std::runtime_error("append name failed");
	}
	if (ArrowArrayAppendInt(array->children[1], 10) != NANOARROW_OK) {
		throw std::runtime_error("append count failed");
	}
	if (ArrowArrayFinishElement(array) != NANOARROW_OK) {
		throw std::runtime_error("finish element 0 failed");
	}
	// Row 1
	ArrowStringView beta {"beta", 4};
	ArrowArrayAppendString(array->children[0], beta);
	ArrowArrayAppendInt(array->children[1], 42);
	ArrowArrayFinishElement(array);

	if (ArrowArrayFinishBuildingDefault(array, &err) != NANOARROW_OK) {
		throw std::runtime_error("finish building failed");
	}
}

} // namespace

TEST_CASE("EncodeIpcStream + DecodeIpcStream roundtrip preserves schema and rows", "[gpl-boundary][arrow_ipc]") {
	ArrowSchema schema {};
	make_schema(&schema);
	ArrowArray array {};
	make_array_two_rows(&schema, &array);

	// Encode
	std::vector<uint8_t> bytes = EncodeIpcStream(&schema, &array, 1);
	REQUIRE_FALSE(bytes.empty());
	// schema is consumed-but-not-released; we still own it for cleanup. Same for array.

	// Decode and confirm we get the same shape back.
	IpcStreamDecoder decoder(bytes.data(), bytes.size());
	ArrowSchema decoded_schema {};
	decoder.GetSchema(&decoded_schema);
	REQUIRE(decoded_schema.n_children == 2);
	REQUIRE(std::string(decoded_schema.children[0]->name) == "name");
	REQUIRE(std::string(decoded_schema.children[1]->name) == "count");
	// Format codes: utf8 = "u", int64 = "l"
	REQUIRE(std::string(decoded_schema.children[0]->format) == "u");
	REQUIRE(std::string(decoded_schema.children[1]->format) == "l");

	ArrowArray decoded_array {};
	REQUIRE(decoder.NextBatch(&decoded_array));
	REQUIRE(decoded_array.length == 2);
	REQUIRE(decoded_array.n_children == 2);

	// Verify int64 column
	const auto *count_buf = static_cast<const int64_t *>(decoded_array.children[1]->buffers[1]);
	REQUIRE(count_buf[0] == 10);
	REQUIRE(count_buf[1] == 42);

	// Verify utf8 column. Layout: buffers[1] = int32 offsets (n+1), buffers[2] = bytes.
	const auto *offsets = static_cast<const int32_t *>(decoded_array.children[0]->buffers[1]);
	const auto *data_bytes = static_cast<const char *>(decoded_array.children[0]->buffers[2]);
	REQUIRE(offsets[0] == 0);
	REQUIRE(std::string(data_bytes + offsets[0], offsets[1] - offsets[0]) == "alpha");
	REQUIRE(std::string(data_bytes + offsets[1], offsets[2] - offsets[1]) == "beta");

	// No more batches.
	ArrowArray empty {};
	REQUIRE_FALSE(decoder.NextBatch(&empty));

	// Cleanup
	if (decoded_schema.release) {
		decoded_schema.release(&decoded_schema);
	}
	if (decoded_array.release) {
		decoded_array.release(&decoded_array);
	}
	if (array.release) {
		array.release(&array);
	}
	if (schema.release) {
		schema.release(&schema);
	}
}

TEST_CASE("EncodeIpcStream produces identical bytes within a single process", "[gpl-boundary][arrow_ipc]") {
	// Two back-to-back calls in the same process MUST produce identical
	// output. This is a much weaker property than cross-build determinism —
	// FlatBuffers builders may have process-state-dependent output across
	// process invocations. Cross-build determinism is verified separately
	// in Phase 6 by diffing committed golden hexdumps after a clean rebuild.
	ArrowSchema schema {};
	make_schema(&schema);
	ArrowArray array_a {};
	make_array_two_rows(&schema, &array_a);
	ArrowArray array_b {};
	make_array_two_rows(&schema, &array_b);

	auto bytes_a = EncodeIpcStream(&schema, &array_a, 1);
	auto bytes_b = EncodeIpcStream(&schema, &array_b, 1);
	REQUIRE(bytes_a == bytes_b);

	if (array_a.release) {
		array_a.release(&array_a);
	}
	if (array_b.release) {
		array_b.release(&array_b);
	}
	if (schema.release) {
		schema.release(&schema);
	}
}

TEST_CASE("DecodeIpcStream rejects garbage bytes", "[gpl-boundary][arrow_ipc]") {
	const std::string garbage = "this is definitely not arrow ipc";
	IpcStreamDecoder decoder(garbage.data(), garbage.size());
	ArrowSchema schema {};
	REQUIRE_THROWS(decoder.GetSchema(&schema));
}

TEST_CASE("DecodeIpcStream rejects truncated bytes", "[gpl-boundary][arrow_ipc]") {
	// Build a real stream and truncate it at a point we know lies inside the
	// batch body (not at the trailing EOS marker — nanoarrow_ipc's decoder
	// tolerates trailing-byte loss as "no more batches"; we need to land
	// inside the actual record-batch payload so the decoder is forced to
	// surface a real error).
	ArrowSchema schema {};
	make_schema(&schema);
	ArrowArray array {};
	make_array_two_rows(&schema, &array);
	auto bytes = EncodeIpcStream(&schema, &array, 1);
	// Empirical sanity: this fixture produces ~700+ bytes total. Cutting at
	// 60% always lands well inside the batch body (after the schema message
	// at the head, before the EOS at the tail).
	REQUIRE(bytes.size() > 200);
	const std::size_t cut = (bytes.size() * 60) / 100;

	std::vector<uint8_t> truncated(bytes.begin(), bytes.begin() + cut);
	IpcStreamDecoder decoder(truncated.data(), truncated.size());

	// Schema must decode (it's at the head of the stream).
	ArrowSchema decoded_schema {};
	REQUIRE_NOTHROW(decoder.GetSchema(&decoded_schema));

	// NextBatch must throw — the batch body is genuinely incomplete here.
	// Unconditional REQUIRE_THROWS — no try/catch dance that could let the
	// test pass vacuously.
	ArrowArray decoded_array {};
	REQUIRE_THROWS(decoder.NextBatch(&decoded_array));

	if (decoded_schema.release) {
		decoded_schema.release(&decoded_schema);
	}
	if (array.release) {
		array.release(&array);
	}
	if (schema.release) {
		schema.release(&schema);
	}
}
