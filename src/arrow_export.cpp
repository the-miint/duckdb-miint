#include "arrow_export.hpp"

#include <cstring>

namespace miint {

namespace {

//! Backing storage for one exported Arrow array, owned by its `private_data`.
//!
//! These exports are always the simplest possible arrays: dense, unsliced and
//! non-null, which is the layout a borrowing consumer can read in place (sc's
//! importer rejects anything else outright). So there is never a validity bitmap
//! to build, and `buffers[0]` is always NULL -- Arrow's encoding for "this array
//! has no nulls at all".
struct BufferBag {
	std::vector<int64_t> i64;
	std::vector<double> f64;
	std::vector<int32_t> offsets;
	std::vector<char> chars;
	const void *buffers[3] = {nullptr, nullptr, nullptr};
};

// function to release c++ memory allocated for arrow array and schema
// called by the consumer of the arrow array and schema
void ReleaseArray(ArrowArray *array) {
	if (!array->release) {
		return;
	}
	// cast it back to the original type and delete it. The consumer doesn't know the type, so it can't delete it
	// directly.
	delete static_cast<BufferBag *>(array->private_data);
	array->private_data = nullptr;
	// Setting release to NULL is how a consumer detects an already-released
	// array; it makes double-release a no-op rather than a double-free.
	array->release = nullptr;
}

// nothing to free since this memory are strings in progmem
// written in at compile time, so we just set the release to nullptr to avoid double free
void ReleaseSchema(ArrowSchema *schema) {
	if (!schema->release) {
		return;
	}
	// `format` points at a static string literal and `name` is NULL, so there
	// is nothing to free — only the released marker to set.
	schema->release = nullptr;
}

void InitSchema(ArrowSchema &schema, const char *format) {
	schema.format = format;
	schema.name = nullptr;
	schema.metadata = nullptr;
	// 0 rather than ARROW_FLAG_NULLABLE: these arrays carry no nulls, and sc
	// rejects any that claim to.
	schema.flags = 0;
	schema.n_children = 0;
	schema.children = nullptr;
	schema.dictionary = nullptr;
	schema.release = ReleaseSchema;
	schema.private_data = nullptr;
}

// void function fills out standard ArrowArray fields and sets the release callback to ReleaseArray
// intakes an ArrowArray object somewhere in memory and mutates it
void InitArray(ArrowArray &array, int64_t length, int64_t n_buffers, BufferBag *bag) {
	array.length = length;
	array.null_count = 0;
	array.offset = 0;
	array.n_buffers = n_buffers;
	array.n_children = 0;
	array.buffers = bag->buffers;
	array.children = nullptr;
	array.dictionary = nullptr;
	array.release = ReleaseArray;
	array.private_data = bag;
}

} // namespace

void ExportInt64(ArrowArray &array, ArrowSchema &schema, std::vector<int64_t> values) {
	auto *bag = new BufferBag();
	bag->i64 = std::move(values);
	bag->buffers[0] = nullptr;
	bag->buffers[1] = bag->i64.data();
	InitArray(array, static_cast<int64_t>(bag->i64.size()), 2, bag);
	InitSchema(schema, "l");
}

//! Export `values` as an Arrow `Float64` array (format "g"): [validity, data].
void ExportFloat64(ArrowArray &array, ArrowSchema &schema, std::vector<double> values) {
	auto *bag = new BufferBag();
	bag->f64 = std::move(values);
	bag->buffers[0] = nullptr;
	bag->buffers[1] = bag->f64.data();
	InitArray(array, static_cast<int64_t>(bag->f64.size()), 2, bag);
	InitSchema(schema, "g");
}

//! Export `values` as an Arrow `Utf8` array (format "u"): [validity, offsets, data].
//!
//! Offsets hold `n + 1` entries and count BYTES, not characters; row `i` is
//! `chars[offsets[i] .. offsets[i + 1]]`.
void ExportUtf8(ArrowArray &array, ArrowSchema &schema, const std::vector<std::string> &values) {
	// one level of indirection to bag sitting on the heap, so that the consumer of the arrow array can free it when
	// done
	auto *bag = new BufferBag();
	bag->offsets.reserve(values.size() + 1); // allocate raw contiguous bytes for the offsets, its always N + 1
	size_t total = 0;
	// &v used in a declaration, not assignment, here it refers
	// to an alias of the string in the vector, not a copy of it. The string is still owned by the vector, so we don't
	// need to free it. this is not getting the address of the string.
	for (const auto &v : values) {
		total += v.size(); // number of bytes in each string
	}
	// since we are reducing std:string type to raw bytes, we need to allocate a contiguous block of memory for the
	// chars
	// we looped through each string in the vector to get its size in bytes and sum them up to get the total size of the
	// chars buffer we need to allocate
	bag->chars.reserve(total);

	int32_t cursor = 0; // int32 is agreed upon arrow type for offsets
	// load offsets with first value 0, we always have the start of a string at byte 0
	bag->offsets.push_back(cursor);
	// loop through std::vector<std::string> &values, and for each string, we insert its string bytes into the chars
	// buffer, and update the cursor to point to the end of the string in bytes, and push that value into the offsets
	// buffer
	for (const auto &v : values) {
		// bag->chars.end() iterator (memory addr) to insert the string bytes into for chars buffer, v.begin() and
		// v.end() are iterators  to the start and end of the string in the values vector, so we are inserting the
		// string bytes into the chars buffer for example if string v is "apple", it goes to the address of 'a', reads
		// all the bytes sequentially until it hits the address of v.end(), and copies those exact ASCII/UTF-8 character
		// bytes
		// ('a', 'p', 'p', 'l', 'e') into the contiguous bag->chars memory block.
		bag->chars.insert(bag->chars.end(), v.begin(), v.end());
		// increment offset by chunk size of the string in bytes, so that the next offset points to the start of the
		// next string in the chars buffer
		cursor += static_cast<int32_t>(v.size());
		// push offset
		bag->offsets.push_back(cursor);
	}
	// no nulls, so validity buffer is absent (NULL pointer). The offsets and chars buffers are always present in sc
	// contract.
	bag->buffers[0] = nullptr;
	// buffers is an array of pointers. offsets.data() returns a pointer to the first element of the offsets vector,
	// which is a contiguous block of memory. We assign that pointer to buffers[1] so that the ArrowArray can access the
	// offsets buffer. arrow contract always calls for fixed width int32 offsets buffer for Utf8 arrays, so on the
	// consumer side they know how to read it.
	bag->buffers[1] = bag->offsets.data();
	// An all-empty-string column allocates nothing; a NULL data buffer with
	// every offset at 0 is well-formed, and sc's `borrow_utf8` short-circuits
	// on length 0 before it reads any pointer.
	bag->buffers[2] = bag->chars.empty() ? nullptr : bag->chars.data();
	InitArray(array, static_cast<int64_t>(values.size()), 3, bag);
	InitSchema(schema, "u");
}

void ReleaseIfLive(ArrowArray &array, ArrowSchema &schema) {
	if (array.release) {
		array.release(&array);
	}
	if (schema.release) {
		schema.release(&schema);
	}
}

} // namespace miint
