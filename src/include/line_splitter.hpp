#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string_view>
#include <vector>

namespace miint {

// Splits a byte stream into lines without holding the whole input in memory.
//
// The source is a callback that writes up to `len` bytes into `dst` and returns how many
// it wrote, or 0 at end of input. A read failure is reported by the callback throwing.
// Nothing here touches DuckDB, so it can be unit tested directly; table functions feed it
// from a DuckDB file handle.
class LineSplitter {
public:
	using ReadFn = std::function<size_t(char *dst, size_t len)>;

	explicit LineSplitter(ReadFn read, size_t buffer_size = 1 << 16);

	// Advance to the next line. Returns false at end of input. `line` excludes the
	// terminating '\n' and a '\r' immediately before it, and stays valid until the next
	// call. A final line with no terminating '\n' is still returned, without a trailing
	// '\r'; an empty input has no lines, and an empty line in the middle of the input is
	// returned as empty.
	bool Next(std::string_view &line);

	// 1-based number of the line most recently returned by Next.
	uint64_t LineNumber() const {
		return line_number_;
	}

private:
	ReadFn read_;
	std::vector<char> buffer_;
	size_t begin_ = 0; // first byte not yet returned as part of a line
	size_t end_ = 0;   // one past the last byte read into the buffer
	bool eof_ = false;
	uint64_t line_number_ = 0;
};

// Split `line` on '\t' into `fields` (cleared first). Empty fields are kept, so a line of
// n tabs always yields n + 1 fields.
void SplitOnTabs(std::string_view line, std::vector<std::string_view> &fields);

} // namespace miint
