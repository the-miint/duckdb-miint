#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "line_splitter.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using miint::LineSplitter;

namespace {

// A source that hands out `input` at most `chunk` bytes per read, so a test controls
// where read boundaries fall relative to line boundaries.
LineSplitter::ReadFn ChunkedSource(std::string input, size_t chunk) {
	auto pos = std::make_shared<size_t>(0);
	auto data = std::make_shared<std::string>(std::move(input));
	return [pos, data, chunk](char *dst, size_t len) -> size_t {
		size_t n = std::min({chunk, len, data->size() - *pos});
		data->copy(dst, n, *pos);
		*pos += n;
		return n;
	};
}

struct Line {
	std::string text;
	uint64_t number;
	bool operator==(const Line &other) const {
		return text == other.text && number == other.number;
	}
};

std::vector<Line> AllLines(const std::string &input, size_t chunk, size_t buffer_size = 1 << 16) {
	LineSplitter splitter(ChunkedSource(input, chunk), buffer_size);
	std::vector<Line> lines;
	std::string_view line;
	while (splitter.Next(line)) {
		lines.push_back({std::string(line), splitter.LineNumber()});
	}
	return lines;
}

} // namespace

TEST_CASE("LineSplitter returns lines and their 1-based numbers", "[line_splitter]") {
	// A missing final newline is common in hand-edited files; dropping that line would lose
	// the last record without any error.
	CHECK(AllLines("a\nbb\nccc", 1024) == std::vector<Line> {{"a", 1}, {"bb", 2}, {"ccc", 3}});
	CHECK(AllLines("a\nbb\n", 1024) == std::vector<Line> {{"a", 1}, {"bb", 2}});
	CHECK(AllLines("", 1024).empty());
}

TEST_CASE("LineSplitter keeps an empty line in the middle of the input", "[line_splitter]") {
	// Readers reject blank lines, and the error has to name the right line, so the empty
	// line must be reported rather than skipped -- and must still advance the count.
	CHECK(AllLines("a\n\nb\n", 1024) == std::vector<Line> {{"a", 1}, {"", 2}, {"b", 3}});
}

TEST_CASE("LineSplitter strips CRLF line endings, and only those", "[line_splitter]") {
	CHECK(AllLines("a\r\nb\r\n", 1024) == std::vector<Line> {{"a", 1}, {"b", 2}});
	// A carriage return inside a line is data, not a line ending.
	CHECK(AllLines("a\rb\n", 1024) == std::vector<Line> {{"a\rb", 1}});
	// A read boundary between '\r' and '\n' must not leave the '\r' on the line.
	CHECK(AllLines("ab\r\ncd\n", 3) == std::vector<Line> {{"ab", 1}, {"cd", 2}});
	// A CRLF file whose last line lost its '\n' still has a '\r' there, which is stripped too.
	CHECK(AllLines("a\r\nb\r", 1024) == std::vector<Line> {{"a", 1}, {"b", 2}});
}

TEST_CASE("LineSplitter output does not depend on where reads or the buffer end", "[line_splitter]") {
	// Lines that straddle a read boundary, and lines longer than the whole buffer, are the
	// cases a buffered splitter gets wrong; compare every small chunk and buffer size
	// against reading the input in one piece.
	const std::string input = "first line\n\nCRLF line\r\na line much longer than the tiny buffer\nlast";
	const auto expected = AllLines(input, input.size());
	REQUIRE(
	    expected ==
	    std::vector<Line> {
	        {"first line", 1}, {"", 2}, {"CRLF line", 3}, {"a line much longer than the tiny buffer", 4}, {"last", 5}});
	for (size_t chunk = 1; chunk <= 7; chunk++) {
		for (size_t buffer : {size_t(1), size_t(2), size_t(4), size_t(16)}) {
			INFO("chunk " << chunk << ", buffer " << buffer);
			CHECK(AllLines(input, chunk, buffer) == expected);
		}
	}
}

TEST_CASE("LineSplitter propagates a read failure", "[line_splitter]") {
	LineSplitter splitter([](char *, size_t) -> size_t { throw std::runtime_error("upstream read failed"); });
	std::string_view line;
	CHECK_THROWS_WITH(splitter.Next(line), ContainsSubstring("upstream read failed"));
}

TEST_CASE("SplitOnTabs keeps empty fields", "[line_splitter]") {
	std::vector<std::string_view> fields;
	miint::SplitOnTabs("a\t\tb", fields);
	CHECK(fields == std::vector<std::string_view> {"a", "", "b"});
	miint::SplitOnTabs("a\t", fields);
	CHECK(fields == std::vector<std::string_view> {"a", ""});
	miint::SplitOnTabs("", fields);
	CHECK(fields == std::vector<std::string_view> {""});
}
