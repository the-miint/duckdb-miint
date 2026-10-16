#include "line_splitter.hpp"

#include <algorithm>
#include <cstring>

namespace miint {

LineSplitter::LineSplitter(ReadFn read, size_t buffer_size)
    : read_(std::move(read)), buffer_(std::max<size_t>(buffer_size, 1)) {
}

bool LineSplitter::Next(std::string_view &line) {
	// Bytes in [begin_, scan_from) are known to hold no '\n', so a refill rescans only
	// what it added.
	size_t scan_from = begin_;
	while (true) {
		const char *base = buffer_.data();
		const void *newline = std::memchr(base + scan_from, '\n', end_ - scan_from);
		if (newline != nullptr || (eof_ && begin_ < end_)) {
			size_t stop = newline != nullptr ? static_cast<size_t>(static_cast<const char *>(newline) - base) : end_;
			size_t line_end = stop;
			if (line_end > begin_ && base[line_end - 1] == '\r') {
				--line_end;
			}
			line = std::string_view(base + begin_, line_end - begin_);
			begin_ = newline != nullptr ? stop + 1 : end_;
			++line_number_;
			return true;
		}
		if (eof_) {
			return false;
		}

		// Move the partial line to the front, grow if it fills the buffer, and refill.
		size_t partial = end_ - begin_;
		if (begin_ > 0) {
			std::memmove(buffer_.data(), base + begin_, partial);
		}
		begin_ = 0;
		end_ = partial;
		scan_from = partial;
		if (end_ == buffer_.size()) {
			buffer_.resize(buffer_.size() * 2);
		}
		size_t n = read_(buffer_.data() + end_, buffer_.size() - end_);
		if (n == 0) {
			eof_ = true;
		}
		end_ += n;
	}
}

void SplitOnTabs(std::string_view line, std::vector<std::string_view> &fields) {
	fields.clear();
	size_t start = 0;
	while (true) {
		size_t tab = line.find('\t', start);
		if (tab == std::string_view::npos) {
			fields.push_back(line.substr(start));
			return;
		}
		fields.push_back(line.substr(start, tab - start));
		start = tab + 1;
	}
}

} // namespace miint
