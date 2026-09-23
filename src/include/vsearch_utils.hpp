#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace miint {

// vsearch's batch and whole-database entry points start opt_threads worker threads
// even when opt_threads is 1. A wasm build without -pthread cannot start any:
// pthread_create fails and vsearch's fatal() calls exit(1). There, the wrappers
// default to the single-threaded stand-ins in vsearch_serial.hpp instead.
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
constexpr bool kVsearchSerialByDefault = true;
#else
constexpr bool kVsearchSerialByDefault = false;
#endif

// vsearch result structs use fixed 1024-char label buffers. Validate at input
// boundaries to throw early rather than silently truncate.
inline void validate_label_length(const std::string &label, const char *context) {
	if (label.size() > 1023) {
		throw std::runtime_error(std::string(context) + ": ID '" + label.substr(0, 40) +
		                         "...' exceeds 1023-character limit (got " + std::to_string(label.size()) + ")");
	}
}

// Convert RNA alphabet (U/u) to DNA (T/t) for vsearch compatibility.
// vsearch operates on DNA internally; RNA sequences must be normalized.
inline std::string normalize_rna(const std::string &seq) {
	std::string out = seq;
	for (auto &c : out) {
		if (c == 'U') {
			c = 'T';
		} else if (c == 'u') {
			c = 't';
		}
	}
	return out;
}

// Prepared batch of sequences for vsearch C batch APIs.
// Owns the normalized sequence data; pointer arrays are valid for its lifetime.
struct VsearchBatchArgs {
	std::vector<std::string> seqs_normalized;
	std::vector<const char *> seq_ptrs;
	std::vector<const char *> head_ptrs;
	std::vector<int> lens;
	std::vector<int> sizes;
	int count = 0;

	VsearchBatchArgs(const std::vector<std::string> &labels, const std::vector<std::string> &sequences) {
		count = static_cast<int>(labels.size());
		seqs_normalized.resize(count);
		seq_ptrs.resize(count);
		head_ptrs.resize(count);
		lens.resize(count);
		sizes.assign(count, 1);
		for (int i = 0; i < count; i++) {
			seqs_normalized[i] = normalize_rna(sequences[i]);
			seq_ptrs[i] = seqs_normalized[i].c_str();
			head_ptrs[i] = labels[i].c_str();
			lens[i] = static_cast<int>(seqs_normalized[i].size());
		}
	}
};

} // namespace miint
