#pragma once
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "SAMRecord.hpp"
#include <htslib-1.22.1/htslib/sam.h>

namespace miint {
struct SAMFileDeleter {
	void operator()(samFile *fp) const {
		if (fp) {
			sam_close(fp);
		}
	}
};

struct SAMHeaderDeleter {
	void operator()(sam_hdr_t *hdr) const {
		if (hdr) {
			sam_hdr_destroy(hdr);
		}
	}
};

struct BAMRecordDeleter {
	void operator()(bam1_t *aln) const {
		if (aln) {
			bam_destroy1(aln);
		}
	}
};

// Type aliases for smart pointers
using SAMFilePtr = std::unique_ptr<samFile, SAMFileDeleter>;
using SAMHeaderPtr = std::unique_ptr<sam_hdr_t, SAMHeaderDeleter>;
using BAMRecordPtr = std::unique_ptr<bam1_t, BAMRecordDeleter>;

// SAMReader: reads SAM/BAM/CRAM files using htslib
// Thread safety: Multiple SAMReader instances can safely read different files concurrently.
// A single SAMReader instance is NOT thread-safe for concurrent calls on the same instance.
class SAMReader {
public:
	// Constructor for SAM files with headers
	explicit SAMReader(const std::string &filename);

	// Constructor for headerless SAM files
	// Creates a synthetic header from the provided reference map
	explicit SAMReader(const std::string &filename, const std::unordered_map<std::string, uint64_t> &references);

	std::vector<SAMRecord> read(const int n) const;

private:
	SAMFilePtr fp;
	SAMHeaderPtr hdr;
	BAMRecordPtr aln;
};
}; // namespace miint
