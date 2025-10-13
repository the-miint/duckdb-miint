#pragma once
#include <vector>
#include <string>
#include <memory>
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

class SAMReader {
public:
	explicit SAMReader(const std::string &filename);
	std::vector<SAMRecord> read(const int n) const;

private:
	SAMFilePtr fp;
	SAMHeaderPtr hdr;
	BAMRecordPtr aln;
};
}; // namespace miint
