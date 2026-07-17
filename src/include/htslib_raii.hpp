#pragma once
//
// RAII smart-pointer wrappers for the htslib handles used by the SAM/BAM/uBAM
// COPY writers. Extracted from copy_sam.cpp so copy_ubam.cpp can reuse the exact
// same lifetime management without duplicating the deleters.

#include <htslib-1.22.1/htslib/sam.h>
#include <memory>

namespace duckdb {

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

using SAMFilePtr = std::unique_ptr<samFile, SAMFileDeleter>;
using SAMHeaderPtr = std::unique_ptr<sam_hdr_t, SAMHeaderDeleter>;
using BAMRecordPtr = std::unique_ptr<bam1_t, BAMRecordDeleter>;

} // namespace duckdb
