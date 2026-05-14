#pragma once
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "SAMRecord.hpp"
#include <htslib-1.22.1/htslib/sam.h>
#include <htslib-1.22.1/htslib/hts_log.h>
#include <htslib-1.22.1/htslib/hfile.h>

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
	// require_references: if true (default), throws when header has no @SQ lines (n_targets == 0)
	//                     if false, accepts headers without @SQ lines (e.g., uBAM files)
	explicit SAMReader(const std::string &filename, bool include_seq_qual = false, bool require_references = true);

	// Constructor for headerless SAM files
	// Creates a synthetic header from the provided reference map
	// Note: Cannot validate that all references in the data are in the reference map due to htslib limitations
	explicit SAMReader(const std::string &filename, const std::unordered_map<std::string, uint64_t> &references,
	                   bool include_seq_qual = false);

	// Constructor for reading SAM/BAM via a pre-opened hFILE (e.g., DuckDB FileHandle for remote files).
	// Takes ownership of the hFILE — hts_close() will close it via the hFILE backend's close callback.
	// require_references: if true (default), throws when header has no @SQ lines
	explicit SAMReader(hFILE *hf, const std::string &name, bool include_seq_qual = false,
	                   bool require_references = true);

	// Constructor for reading headerless SAM via a pre-opened hFILE.
	// Creates a synthetic header from the provided reference map.
	explicit SAMReader(hFILE *hf, const std::string &name, const std::unordered_map<std::string, uint64_t> &references,
	                   bool include_seq_qual = false);

#ifndef _WIN32
	// Constructor for reading SAM from a file descriptor (e.g., pipe from subprocess)
	// Takes ownership of the file descriptor - it will be closed when the reader is destroyed
	// The SAM stream must include a header (POSIX only).
	explicit SAMReader(int fd, const std::string &name, bool include_seq_qual = false);
#endif

	// Read up to n records into a batch
	SAMRecordBatch read(const int n);

	// Number of reference sequences in the header (@SQ lines).
	// Returns 0 if no references (headerless/uBAM without synthetic header).
	int n_targets() const {
		return hdr ? hdr->n_targets : 0;
	}

	// Detect format (true if BAM, false if SAM/other).
	bool is_bam() const {
		if (!fp) {
			return false;
		}
		const htsFormat *fmt = hts_get_format(fp.get());
		return (fmt && fmt->format == bam);
	}

private:
	SAMFilePtr fp;
	SAMHeaderPtr hdr;
	BAMRecordPtr aln;
	bool include_seq_qual;
};
}; // namespace miint
