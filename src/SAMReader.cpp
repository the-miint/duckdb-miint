#include <SAMReader.hpp>
#include <htslib-1.22.1/htslib/sam.h>
#include <htslib-1.22.1/htslib/hfile.h>
#include <sys/resource.h>
#include <regex>
#include <unistd.h>

namespace miint {
// Constructor for SAM files with headers
SAMReader::SAMReader(const std::string &filename, bool include_seq_qual)
    : fp(sam_open(filename.c_str(), "r")), hdr(sam_hdr_read(fp.get())), aln(bam_init1()),
      include_seq_qual(include_seq_qual) {
	if (!fp) {
		throw std::runtime_error("Failed to open SAM file");
	}

	// Validate header exists and is not empty
	if (!hdr || hdr->n_targets == 0) {
		throw std::runtime_error("SAM file missing required header");
	}

	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

// Constructor for headerless SAM files
SAMReader::SAMReader(const std::string &filename, const std::unordered_map<std::string, uint64_t> &references,
                     bool include_seq_qual)
    : fp(sam_open(filename.c_str(), "r")), aln(bam_init1()), include_seq_qual(include_seq_qual) {
	// Open file
	if (!fp) {
		throw std::runtime_error("Failed to open SAM file");
	}

	// Note: We don't call sam_hdr_read() to detect headers because it consumes the file position.
	// The user must use the correct constructor based on whether the file has a header.
	// The header-based constructor will validate that a header exists.

	// Validate reference map
	if (references.empty()) {
		throw std::runtime_error("Reference map cannot be empty");
	}

	for (const auto &[name, length] : references) {
		// Check for empty name
		if (name.empty()) {
			throw std::runtime_error("Reference name cannot be empty");
		}

		// Check for zero length
		if (length == 0) {
			throw std::runtime_error("Reference length cannot be zero");
		}

		// Check first character restrictions
		if (name[0] == '*' || name[0] == '=') {
			throw std::runtime_error("Reference name cannot start with '*' or '='");
		}

		// Check for invalid characters (tab, newline)
		if (name.find('\t') != std::string::npos || name.find('\n') != std::string::npos) {
			throw std::runtime_error("Reference name contains invalid characters");
		}

		// Check length limit
		if (name.length() > 1024) {
			throw std::runtime_error("Reference name exceeds maximum length of 1024 characters");
		}

		// Check for position-like pattern at end: :<digits> or :<digits>-<digits>
		// Regex pattern: :[0-9]+(-[0-9]+)?$
		std::regex position_pattern(":[0-9]+(-[0-9]+)?$");
		if (std::regex_search(name, position_pattern)) {
			throw std::runtime_error(
			    "Reference name ends with position-like pattern (:<digits> or :<digits>-<digits>)");
		}
	}

	// Create synthetic header
	// Build header text as string for efficiency with large reference sets
	std::string header_text;
	for (const auto &[name, length] : references) {
		header_text += "@SQ\tSN:" + name + "\tLN:" + std::to_string(length) + "\n";
	}

	// Parse header text once (much faster than repeated sam_hdr_add_line calls)
	hdr.reset(sam_hdr_parse(header_text.length(), header_text.c_str()));
	if (!hdr) {
		throw std::runtime_error("Failed to parse SAM header");
	}

	// Initialize alignment record
	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

// Constructor for reading SAM from a file descriptor (e.g., pipe from subprocess)
SAMReader::SAMReader(int fd, const std::string &name, bool include_seq_qual)
    : aln(bam_init1()), include_seq_qual(include_seq_qual) {

	// Wrap fd in HTSlib's hFILE
	hFILE *hfile = hdopen(fd, "r");
	if (!hfile) {
		close(fd);
		throw std::runtime_error("Failed to wrap file descriptor in hFILE");
	}

	// Wrap hFILE in htsFile for SAM parsing
	// "r" mode = auto-detect format (SAM/BAM/CRAM)
	htsFile *hts_fp = hts_hopen(hfile, name.c_str(), "r");
	if (!hts_fp) {
		hclose(hfile);
		throw std::runtime_error("Failed to open file descriptor as SAM stream");
	}
	fp.reset(hts_fp);

	// Read header
	hdr.reset(sam_hdr_read(fp.get()));
	if (!hdr) {
		throw std::runtime_error("Failed to read SAM header from stream");
	}

	// Initialize alignment record
	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

SAMRecordBatch SAMReader::read(const int n) {
	SAMRecordBatch batch;
	batch.reserve(n);

	// Read up to n records from the SAM file
	for (int i = 0; i < n && sam_read1(fp.get(), hdr.get(), aln.get()) >= 0; ++i) {
		// Note: We cannot validate that references in headerless files match the expected set
		// because htslib automatically marks reads with unknown references as unmapped (tid=-1, FLAG 0x4)
		// making them indistinguishable from genuinely unmapped reads. Users must ensure their
		// reference_lengths table includes all references present in the data files.
		sam_utils::parse_record_to_batch(aln.get(), hdr.get(), batch, include_seq_qual);
	}
	return batch;
}
}; // namespace miint
