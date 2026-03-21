#include <SAMReader.hpp>
#include <htslib-1.22.1/htslib/sam.h>
#include <htslib-1.22.1/htslib/hfile.h>
#include <regex>
#ifndef _WIN32
#include <unistd.h>
#endif

namespace miint {

// Validate a reference map for SAM spec compliance.
static void ValidateReferenceMap(const std::unordered_map<std::string, uint64_t> &references) {
	if (references.empty()) {
		throw std::runtime_error("Reference map cannot be empty");
	}

	static const std::regex position_pattern(":[0-9]+(-[0-9]+)?$");

	for (const auto &[name, length] : references) {
		if (name.empty()) {
			throw std::runtime_error("Reference name cannot be empty");
		}
		if (length == 0) {
			throw std::runtime_error("Reference length cannot be zero");
		}
		if (name[0] == '*' || name[0] == '=') {
			throw std::runtime_error("Reference name cannot start with '*' or '='");
		}
		if (name.find('\t') != std::string::npos || name.find('\n') != std::string::npos) {
			throw std::runtime_error("Reference name contains invalid characters");
		}
		if (name.length() > 1024) {
			throw std::runtime_error("Reference name exceeds maximum length of 1024 characters");
		}
		if (std::regex_search(name, position_pattern)) {
			throw std::runtime_error(
			    "Reference name ends with position-like pattern (:<digits> or :<digits>-<digits>)");
		}
	}
}

// Read any existing header from the file, then build or augment it with the provided references.
// Returns a valid SAMHeaderPtr with @SQ lines.
static SAMHeaderPtr BuildHeaderFromReferences(samFile *fp,
                                              const std::unordered_map<std::string, uint64_t> &references) {
	// Try to read any existing header (even partial headers like @HD without @SQ).
	// This consumes header lines and positions the file pointer at the first alignment record.
	SAMHeaderPtr existing_hdr(sam_hdr_read(fp));

	if (existing_hdr && existing_hdr->n_targets == 0) {
		// Partial header exists (e.g., @HD line) but no @SQ lines — add them
		for (const auto &[name, length] : references) {
			std::string length_str = std::to_string(length);
			if (sam_hdr_add_line(existing_hdr.get(), "SQ", "SN", name.c_str(), "LN", length_str.c_str(), NULL) != 0) {
				throw std::runtime_error("Failed to add @SQ line to existing header");
			}
		}
		return existing_hdr;
	} else if (!existing_hdr || existing_hdr->n_targets == 0) {
		// No header — create synthetic
		std::string header_text;
		for (const auto &[name, length] : references) {
			header_text += "@SQ\tSN:" + name + "\tLN:" + std::to_string(length) + "\n";
		}
		SAMHeaderPtr hdr(sam_hdr_parse(header_text.length(), header_text.c_str()));
		if (!hdr) {
			throw std::runtime_error("Failed to parse SAM header");
		}
		return hdr;
	} else {
		// File has a complete header with @SQ lines — use as-is
		return existing_hdr;
	}
}

// Constructor for SAM files with headers
SAMReader::SAMReader(const std::string &filename, bool include_seq_qual, bool require_references)
    : fp(sam_open(filename.c_str(), "r")), hdr(sam_hdr_read(fp.get())), aln(bam_init1()),
      include_seq_qual(include_seq_qual) {
	if (!fp) {
		throw std::runtime_error("Failed to open SAM file");
	}

	// Validate header exists
	if (!hdr) {
		throw std::runtime_error("SAM file missing required header");
	}

	// Validate header has reference sequences (unless caller opts out, e.g. for uBAM files)
	if (require_references && hdr->n_targets == 0) {
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
	if (!fp) {
		throw std::runtime_error("Failed to open SAM file");
	}

	ValidateReferenceMap(references);
	hdr = BuildHeaderFromReferences(fp.get(), references);

	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

// Constructor for reading SAM/BAM via a pre-opened hFILE.
// Takes ownership of the hFILE — hts_close() will close it.
SAMReader::SAMReader(hFILE *hf, const std::string &name, bool include_seq_qual, bool require_references)
    : aln(bam_init1()), include_seq_qual(include_seq_qual) {
	if (!hf) {
		throw std::runtime_error("hFILE is null");
	}

	htsFile *hts_fp = hts_hopen(hf, name.c_str(), "r");
	if (!hts_fp) {
		hclose(hf);
		throw std::runtime_error("Failed to open hFILE as SAM/BAM stream");
	}
	fp.reset(hts_fp);

	hdr.reset(sam_hdr_read(fp.get()));
	if (!hdr) {
		throw std::runtime_error("SAM file missing required header");
	}

	if (require_references && hdr->n_targets == 0) {
		throw std::runtime_error("SAM file missing required header");
	}

	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

// Constructor for reading headerless SAM via a pre-opened hFILE.
SAMReader::SAMReader(hFILE *hf, const std::string &name, const std::unordered_map<std::string, uint64_t> &references,
                     bool include_seq_qual)
    : aln(bam_init1()), include_seq_qual(include_seq_qual) {
	if (!hf) {
		throw std::runtime_error("hFILE is null");
	}

	htsFile *hts_fp = hts_hopen(hf, name.c_str(), "r");
	if (!hts_fp) {
		hclose(hf);
		throw std::runtime_error("Failed to open hFILE as SAM/BAM stream");
	}
	fp.reset(hts_fp);

	ValidateReferenceMap(references);
	hdr = BuildHeaderFromReferences(fp.get(), references);

	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

#ifndef _WIN32
// Constructor for reading SAM from a file descriptor (e.g., pipe from subprocess).
// Only available on POSIX — used by Bowtie2Aligner for subprocess piping.
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
#endif

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
