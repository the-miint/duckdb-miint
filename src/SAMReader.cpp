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
		(void)hclose(hf);
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
		(void)hclose(hf);
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
// POSIX-only path for subprocess pipe input. No in-tree consumers today.
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

void SAMReader::set_threads(int n) {
	if (n > 1 && fp) {
		// Adds an HTSlib worker pool that decompresses BGZF blocks ahead of the parser.
		// Blocks are still delivered in order, so read() output is unchanged. No-op for
		// uncompressed SAM (nothing to decompress).
		hts_set_threads(fp.get(), n);
	}
}

// sam_read1 returns >= 0 on success, -1 at end of stream, and < -1 on a read error (truncated
// file, unparseable record). Both read paths used to test only `>= 0`, which conflated the two:
// a malformed record ended the batch, and since the table functions call read() again for the
// next batch, the scan RESUMED after it -- the bad record was silently dropped and the query
// succeeded. samtools exits 1 on the same input. A corrupt file must not read back as a subset
// of itself with nothing to signal the loss.
// True when sam_read1 can actually read records out of this format. For anything else, its -3
// return is a statement about the FILE TYPE rather than about a record.
static bool FormatCarriesRecords(const htsFile *fp) {
	switch (fp->format.format) {
	case bam:
	case cram:
	case sam:
	case fasta_format:
	case fastq_format:
		return true;
	default:
		// empty_format, text_format, unknown_format, and every index/variant format.
		return false;
	}
}

static void CheckReadError(htsFile *fp, int ret) {
	if (ret >= -1) {
		return; // -1: end of stream
	}
	// -3 is overloaded, and getting this wrong silently truncates files.
	//
	// For a format sam_read1 cannot read records from, -3 is decided before any record is examined
	// -- empty_format (EPIPE) for a zero-byte or empty-after-decompression file, EFTYPE for
	// anything unidentifiable. Those are valid inputs holding zero records, and the headerless
	// constructor reads them deliberately, since its header comes from the reference map rather
	// than from the file (test_SAMReader.cpp pins empty and whitespace-only input that way).
	//
	// But once the format IS readable, -3 also reports genuine corruption: bam_read1 returns -3
	// for a record truncated between its length prefix and its 32-byte core, and sam_read1_bam for
	// a tid/mtid outside the header's range. Exempting the bare return code therefore swallowed
	// truncated BAM -- precisely the corruption this check exists to catch, and invisible to the
	// text-SAM tests because those fail with -2. Gate on the format, which is what the exemption
	// is actually about.
	if (ret == -3 && (!fp || !FormatCarriesRecords(fp))) {
		return;
	}
	// fp->lineno is htslib's own locator for text SAM -- it reports "Parse error at line %lld" from
	// this same failure -- and htsFile is a deliberately public struct whose fields samtools reads
	// directly, so this is portable (including under Emscripten; no POSIX involved).
	//
	// Only trust it for single-threaded text SAM. It is 0 for BAM, and stale under the
	// multi-threaded SAM parser: hts_set_threads on a text SAM routes to sam_set_threads (not to
	// bgzf), whose dispatcher consumes raw blocks without going through hts_getline, leaving lineno
	// frozen near the header. Reporting "at line 3" for a bad record at line 900,000 is worse than
	// reporting no line at all. fp->state is that parser's SAM_state, so it doubles as the signal.
	std::string where;
	if (fp && fp->format.format == sam && !fp->state && fp->lineno > 0) {
		where = " at line " + std::to_string(fp->lineno);
	}
	throw std::runtime_error("Failed to read SAM/BAM record" + where + ": the file is truncated or malformed" +
	                         " (htslib error " + std::to_string(ret) + ")");
}

const bam1_t *SAMReader::read_raw() {
	int ret = sam_read1(fp.get(), hdr.get(), aln.get());
	CheckReadError(fp.get(), ret);
	if (ret >= 0) {
		return aln.get();
	}
	return nullptr;
}

SAMRecordBatch SAMReader::read(const int n) {
	SAMRecordBatch batch;
	batch.reserve(n);

	// Read up to n records from the SAM file
	for (int i = 0; i < n; ++i) {
		int ret = sam_read1(fp.get(), hdr.get(), aln.get());
		CheckReadError(fp.get(), ret);
		if (ret < 0) {
			break; // -1: end of stream
		}
		// Note: We cannot validate that references in headerless files match the expected set
		// because htslib automatically marks reads with unknown references as unmapped (tid=-1, FLAG 0x4)
		// making them indistinguishable from genuinely unmapped reads. Users must ensure their
		// reference_lengths table includes all references present in the data files.
		sam_utils::parse_record_to_batch(aln.get(), hdr.get(), batch, include_seq_qual);
	}
	return batch;
}
}; // namespace miint
