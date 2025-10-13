#include <SAMReader.hpp>
#include <htslib-1.22.1/htslib/sam.h>
#include <sys/resource.h>

namespace miint {
// Custom deleters for htslib types
SAMReader::SAMReader(const std::string &filename)
    : fp(sam_open(filename.c_str(), "r")), hdr(sam_hdr_read(fp.get())), aln(bam_init1()) {
	if (!fp) {
		// better error?
		throw std::runtime_error("Failed to open SAM file");
	}
	// we can get the end position from cigar at parse
	// hts_pos_t bam_endpos(const bam1_t *b);
	// 
	// sam_hdr_t* hdr = sam_hdr_init();

	// Add a reference sequence
	// sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "248956422", NULL);
	// sam_hdr_add_line(hdr, "SQ", "SN", "chr2", "LN", "242193529", NULL);
	if (!hdr) {
		throw std::runtime_error("No SAM header");
	}

	if (!aln) {
		throw std::runtime_error("Cannot initialize BAM record");
	}
}

std::vector<SAMRecord> SAMReader::read(const int n) const {
	std::vector<SAMRecord> records;
	records.reserve(n);

	// Read n records from the SAM file
	while (sam_read1(fp.get(), hdr.get(), aln.get()) >= 0) {
		auto rec = SAMRecord(aln.get(), hdr.get());
		records.emplace_back(std::move(rec));
	}
	return records;
}
}; // namespace miint
