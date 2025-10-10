#pragma once
#include <cstdint>
#include <string>
#include <sstream>
#include <htslib-1.22.1/htslib/sam.h>

namespace miint {

/* Bowtie2 optional flags (https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output)

 AS:i:<N>
 Alignment score. Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode).
 Only present if SAM record is for an aligned read.

 XS:i:<N>
 Alignment score for the best-scoring alignment found other than the alignment reported. Can be
 negative. Can be greater than 0 in --local mode (but not in --end-to-end mode). Only present if the SAM
 record is for an aligned read and more than one alignment was found for the read. Note that, when the
 read is part of a concordantly-aligned pair, this score could be greater than AS:i.

 YS:i:<N>
 Alignment score for opposite mate in the paired-end alignment. Only present if the SAM record is for a
 read that aligned as part of a paired-end alignment.

 XN:i:<N>
 The number of ambiguous bases in the reference covering this alignment. Only present if SAM record is
 for an aligned read.

 XM:i:<N>
 The number of mismatches in the alignment. Only present if SAM record is for an aligned read.

 XO:i:<N>
 The number of gap opens, for both read and reference gaps, in the alignment. Only present if SAM record
 is for an aligned read.

 XG:i:<N>
 The number of gap extensions, for both read and reference gaps, in the alignment. Only present if SAM
 record is for an aligned read.

 NM:i:<N>
 The edit distance; that is, the minimal number of one-nucleotide edits (substitutions, insertions and
 deletions) needed to transform the read string into the reference string. Only present if SAM record is
 for an aligned read.

 YF:Z:<S>  NOTE: Not included in parser as we are representing _aligned_ reads
 String indicating reason why the read was filtered out. See also: Filtering. Only appears for reads
 that were filtered out.

 YT:Z:<S>
 Value of UU indicates the read was not part of a pair. Value of CP indicates the read was part of a
 pair and the pair aligned concordantly. Value of DP indicates the read was part of a pair and the pair
 aligned discordantly. Value of UP indicates the read was part of a pair but the pair failed to aligned
 either concordantly or discordantly.

 MD:Z:<S>
 A string representation of the mismatched reference bases in the alignment. See SAM Tags format
 specification for details. Only present if SAM record is for an aligned read.

 minimap2 flags (https://lh3.github.io/minimap2/minimap2.html)
 NOTE: these are described as relative to PAF, so we are only representing those below which are
 also part of the SAM spec. Only listing the entries below which do not appear in the bowtie2
 listing

 SA	Z	List of other supplementary alignments (with approximate CIGAR strings)
*/
enum class SAMRecordField {
	READ_ID = 0,
	FLAGS,
	REFERENCE,
	POSITION,
	MAPQ,
	CIGAR,
	MATE_REFERENCE,
	MATE_POSITION,
	TEMPLATE_LENGTH,
	TAG_AS,
	TAG_XS,
	TAG_YS,
	TAG_XN,
	TAG_XM,
	TAG_XO,
	TAG_XG,
	TAG_NM,
	TAG_YT,
	TAG_MD,
	TAG_SA
};

class SAMRecord {
private:
	int64_t get_int_tag(const bam1_t *aln, const char *tag, int64_t default_val = -1) {
		uint8_t *aux = bam_aux_get(aln, tag);
		if (!aux) {
			return default_val;
		}
		return bam_aux2i(aux);
	}

	// Helper to get string tag value, returns empty string if not found
	std::string get_string_tag(const bam1_t *aln, const char *tag) {
		uint8_t *aux = bam_aux_get(aln, tag);
		if (!aux) {
			return "";
		}
		return std::string(bam_aux2Z(aux));
	}

	// Helper to convert CIGAR to string
	std::string cigar_to_string(const bam1_t *aln) {
		std::ostringstream oss;
		uint32_t *cigar = reinterpret_cast<uint32_t *>(aln->data + aln->core.l_qname);
		for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
			uint32_t op_len = bam_cigar_oplen(cigar[i]);
			char op_chr = bam_cigar_opchr(cigar[i]);
			oss << op_len << op_chr;
		}
		return oss.str();
	}

public:
	std::string read_id;
	uint16_t flags;
	std::string reference;
	int64_t position;
	uint8_t mapq;
	std::string cigar;
	std::string mate_reference;
	int64_t mate_position;
	int64_t template_length;
	int64_t tag_as;
	int64_t tag_xs;
	int64_t tag_ys;
	int64_t tag_xn;
	int64_t tag_xm;
	int64_t tag_xo;
	int64_t tag_xg;
	int64_t tag_nm;
	std::string tag_yt;
	std::string tag_md;
	std::string tag_sa;

	explicit SAMRecord(const bam1_t *aln, const sam_hdr_t *hdr) noexcept {
		// bam_get_qname performs a c-style cast
		// recommended approach from claude for a reinterpret cast
		read_id = std::string(reinterpret_cast<char *>(aln->data), aln->core.l_qname - 1);
		flags = aln->core.flag;

		if (aln->core.tid >= 0) {
			reference = sam_hdr_tid2name(hdr, aln->core.tid);
		} else {
			reference = "*";
		}
		// Position (POS) - convert from 0-based to 1-based
		position = aln->core.pos >= 0 ? aln->core.pos + 1 : 0;

		// Mapping quality (MAPQ)
		mapq = aln->core.qual;

		// CIGAR string
		cigar = cigar_to_string(aln);

		// Mate reference name (RNEXT)
		if (aln->core.mtid >= 0) {
			if (aln->core.mtid == aln->core.tid) {
				mate_reference = "=";
			} else {
				mate_reference = sam_hdr_tid2name(hdr, aln->core.mtid);
			}
		} else {
			mate_reference = "*";
		}

		// Mate position (PNEXT) - convert from 0-based to 1-based
		mate_position = aln->core.mpos >= 0 ? aln->core.mpos + 1 : 0;

		// Template length (TLEN)
		template_length = aln->core.isize;

		// Optional tags
		tag_as = get_int_tag(aln, "AS");
		tag_xs = get_int_tag(aln, "XS");
		tag_ys = get_int_tag(aln, "YS");
		tag_xn = get_int_tag(aln, "XN");
		tag_xm = get_int_tag(aln, "XM");
		tag_xo = get_int_tag(aln, "XO");
		tag_xg = get_int_tag(aln, "XG");
		tag_nm = get_int_tag(aln, "NM");
		tag_yt = get_string_tag(aln, "YT");
		tag_md = get_string_tag(aln, "MD");
		tag_sa = get_string_tag(aln, "SA");
	}

	const std::string &GetString(const SAMRecordField field) const {
		switch (field) {
		case SAMRecordField::READ_ID:
			return read_id;
		case SAMRecordField::REFERENCE:
			return reference;
		case SAMRecordField::CIGAR:
			return cigar;
		case SAMRecordField::MATE_REFERENCE:
			return mate_reference;
		case SAMRecordField::TAG_YT:
			return tag_yt;
		case SAMRecordField::TAG_MD:
			return tag_md;
		case SAMRecordField::TAG_SA:
			return tag_sa;
		default:
			throw std::invalid_argument("Invalid field");
		}
	}

	const int64_t GetInt64(const SAMRecordField field) const {
		switch (field) {
		case SAMRecordField::POSITION:
			return position;
		case SAMRecordField::MATE_POSITION:
			return mate_position;
		case SAMRecordField::TEMPLATE_LENGTH:
			return template_length;
		case SAMRecordField::TAG_AS:
			return tag_as;
		case SAMRecordField::TAG_XS:
			return tag_xs;
		case SAMRecordField::TAG_YS:
			return tag_ys;
		case SAMRecordField::TAG_XN:
			return tag_xn;
		case SAMRecordField::TAG_XM:
			return tag_xm;
		case SAMRecordField::TAG_XO:
			return tag_xo;
		case SAMRecordField::TAG_XG:
			return tag_xg;
		case SAMRecordField::TAG_NM:
			return tag_nm;
		default:
			throw std::invalid_argument("Invalid field");
		}
	}

	const uint16_t GetUInt16(const SAMRecordField field) const {
		switch (field) {
		case SAMRecordField::FLAGS:
			return flags;
		default:
			throw std::invalid_argument("Invalid field");
		}
	}

	const uint8_t GetUInt8(const SAMRecordField field) const {
		switch (field) {
		case SAMRecordField::MAPQ:
			return mapq;
		default:
			throw std::invalid_argument("Invalid field");
		}
	}
};
}; // namespace miint
