#pragma once
#include <cstdint>
#include <string>
#include <sstream>
#include <vector>
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
	STOP_POSITION,
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

// SOA (Struct of Arrays) layout for efficient batch processing
// Stores a batch of SAM records with fields organized by column for better cache locality
struct SAMRecordBatch {
	std::vector<std::string> read_ids;
	std::vector<uint16_t> flags;
	std::vector<std::string> references;
	std::vector<int64_t> positions;
	std::vector<int64_t> stop_positions;
	std::vector<uint8_t> mapqs;
	std::vector<std::string> cigars;
	std::vector<std::string> mate_references;
	std::vector<int64_t> mate_positions;
	std::vector<int64_t> template_lengths;
	std::vector<int64_t> tag_as_values;
	std::vector<int64_t> tag_xs_values;
	std::vector<int64_t> tag_ys_values;
	std::vector<int64_t> tag_xn_values;
	std::vector<int64_t> tag_xm_values;
	std::vector<int64_t> tag_xo_values;
	std::vector<int64_t> tag_xg_values;
	std::vector<int64_t> tag_nm_values;
	std::vector<std::string> tag_yt_values;
	std::vector<std::string> tag_md_values;
	std::vector<std::string> tag_sa_values;

	size_t size() const {
		return read_ids.size();
	}

	bool empty() const {
		return read_ids.empty();
	}

	void reserve(size_t n) {
		read_ids.reserve(n);
		flags.reserve(n);
		references.reserve(n);
		positions.reserve(n);
		stop_positions.reserve(n);
		mapqs.reserve(n);
		cigars.reserve(n);
		mate_references.reserve(n);
		mate_positions.reserve(n);
		template_lengths.reserve(n);
		tag_as_values.reserve(n);
		tag_xs_values.reserve(n);
		tag_ys_values.reserve(n);
		tag_xn_values.reserve(n);
		tag_xm_values.reserve(n);
		tag_xo_values.reserve(n);
		tag_xg_values.reserve(n);
		tag_nm_values.reserve(n);
		tag_yt_values.reserve(n);
		tag_md_values.reserve(n);
		tag_sa_values.reserve(n);
	}

	void clear() {
		read_ids.clear();
		flags.clear();
		references.clear();
		positions.clear();
		stop_positions.clear();
		mapqs.clear();
		cigars.clear();
		mate_references.clear();
		mate_positions.clear();
		template_lengths.clear();
		tag_as_values.clear();
		tag_xs_values.clear();
		tag_ys_values.clear();
		tag_xn_values.clear();
		tag_xm_values.clear();
		tag_xo_values.clear();
		tag_xg_values.clear();
		tag_nm_values.clear();
		tag_yt_values.clear();
		tag_md_values.clear();
		tag_sa_values.clear();
	}
};

// Helper functions for extracting fields from bam1_t and populating batches
namespace sam_utils {
	inline int64_t get_int_tag(const bam1_t *aln, const char *tag, int64_t default_val = -1) {
		uint8_t *aux = bam_aux_get(aln, tag);
		if (!aux) {
			return default_val;
		}
		return bam_aux2i(aux);
	}

	inline std::string get_string_tag(const bam1_t *aln, const char *tag) {
		const uint8_t *aux = bam_aux_get(aln, tag);
		if (!aux) {
			return "";
		}
		return std::string(bam_aux2Z(aux));
	}

	inline std::string cigar_to_string(const bam1_t *aln) {
		std::ostringstream oss;
		const uint32_t *cigar = reinterpret_cast<uint32_t *>(aln->data + aln->core.l_qname);
		for (uint32_t i = 0; i < aln->core.n_cigar; i++) {
			uint32_t op_len = bam_cigar_oplen(cigar[i]);
			char op_chr = bam_cigar_opchr(cigar[i]);
			oss << op_len << op_chr;
		}
		return oss.str();
	}

	// Parse a single record into a batch (append)
	inline void parse_record_to_batch(const bam1_t *aln, const sam_hdr_t *hdr, SAMRecordBatch &batch) {
		batch.read_ids.emplace_back(reinterpret_cast<char *>(aln->data));
		batch.flags.push_back(aln->core.flag);

		if (aln->core.tid >= 0) {
			batch.references.emplace_back(sam_hdr_tid2name(hdr, aln->core.tid));
		} else {
			batch.references.emplace_back("*");
		}

		batch.positions.push_back(aln->core.pos >= 0 ? aln->core.pos + 1 : 0);

		if (aln->core.flag & 0x4) { // BAM_FUNMAP
			batch.stop_positions.push_back(0);
		} else {
			hts_pos_t end_pos = bam_endpos(aln);
			batch.stop_positions.push_back(end_pos >= 0 ? end_pos + 1 : 0);
		}

		batch.mapqs.push_back(aln->core.qual);
		batch.cigars.emplace_back(cigar_to_string(aln));

		if (aln->core.mtid >= 0) {
			if (aln->core.mtid == aln->core.tid) {
				batch.mate_references.emplace_back("=");
			} else {
				batch.mate_references.emplace_back(sam_hdr_tid2name(hdr, aln->core.mtid));
			}
		} else {
			batch.mate_references.emplace_back("*");
		}

		batch.mate_positions.push_back(aln->core.mpos >= 0 ? aln->core.mpos + 1 : 0);
		batch.template_lengths.push_back(aln->core.isize);

		batch.tag_as_values.push_back(get_int_tag(aln, "AS"));
		batch.tag_xs_values.push_back(get_int_tag(aln, "XS"));
		batch.tag_ys_values.push_back(get_int_tag(aln, "YS"));
		batch.tag_xn_values.push_back(get_int_tag(aln, "XN"));
		batch.tag_xm_values.push_back(get_int_tag(aln, "XM"));
		batch.tag_xo_values.push_back(get_int_tag(aln, "XO"));
		batch.tag_xg_values.push_back(get_int_tag(aln, "XG"));
		batch.tag_nm_values.push_back(get_int_tag(aln, "NM"));
		batch.tag_yt_values.emplace_back(get_string_tag(aln, "YT"));
		batch.tag_md_values.emplace_back(get_string_tag(aln, "MD"));
		batch.tag_sa_values.emplace_back(get_string_tag(aln, "SA"));
	}
}
}; // namespace miint
