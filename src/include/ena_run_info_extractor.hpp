#pragma once

#include "aspera_utils.hpp"
#include "ena_parser.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace miint {

enum class ENASequenceFormat { FASTX, SFF };

struct ENARunInfo {
	std::string run_accession;
	std::string sample_accession;
	std::string experiment_accession;
	std::vector<std::string> fastq_urls;
	// ENA's per-file md5, aligned 1:1 with fastq_urls (same index i = same
	// file). Empty when ENA did not report md5s for this run, or when the
	// reported count didn't line up with fastq_urls after paired-file
	// filtering -- see FromTSVRow. Never partially/misaligned: either fully
	// populated and aligned, or fully empty.
	std::vector<std::string> fastq_md5;
	std::vector<AsperaPath> aspera_paths;
	bool is_paired = false;
	uint64_t total_bytes = 0;
	ENASequenceFormat format = ENASequenceFormat::FASTX;
	std::string sff_url;
};

struct ENAColumnIndexMap {
	int run_accession = -1;
	int sample_accession = -1;
	int experiment_accession = -1;
	int study_accession = -1;
	int fastq_ftp = -1;
	int fastq_aspera = -1;
	int fastq_bytes = -1;
	int fastq_md5 = -1;
	int library_layout = -1;
	int submitted_ftp = -1;
	int submitted_aspera = -1;
	int submitted_bytes = -1;
	int submitted_format = -1;

	static ENAColumnIndexMap FromHeader(const std::vector<std::string> &header);

	static std::string Get(const std::vector<std::string> &row, int col);
};

class ENARunInfoExtractor {
public:
	// Pure per-row extractor. Returns zero or more RunInfo records from a single TSV row.
	//   FASTQ runs yield exactly one RunInfo.
	//   SFF runs flatten to one RunInfo per SFF file.
	//   Rows with no usable format return an empty vector.
	//
	// prefer_format ∈ {"auto", "fastq", "sff"}:
	//   "auto"  — FASTQ if available, else SFF
	//   "fastq" — FASTQ only
	//   "sff"   — SFF if available, else fall back to FASTQ
	static std::vector<ENARunInfo> FromTSVRow(const std::vector<std::string> &row, const ENAColumnIndexMap &cols,
	                                          const std::string &prefer_format);
};

} // namespace miint
