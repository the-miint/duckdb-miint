#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct BlastSubmitResult {
	std::string rid;
	int rtoe = 0;
	std::string error;
};

enum class BlastStatus { WAITING, READY, UNKNOWN };

struct BlastHit {
	std::string query_id;
	std::string subject_id;
	double pct_identity;
	int32_t alignment_length;
	int32_t mismatches;
	int32_t gap_opens;
	int64_t query_start;
	int64_t query_end;
	int64_t subject_start;
	int64_t subject_end;
	double evalue;
	double bit_score;
};

class BlastParser {
public:
	static BlastSubmitResult ParseSubmitResponse(const std::string &html);
	static BlastStatus ParseStatusResponse(const std::string &text);
	static std::vector<BlastHit> ParseTabularResults(const std::string &tabular);
	static bool ValidateProgram(const std::string &program);
	static std::string BuildFastaPayload(const std::vector<std::string> &ids,
	                                     const std::vector<std::string> &sequences);
	static std::string UrlEncode(const std::string &value);
	static std::string BuildSubmitBody(const std::string &program, const std::string &database,
	                                   const std::string &fasta_query, double evalue, int max_targets, bool megablast);

	struct SequenceBatch {
		std::vector<std::string> ids;
		std::vector<std::string> sequences;
	};
	static std::vector<SequenceBatch> SplitIntoBatches(const std::vector<std::string> &ids,
	                                                   const std::vector<std::string> &sequences, size_t max_bytes);
};

} // namespace miint
