#pragma once
#include <optional>
#include <vector>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>
#include "SequenceRecord.hpp"

namespace miint {
class SequenceReader {
public:
	explicit SequenceReader(const std::string &path1, const std::optional<std::string> &path2 = std::nullopt);
	SequenceRecordBatch read(const int n);

private:
	using SeqStreamIn = klibpp::SeqStreamIn;

	std::unique_ptr<SeqStreamIn> sequence1_reader_;
	std::optional<std::unique_ptr<SeqStreamIn>> sequence2_reader_;

	bool paired_;
	bool first_read_; // Track if we need to return buffered data
	std::vector<klibpp::KSeq> buffered_read1_;
	std::vector<klibpp::KSeq> buffered_read2_;

	SequenceRecordBatch read_se(const int n);
	SequenceRecordBatch read_pe(const int n);
};
}; // namespace miint
