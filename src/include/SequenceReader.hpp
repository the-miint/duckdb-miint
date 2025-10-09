#pragma once
#include <optional>
#include <vector>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>
#include <thread>
#include "SequenceRecord.hpp"

namespace miint {
class SequenceReader {
public:
	explicit SequenceReader(const std::string &path1, const std::optional<std::string> &path2 = std::nullopt);
	std::vector<SequenceRecord> read(const int n) const;

private:
	using SeqStreamIn = klibpp::SeqStreamIn;

	std::unique_ptr<SeqStreamIn> reader1_;
	std::optional<std::unique_ptr<SeqStreamIn>> reader2_;

	bool paired_;

	std::vector<SequenceRecord> read_se(const int n) const;
	std::vector<SequenceRecord> read_pe(const int n) const;
};
}; // namespace miint
