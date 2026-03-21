#pragma once
#include <optional>
#include <variant>
#include <vector>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>
#include "SequenceRecord.hpp"

#ifdef MIINT_STATIC_BUILD
#include "duckdb_seq_stream.hpp"
#endif

namespace miint {
class SequenceReader {
public:
	explicit SequenceReader(const std::string &path1, const std::optional<std::string> &path2 = std::nullopt);

#ifdef MIINT_STATIC_BUILD
	// Streaming constructor for remote files via DuckDB FileHandle.
	// Takes ownership of stream pointers (freed by kseq++ close callback).
	SequenceReader(DuckDBSeqStream *stream1, DuckDBSeqStream *stream2_or_null);
#endif

	SequenceRecordBatch read(const int n);

private:
	using SeqStreamIn = klibpp::SeqStreamIn;

#ifdef MIINT_STATIC_BUILD
	using StreamVar = std::variant<std::unique_ptr<SeqStreamIn>, std::unique_ptr<DuckDBSeqStreamIn>>;
#else
	using StreamVar = std::variant<std::unique_ptr<SeqStreamIn>>;
#endif

	static std::vector<klibpp::KSeq> read_stream(StreamVar &var, int n);

	StreamVar sequence1_reader_;
	std::optional<StreamVar> sequence2_reader_;

	bool paired_;
	bool first_read_; // Track if we need to return buffered data
	std::vector<klibpp::KSeq> buffered_read1_;
	std::vector<klibpp::KSeq> buffered_read2_;

	SequenceRecordBatch read_se(const int n);
	SequenceRecordBatch read_pe(const int n);
};
}; // namespace miint
