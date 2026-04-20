#pragma once
#include <climits>
#include <cstddef>
#include <optional>
#include <variant>
#include <vector>
#include <string>
#include <memory>
#include <kseq++/seqio.hpp>
#include "SequenceRecord.hpp"

#ifdef MIINT_STATIC_BUILD
#include "duckdb_seq_stream.hpp"
#include "aspera_stream.hpp"
#endif

namespace miint {
class SequenceReader {
public:
	explicit SequenceReader(const std::string &path1, const std::optional<std::string> &path2 = std::nullopt,
	                        bool allow_format_mismatch = false);

#ifdef MIINT_STATIC_BUILD
	// Streaming constructor for remote files via DuckDB FileHandle.
	// Takes ownership of stream pointers (freed by kseq++ close callback).
	SequenceReader(DuckDBSeqStream *stream1, DuckDBSeqStream *stream2_or_null, bool allow_format_mismatch = false);

#if MIINT_ASPERA_SUPPORTED
	// Streaming constructor for Aspera pipe-backed streams.
	SequenceReader(AsperaSeqStream *stream1, AsperaSeqStream *stream2_or_null, bool allow_format_mismatch = false);

	// Mixed constructor: DuckDB stream for s1 (e.g., temp file), Aspera stream for s2 (live pipe).
	// Used for paired-end Aspera downloads where R1 is buffered to temp file.
	SequenceReader(DuckDBSeqStream *stream1, AsperaSeqStream *stream2, bool allow_format_mismatch = false);
#endif // MIINT_ASPERA_SUPPORTED
#endif // MIINT_STATIC_BUILD

	// Read up to `n` records or until cumulative bytes across the accumulated batch exceed
	// `max_bytes`, whichever comes first. Byte accounting sums `name + comment + seq + qual`
	// across both mates in paired-end mode. At least one record is always returned if the
	// stream has data, even when a single record alone exceeds `max_bytes` (starvation guard).
	// Defaults to SIZE_MAX, preserving the row-only behavior for existing callers.
	SequenceRecordBatch read(const int n, const size_t max_bytes = SIZE_MAX);

private:
	using SeqStreamIn = klibpp::SeqStreamIn;

#if defined(MIINT_STATIC_BUILD) && MIINT_ASPERA_SUPPORTED
	using StreamVar = std::variant<std::unique_ptr<SeqStreamIn>, std::unique_ptr<DuckDBSeqStreamIn>,
	                               std::unique_ptr<AsperaSeqStreamIn>>;
#elif defined(MIINT_STATIC_BUILD)
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

	SequenceRecordBatch read_se(const int n, const size_t max_bytes);
	SequenceRecordBatch read_pe(const int n, const size_t max_bytes);
};
}; // namespace miint
