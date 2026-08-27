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

// A local-file sequence stream that checks gzopen() actually succeeded.
//
// klibpp::SeqStreamIn hands gzopen()'s result straight to KStream without looking at it.
// A NULL handle does not fail loudly: it produces a stream that yields zero records, which
// SequenceReader cannot tell apart from a genuinely empty file, so every failed open --
// EACCES, EMFILE/ENFILE on descriptor exhaustion, ENOMEM -- was reported as
// "Empty file: <path>" naming a file that is fine. ext/kseq++ is vendored unpatched, so the
// check lives here instead. Deriving from SeqStreamIn's own base lets us open the file
// ourselves and read errno with nothing running in between; SeqStreamIn adds only
// constructors, so behaviour is otherwise identical.
class CheckedSeqStreamIn : public klibpp::SeqStreamIn::base_type {
public:
	explicit CheckedSeqStreamIn(const std::string &path)
	    : klibpp::SeqStreamIn::base_type(OpenOrThrow(path), gzread, gzclose) {
	}

private:
	static gzFile OpenOrThrow(const std::string &path);
};

class SequenceReader {
public:
	explicit SequenceReader(const std::string &path1, const std::optional<std::string> &path2 = std::nullopt,
	                        bool allow_format_mismatch = false);

#ifdef MIINT_STATIC_BUILD
	// Streaming constructor for remote files via DuckDB FileHandle. Ownership of
	// the streams is transferred via DuckDBSeqStreamHandle: pass a nullptr-valued
	// handle (still constructed with the duckdb_seq_close deleter) for the
	// second slot when the run is single-end.
	SequenceReader(DuckDBSeqStreamHandle stream1, DuckDBSeqStreamHandle stream2_or_null,
	               bool allow_format_mismatch = false);

#if MIINT_ASPERA_SUPPORTED
	// Streaming constructor for Aspera pipe-backed streams.
	SequenceReader(AsperaSeqStreamHandle stream1, AsperaSeqStreamHandle stream2_or_null,
	               bool allow_format_mismatch = false);
#endif // MIINT_ASPERA_SUPPORTED
#endif // MIINT_STATIC_BUILD

	// Read up to `n` records or until cumulative bytes across the accumulated batch exceed
	// `max_bytes`, whichever comes first. Byte accounting sums `name + comment + seq + qual`
	// across both mates in paired-end mode. At least one record is always returned if the
	// stream has data, even when a single record alone exceeds `max_bytes` (starvation guard).
	// Defaults to SIZE_MAX, preserving the row-only behavior for existing callers.
	SequenceRecordBatch read(const int n, const size_t max_bytes = SIZE_MAX);

	// Test seam: the largest number of records materialized by any single underlying stream
	// poll over this reader's lifetime. Lets tests assert that large records do not trigger
	// oversized prefetches (see dynamic poll sizing in read_se/read_pe).
	size_t MaxPollCount() const {
		return max_poll_count_;
	}

private:
	using SeqStreamIn = CheckedSeqStreamIn;

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

	// Running estimate of bytes per record (per pair in PE), used to size each poll so a
	// single read_stream call never materializes far more than the byte budget. 0 until the
	// first record is accepted.
	size_t observed_record_bytes_ = 0;
	size_t max_poll_count_ = 0; // test seam; see MaxPollCount()

	SequenceRecordBatch read_se(const int n, const size_t max_bytes);
	SequenceRecordBatch read_pe(const int n, const size_t max_bytes);
};
}; // namespace miint
