#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

// HMMER profile HMMs (HMMER 3.4, Easel 0.49), read one model at a time. Nothing
// here exposes a HMMER or Easel type: code that includes only this header needs
// neither library's headers nor their include path.
namespace miint {
namespace hmmer {

// Prepares HMMER for use in this process. Call it once at extension load,
// before anything can call into HMMER; later calls do nothing, from any thread.
//
// - Installs Easel's exception handler. It is one process-global pointer, and
//   the default handler aborts the process. This one records the message for
//   the thread that raised it (TakeExceptionMessage) and returns, so the HMMER
//   call that raised returns its error code instead.
// - Fills HMMER's logsum lookup table. HMMER otherwise fills it on first use,
//   and two threads getting there at once race on it.
void Initialize();

// The message of the latest Easel exception raised on this thread, or "" if
// none was raised since the last call. Clears it.
std::string TakeExceptionMessage();

// The header of one model, which is one row of read_hmm. A field the model does
// not have is nullopt: HMMER's own flags decide, never its placeholder values.
struct ProfileHeader {
	std::string name;                              // NAME
	std::optional<std::string> accession;          // ACC
	std::optional<std::string> description;        // DESC
	int64_t length = 0;                            // LENG: the number of match states
	std::string alphabet;                          // ALPH, under Easel's names: amino, DNA, RNA, ...
	std::string format;                            // HMMER3/a ... HMMER3/f, or HMMER2.0
	std::optional<int64_t> max_length;             // MAXL
	std::optional<int64_t> sequence_count;         // NSEQ
	std::optional<float> effective_sequence_count; // EFFN
	std::optional<uint32_t> checksum;              // CKSUM
	std::optional<float> gathering_sequence;       // GA, per-sequence
	std::optional<float> gathering_domain;         // GA, per-domain
	std::optional<float> trusted_sequence;         // TC, per-sequence
	std::optional<float> trusted_domain;           // TC, per-domain
	std::optional<float> noise_sequence;           // NC, per-sequence
	std::optional<float> noise_domain;             // NC, per-domain
	std::string consensus;                         // the CONS column, or computed by HMMER when the model has none
	std::optional<float> msv_mu;                   // STATS LOCAL MSV
	std::optional<float> msv_lambda;               // STATS LOCAL MSV
	std::optional<float> viterbi_mu;               // STATS LOCAL VITERBI
	std::optional<float> viterbi_lambda;           // STATS LOCAL VITERBI
	std::optional<float> forward_tau;              // STATS LOCAL FORWARD
	std::optional<float> forward_lambda;           // STATS LOCAL FORWARD
};

// A model HMMER could not read. what() is HMMER's message, or this layer's
// where HMMER gives none. Line() is the line of the record, counting from 1,
// at which the reader stopped.
class ProfileParseError : public std::runtime_error {
public:
	ProfileParseError(const std::string &message, int64_t line) : std::runtime_error(message), line(line) {
	}
	int64_t Line() const {
		return line;
	}

private:
	int64_t line;
};

// How HMMER reads the lines of a model file, for a caller splitting a file into
// records. HMMER takes each line's first word, its first whitespace-delimited
// token: it skips a line with none (blank) or whose first word starts with '#'
// (a comment); a model's first line it doesn't skip starts with the model's
// format tag; and a model ends at a line whose first word is "//".
std::string_view FirstWord(std::string_view line);
bool IsBlankOrComment(std::string_view first_word);
bool IsFormatTag(std::string_view first_word); // HMMER3/a ... HMMER3/f, or HMMER2.0
bool IsModelEnd(std::string_view first_word);

// The largest record ParseProfile reads: HMMER takes a buffer's length as an int.
constexpr size_t MAX_RECORD_BYTES = 2147483647;

// Reads one model: the text of one record, from its format line ("HMMER3/f
// ...", "HMMER2.0 ...") through its closing "//". Text after the "//" is not
// read; splitting a file into records is the caller's job.
//
// For a HMMER2 model, HMMER calibrates the E-value statistics as it reads, as
// the HMMER programs do. Its random seed is fixed, so a model always gets the
// same values.
//
// Throws ProfileParseError when HMMER rejects the model, and std::bad_alloc
// when HMMER runs out of memory. Initialize() must have run.
ProfileHeader ParseProfile(std::string_view record);

namespace detail {
// What ParseProfile makes of p7_hmmfile_Read's result: returns when the read
// gave a model and throws as ParseProfile documents otherwise. <status> is the
// HMMER return code, <errbuf> HMMER's error buffer, and <line> the line the
// reader stopped at. Takes this thread's exception message.
//
// Exposed for the tests: one outcome, a success code with no model (HMMER
// returns it when allocating a HMMER2 model's background fails), can't be
// produced from any input.
void ThrowIfReadFailed(int status, bool have_model, const char *errbuf, int64_t line);
} // namespace detail

} // namespace hmmer
} // namespace miint
