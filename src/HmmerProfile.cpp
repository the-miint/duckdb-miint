#include "HmmerProfile.hpp"

#include <cctype>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <memory>
#include <mutex>
#include <new>

extern "C" {
#include "hmmer.h"
}

namespace miint {
namespace hmmer {

namespace {

// The latest Easel exception's message on this thread. A fixed buffer, so the
// handler, which runs inside C frames, neither allocates nor can throw.
thread_local char exception_message[512] = "";

// Easel's exception handler. Records the message and returns, so the call that
// raised returns its error code. The source file and line, which Easel's
// default handler prints, are left out: they locate HMMER's code, not the
// user's input. So is errno, which Easel passes for failed system calls: the
// read path here makes none, as it reads from memory.
void RecordException(int, int, char *, int, char *format, va_list argp) {
	if (std::vsnprintf(exception_message, sizeof(exception_message), format, argp) < 0) {
		std::snprintf(exception_message, sizeof(exception_message), "%s",
		              "Easel exception with an unprintable message");
	}
}

struct HmmFileCloser {
	void operator()(P7_HMMFILE *hfp) const {
		p7_hmmfile_Close(hfp);
	}
};
struct AlphabetDestroyer {
	void operator()(ESL_ALPHABET *abc) const {
		esl_alphabet_Destroy(abc);
	}
};
struct ModelDestroyer {
	void operator()(P7_HMM *hmm) const {
		p7_hmm_Destroy(hmm);
	}
};

// Every format p7_hmmfile_OpenBuffer reads from text, each named by its tag.
constexpr int TEXT_FORMATS[] = {p7_HMMFILE_20, p7_HMMFILE_3a, p7_HMMFILE_3b, p7_HMMFILE_3c,
                                p7_HMMFILE_3d, p7_HMMFILE_3e, p7_HMMFILE_3f};

const char *FormatName(int format) {
	switch (format) {
	case p7_HMMFILE_20:
		return "HMMER2.0";
	case p7_HMMFILE_3a:
		return "HMMER3/a";
	case p7_HMMFILE_3b:
		return "HMMER3/b";
	case p7_HMMFILE_3c:
		return "HMMER3/c";
	case p7_HMMFILE_3d:
		return "HMMER3/d";
	case p7_HMMFILE_3e:
		return "HMMER3/e";
	case p7_HMMFILE_3f:
		return "HMMER3/f";
	default:
		// p7_hmmfile_OpenBuffer sets one of the above or fails.
		throw std::logic_error("HMMER profile: unexpected format code " + std::to_string(format));
	}
}

// p7_hmmfile_OpenBuffer failed with a format error, which it reports with no
// message and no line: its first token is not a HMMER format tag, or there is
// none. Names the token, found as Easel's parser finds it: the first word of
// the first line that is neither blank nor a comment.
[[noreturn]] void ThrowNotAModel(std::string_view record) {
	int64_t line = 0;
	size_t start = 0;
	while (start < record.size()) {
		size_t end = record.find('\n', start);
		if (end == std::string_view::npos) {
			end = record.size();
		}
		line++;
		const auto word = FirstWord(record.substr(start, end - start));
		if (!IsBlankOrComment(word)) {
			throw ProfileParseError("Format tag is '" + std::string(word) +
			                            "': unrecognized (expected HMMER3/a to HMMER3/f, or HMMER2.0)",
			                        line);
		}
		start = end + 1;
	}
	// No line to point at: the record's first is reported.
	throw ProfileParseError("Empty record: no HMMER format tag", 1);
}

} // namespace

std::string_view FirstWord(std::string_view line) {
	size_t begin = 0;
	while (begin < line.size() && std::isspace(static_cast<unsigned char>(line[begin]))) {
		begin++;
	}
	size_t end = begin;
	while (end < line.size() && !std::isspace(static_cast<unsigned char>(line[end]))) {
		end++;
	}
	return line.substr(begin, end - begin);
}

bool IsBlankOrComment(std::string_view first_word) {
	return first_word.empty() || first_word[0] == '#';
}

bool IsFormatTag(std::string_view first_word) {
	for (int format : TEXT_FORMATS) {
		if (first_word == FormatName(format)) {
			return true;
		}
	}
	return false;
}

bool IsModelEnd(std::string_view first_word) {
	return first_word == "//";
}

void Initialize() {
	static std::once_flag initialized;
	std::call_once(initialized, [] {
		esl_exception_SetHandler(RecordException);
		p7_FLogsumInit();
	});
}

std::string TakeExceptionMessage() {
	std::string message(exception_message);
	exception_message[0] = '\0';
	return message;
}

ProfileHeader ParseProfile(std::string_view record) {
	// Anything an earlier call on this thread left must not be taken as this
	// read's.
	TakeExceptionMessage();

	static_assert(MAX_RECORD_BYTES == static_cast<size_t>(std::numeric_limits<int>::max()),
	              "HMMER's buffer length is an int");
	if (record.size() > MAX_RECORD_BYTES) {
		throw ProfileParseError(
		    "Model record of " + std::to_string(record.size()) + " bytes is larger than HMMER can read (2 GiB)", 1);
	}
	// HMMER parses the buffer in place: it is not copied, and must outlive the
	// handle.
	P7_HMMFILE *opened = nullptr;
	int status = p7_hmmfile_OpenBuffer(record.empty() ? "" : record.data(), static_cast<int>(record.size()), &opened);
	if (status == eslEMEM) {
		throw std::bad_alloc();
	}
	if (status != eslOK) {
		ThrowNotAModel(record);
	}
	std::unique_ptr<P7_HMMFILE, HmmFileCloser> hfp(opened);

	// One read per handle, and so one model per record. A second read from a
	// buffer handle crashes in HMMER 3.4 (ftello on its NULL FILE *).
	//
	// A fresh alphabet each time: HMMER creates it from the model's ALPH line,
	// so a file may mix alphabets. Declared before the model, which borrows it
	// (hmm->abc), so the model is destroyed first.
	ESL_ALPHABET *abc_read = nullptr;
	P7_HMM *hmm_read = nullptr;
	status = p7_hmmfile_Read(hfp.get(), &abc_read, &hmm_read);
	std::unique_ptr<ESL_ALPHABET, AlphabetDestroyer> abc(abc_read);
	std::unique_ptr<P7_HMM, ModelDestroyer> hmm(hmm_read);
	detail::ThrowIfReadFailed(status, hmm != nullptr, hfp->errbuf, hfp->efp->linenumber);

	// HMMER sets the consensus itself only for formats older than HMMER3/e; a
	// later file can say CONS no.
	if (!(hmm->flags & p7H_CONS) && p7_hmm_SetConsensus(hmm.get(), nullptr) != eslOK) {
		throw std::bad_alloc(); // its only failure
	}

	ProfileHeader h;
	h.name = hmm->name; // HMMER rejects a model without one
	if ((hmm->flags & p7H_ACC) && hmm->acc) {
		h.accession = hmm->acc;
	}
	if ((hmm->flags & p7H_DESC) && hmm->desc) {
		h.description = hmm->desc;
	}
	h.length = hmm->M;
	h.alphabet = esl_abc_DecodeType(hmm->abc->type);
	h.format = FormatName(hfp->format);
	// No flag marks these three: they are present when HMMER would write them
	// back out (p7_hmmfile_WriteASCII).
	if (hmm->max_length > 0) {
		h.max_length = hmm->max_length;
	}
	if (hmm->nseq > 0) {
		h.sequence_count = hmm->nseq;
	}
	if (hmm->eff_nseq >= 0) {
		h.effective_sequence_count = hmm->eff_nseq;
	}
	if (hmm->flags & p7H_CHKSUM) {
		h.checksum = hmm->checksum;
	}
	if (hmm->flags & p7H_GA) {
		h.gathering_sequence = hmm->cutoff[p7_GA1];
		h.gathering_domain = hmm->cutoff[p7_GA2];
	}
	if (hmm->flags & p7H_TC) {
		h.trusted_sequence = hmm->cutoff[p7_TC1];
		h.trusted_domain = hmm->cutoff[p7_TC2];
	}
	if (hmm->flags & p7H_NC) {
		h.noise_sequence = hmm->cutoff[p7_NC1];
		h.noise_domain = hmm->cutoff[p7_NC2];
	}
	h.consensus = hmm->consensus + 1; // [0] is a blank; [1..M] are the nodes
	if (hmm->flags & p7H_STATS) {
		h.msv_mu = hmm->evparam[p7_MMU];
		h.msv_lambda = hmm->evparam[p7_MLAMBDA];
		h.viterbi_mu = hmm->evparam[p7_VMU];
		h.viterbi_lambda = hmm->evparam[p7_VLAMBDA];
		h.forward_tau = hmm->evparam[p7_FTAU];
		h.forward_lambda = hmm->evparam[p7_FLAMBDA];
	}
	return h;
}

namespace detail {
void ThrowIfReadFailed(int status, bool have_model, const char *errbuf, int64_t line) {
	// Taken whatever the outcome, so it cannot be blamed on a later read.
	std::string exception = TakeExceptionMessage();
	if (status == eslOK) {
		if (have_model) {
			return;
		}
		throw ProfileParseError("HMMER reported success but returned no model", line);
	}
	if (status == eslEMEM) {
		throw std::bad_alloc();
	}
	// The parser's own message says where in the model it stopped; an
	// exception's comes from a function the parser called.
	if (errbuf && *errbuf) {
		throw ProfileParseError(errbuf, line);
	}
	// The record ended before the header did. HMMER takes that as the end of
	// the file, with no message.
	if (status == eslEOF) {
		throw ProfileParseError("Premature end of data in the header: no HMM line", line);
	}
	if (!exception.empty()) {
		throw ProfileParseError(exception, line);
	}
	throw ProfileParseError("HMMER failed with status " + std::to_string(status) + " and no message", line);
}
} // namespace detail

} // namespace hmmer
} // namespace miint
