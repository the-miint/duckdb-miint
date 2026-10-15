#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace miint {

// Policy and message formatting for a run that fails: at open, after its one
// retry, or at completion (PerRunReader::Finish()), with an md5 digest mismatch
// at true EOF or a bad ascp exit code.
//
// ShouldSkipFailedRun returns true if the run should be
// skipped-with-warning (recorded in skipped_runs, surfaced via
// miint_warnings()) rather than thrown, which would abort the whole scan.
//
// Skip only when there is a sibling run to protect. A single-run scan keeps the
// hard throw: there is nothing else to discard, and a downstream that resolves a
// scalar accession one run per call (the Qiita ENA ingest job) relies on that
// error to fail the run. A multi-run scan -- a `varchar[]` of accessions, or a
// project accession that expands to many runs -- must NOT let one failed run
// abort the query and throw away every sibling already downloaded, so it skips
// the bad run and lets the rest land. (By true EOF the skipped run's rows were
// already emitted downstream; the warning flags them as unverified, the same
// contract as the mid-stream-truncation skip.) When every run is skipped, see
// EveryRunFailed.
inline bool ShouldSkipFailedRun(size_t total_runs) {
	return total_runs > 1;
}

// True when every run of a literal scan was skipped. There is then no sibling
// left to protect, so the scan raises BuildEveryRunFailedMessage instead of
// returning what the failed runs left: no rows, or one run's partial rows.
inline bool EveryRunFailed(size_t skipped_runs, size_t total_runs) {
	return skipped_runs > 0 && skipped_runs == total_runs;
}

// `accessions` in the order the runs failed; `first_error` is the first failure
// recorded, which with parallel workers need not be the first run's.
inline std::string BuildEveryRunFailedMessage(const std::vector<std::string> &accessions,
                                              const std::string &first_error) {
	std::string msg = "read_ena_sequences: every run failed (" + std::to_string(accessions.size()) + "):";
	for (const auto &acc : accessions) {
		msg += " " + acc;
	}
	return msg + "; first error: " + first_error;
}

// The mismatch message is a published contract (docs/insdc_ena.md, "Message
// contract"), pinned by test/cpp/test_ena_md5.cpp; rewording it is breaking.
inline std::string BuildMd5MismatchMessage(const std::string &label, const std::string &expected,
                                           const std::string &actual) {
	return "read_ena_sequences: md5 mismatch for '" + label + "': ENA reported " + expected +
	       " but downloaded bytes hash to " + actual;
}

} // namespace miint
