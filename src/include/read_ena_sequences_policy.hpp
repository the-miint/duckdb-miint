#pragma once

#include <cstddef>

namespace miint {

// Policy for a run's integrity failure detected at completion (PerRunReader::
// Finish()): an md5 digest mismatch at true EOF, or a bad ascp exit code.
// Returns true if the run should be skipped-with-warning (recorded in
// skipped_runs, surfaced via miint_warnings()) rather than thrown, which would
// abort the whole scan.
//
// Skip only when there is a sibling run to protect. A single-run scan keeps the
// hard throw: there is nothing else to discard, and a downstream that resolves a
// scalar accession one run per call (the Qiita ENA ingest job) relies on that
// error to fail the run. A multi-run scan -- a `varchar[]` of accessions, or a
// project accession that expands to many runs -- must NOT let one corrupt run
// abort the query and throw away every sibling already downloaded, so it skips
// the bad run and lets the rest land. (By true EOF the skipped run's rows were
// already emitted downstream; the warning flags them as unverified, the same
// contract as the mid-stream-truncation skip.)
inline bool ShouldSkipRunIntegrityFailure(size_t total_runs) {
	return total_runs > 1;
}

} // namespace miint
