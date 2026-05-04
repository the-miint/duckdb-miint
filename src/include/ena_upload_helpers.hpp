// SPDX-License-Identifier: MIT
//
// Pure-data helpers for the `ena_upload_reads` table function:
//
//   - Per-sample layout detection (single / paired / paired-interleaved /
//     mixed-error)
//   - Filename emission per layout
//   - Upload target URL parsing (aspera:// vs file://)
//   - `ascp --mode=send` argv construction
//
// All routines here are DuckDB-free — the unit-test target links the .cpp
// directly without pulling libduckdb. Same convention as `fastq_encoder`.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// User-requested or auto-detected per-sample FASTQ layout. AUTO is a
// bind-time placeholder; ResolveLayout collapses it to SINGLE or PAIRED based
// on the actual rows. PAIRED_INTERLEAVED is a user override (the only layout
// the encoder side cannot infer from data alone) that emits one file with R1
// and R2 records alternating.
enum class FastqLayoutMode : uint8_t {
	AUTO,
	SINGLE,
	PAIRED,
	PAIRED_INTERLEAVED,
};

// Lower-case canonical names for error messages and option parsing
// (`auto`, `single`, `paired`, `paired_interleaved`).
const char *FastqLayoutModeName(FastqLayoutMode mode);

// Parse a user-supplied layout string (case-insensitive). Throws
// std::runtime_error with the offending value when the name isn't recognised.
FastqLayoutMode ParseFastqLayoutMode(const std::string &name);

// Resolve the per-sample layout from rows. `has_r2[i]` records whether row i
// in the sample's grouped chunk supplied a non-null R2. Mixed rows (some with,
// some without) raise std::runtime_error naming the sample_ref and the
// 0-based index of the first inconsistent row. PAIRED_INTERLEAVED is honored
// only when every row has a non-null R2; otherwise behaves like PAIRED.
FastqLayoutMode ResolveLayout(const std::string &sample_ref, FastqLayoutMode requested,
                              const std::vector<bool> &has_r2);

// Output FASTQ filenames for one sample given its resolved layout.
//   SINGLE             → ["{sample_ref}.fastq.gz"]
//   PAIRED             → ["{sample_ref}_1.fastq.gz", "{sample_ref}_2.fastq.gz"]
//   PAIRED_INTERLEAVED → ["{sample_ref}.fastq.gz"]
// Passing AUTO throws std::runtime_error — the caller must resolve first.
std::vector<std::string> OutputFilenames(const std::string &sample_ref, FastqLayoutMode layout);

// ---- Upload target URL parsing -------------------------------------------

enum class UploadTransport : uint8_t {
	ASPERA,     // aspera://host/path/  → ascp --mode=send (temp file)
	LOCAL_FILE, // file://path/         → direct write
	CURL,       // ftp/ftps/http/https://host/path/ → libcurl streaming
};

// The libcurl transport spans four URL schemes. ScheamName() echoes the
// original wire scheme so error messages can be specific.
struct UploadTargetURL {
	UploadTransport transport;
	std::string scheme;       // lower-case original scheme ("ftp", "https", "aspera", "file", ...)
	std::string host;         // empty for LOCAL_FILE
	std::string remote_dir;   // always ends with '/'; "/" if no path component
	std::string url_for_curl; // populated only for CURL transport: full URL with trailing slash
};

// Parse one of:
//   aspera://host[/path/]               → ASPERA
//   file://[abs-or-relative-path]/      → LOCAL_FILE
//   ftp://host[/path/]   ftps://...     → CURL
//   http://host[/path/]  https://...    → CURL
// The remote_dir is normalised to always have a trailing '/'. Unsupported
// schemes or malformed input throw std::runtime_error.
UploadTargetURL ParseUploadTargetURL(const std::string &url);

// ---- ascp argv builder ---------------------------------------------------

struct AsperaSendOptions {
	std::string ascp_path;  // resolved binary path
	std::string key_path;   // resolved openssh private key path
	std::string user;       // Webin-XXXXX
	std::string host;       // e.g. "webin2.ebi.ac.uk"
	int port = 33001;       // ENA Aspera default
	std::string max_rate;   // e.g. "300m"; empty omits the -l flag
	std::string local_path; // already-gzipped file to upload
	std::string remote_dir; // destination directory under the user's webin area
};

// Build argv for `ascp --mode=send`. Layout:
//   [ascp_path, --mode=send, --user=USER, --host=HOST, -P, port,
//    -i, key_path, -Q, -d, (-l, max_rate)?, local_path, remote_dir]
//
// Notes:
//   - The ENA password is NOT in argv — it goes via ASPERA_SCP_PASS in the
//     child process environment so it never lands in /proc/self/cmdline.
//   - `-d` tells ascp to create the destination directory if missing.
//   - `-Q` enables adaptive rate control; `-l` caps it when supplied.
std::vector<std::string> BuildAscpSendArgv(const AsperaSendOptions &opts);

} // namespace duckdb
