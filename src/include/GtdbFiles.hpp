#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// The file rules of the GTDB readers, apart from DuckDB: release names, where GTDB
// serves a release's files and where the cache keeps them, which files a directory
// holds, and the first line of a file as its bytes arrive. Errors are std::runtime_error
// (std::invalid_argument for a release), with messages meant for the user; the readers
// add their own name and the path.
namespace miint {
namespace gtdb {

// A GTDB release as its directories name it: "232.0", "214.1".
struct Release {
	std::string major;
	std::string minor;

	std::string ToString() const {
		return major + "." + minor;
	}
};

// Parses a release as a user writes it: digits, optionally a dot and more digits. A bare
// major number is that release's first, ".0". Throws std::invalid_argument otherwise.
Release ParseRelease(std::string_view release);

// The files of a release the readers use.
enum class ReleaseFile { BAC120_METADATA, AR53_METADATA, QC_FAILED };

// The file's name as GTDB serves it, such as "bac120_metadata_r232.tsv.gz" or
// "qc_failed_r232.tsv". Only releases from 214.1 serve these names; earlier ones ship
// .tar.gz archives.
std::string ReleaseFileName(ReleaseFile file, const Release &release);

// Where GTDB serves the file under `base_url`:
// <base>/release232/232.0/bac120_metadata_r232.tsv.gz, and qc_failed under
// auxillary_files/ (GTDB's spelling).
std::string ReleaseFileUrl(std::string_view base_url, ReleaseFile file, const Release &release);

// Where the cache keeps the file: <cache_dir>/232.0/<its name>. Every file of a release
// sits in one directory, so that directory is itself a source the readers accept.
std::string ReleaseCachePath(std::string_view cache_dir, ReleaseFile file, const Release &release);

// The last component of a path or URL: what follows the last '/' or '\', and for an
// http(s) URL, without its query or fragment.
std::string_view FileName(std::string_view path);

// The domain a metadata file holds, from its name: "Bacteria" for bac120_*, "Archaea"
// for ar53_* and for ar122_*, the archaeal marker set before release 207. nullopt for
// any other name.
std::optional<std::string> MetadataDomain(std::string_view file_name);

// The metadata files among a directory's file names: at most one bacterial
// (bac120_metadata_r<N>.tsv or .tsv.gz) and one archaeal (ar53_ or ar122_), at least one
// of them, of one release. Bacterial first. Throws when the names do not fit.
std::vector<std::string> MatchMetadataFiles(const std::vector<std::string> &names);

// The one qc_failed_r<N>.tsv or .tsv.gz among a directory's file names. Throws when there
// is none or more than one.
std::string MatchQcFailedFile(const std::vector<std::string> &names);

// Splits a header line on tabs into column names, exactly as written. Throws for an
// empty name, or a name that repeats another ignoring ASCII case, as DuckDB's column
// names do.
std::vector<std::string> SplitHeader(std::string_view line);

// Reads a file's first line from its bytes as they arrive, so a caller can stop reading
// (or downloading) once it has the header. Plain text and gzip are both read, told apart
// by gzip's magic bytes.
class HeaderLineReader {
public:
	static constexpr size_t DEFAULT_MAX_LINE = size_t(1) << 20;

	explicit HeaderLineReader(size_t max_line = DEFAULT_MAX_LINE);
	~HeaderLineReader();
	HeaderLineReader(HeaderLineReader &&) noexcept;
	HeaderLineReader &operator=(HeaderLineReader &&) noexcept;

	// Hands over the next bytes. Returns true once the first line is complete; bytes after
	// that are ignored. Throws std::runtime_error for a line longer than max_line, a corrupt
	// gzip stream, or a gzip stream that ends before its first newline.
	bool Feed(const char *data, size_t size);

	// Ends the input and returns the first line without its "\n" or "\r\n". Throws
	// std::runtime_error when the input ended before the line did: no newline, or a gzip
	// stream cut short.
	std::string Finish();

	// Whether the bytes were gzip. Known once Feed returned true or Finish returned.
	bool IsGzip() const {
		return gzip_;
	}

private:
	struct Inflater;

	bool Decide();
	bool Append(const char *data, size_t size);
	bool Inflate(const char *data, size_t size);

	size_t max_line_;
	std::string pending_; // the first bytes, until there are enough to tell gzip apart
	bool decided_ = false;
	bool gzip_ = false;
	bool done_ = false;
	std::string line_;
	std::unique_ptr<Inflater> inflater_;
};

} // namespace gtdb
} // namespace miint
