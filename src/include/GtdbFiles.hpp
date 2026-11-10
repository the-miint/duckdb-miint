#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <vector>

// The file rules of the GTDB readers: release names, where GTDB serves a release's files
// and where the cache keeps them, which files a directory holds, and the header's column
// names. They use only DuckDB's header-only helpers, not its library, so their tests run
// in the `tests` binary, which links no DuckDB and is built for WASM too. Errors are
// std::runtime_error (std::invalid_argument for a release), with messages meant for the
// user; the readers add their own name and the path.
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

// Parses a release as a user writes it: a number, optionally a dot and another, each
// without a leading zero ("232", "214.1"). A bare major number is that release's first,
// ".0". Throws std::invalid_argument otherwise.
Release ParseRelease(std::string_view release);

// The files of a release the readers use.
enum class ReleaseFile { BAC120_METADATA, AR53_METADATA, QC_FAILED };

// The file's name as GTDB serves it, such as "bac120_metadata_r232.tsv.gz" or
// "qc_failed_r232.tsv". Only releases from 214.1 serve these names; earlier ones ship
// .tar.gz archives.
std::string ReleaseFileName(ReleaseFile file, const Release &release);

// Where GTDB serves the file, under its releases URL:
// release232/232.0/bac120_metadata_r232.tsv.gz, and qc_failed under auxillary_files/
// (GTDB's spelling). The reader joins it to the URL with DuckDB's FileSystem::JoinPath.
std::string ServedPath(ReleaseFile file, const Release &release);

// Where the cache keeps the file, under the cache directory: 232.0/<its name>. Every file
// of a release sits in one directory, so that directory is itself a source the readers
// accept.
std::string CachedPath(ReleaseFile file, const Release &release);

// The last component of a path or URL: what follows the last '/' or '\', and for a
// remote URL (any scheme but file://), without its query or fragment.
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
// empty name, a name that repeats another ignoring ASCII case, as DuckDB's column names
// do, or a name holding '"': the rows are read as CSV, where it quotes, and GTDB never
// writes one in a name.
std::vector<std::string> SplitHeader(std::string_view line);

} // namespace gtdb
} // namespace miint
