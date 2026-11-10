#include "GtdbFiles.hpp"

#include "remote_file_helper.hpp"
#include "duckdb/common/string_util.hpp"

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <unordered_map>

namespace miint {
namespace gtdb {

namespace {

// ASCII digits, by DuckDB's test (std::isdigit is undefined for a negative char).
bool IsDigits(std::string_view s) {
	return !s.empty() && std::all_of(s.begin(), s.end(), duckdb::StringUtil::CharacterIsDigit);
}

// A number as GTDB writes one in a release's name: digits, without a leading zero.
bool IsReleaseNumber(std::string_view s) {
	return IsDigits(s) && (s.size() == 1 || s.front() != '0');
}

// GTDB's name for a file of a release: <prefix>_r<major>.tsv[.gz].
struct NamedFile {
	std::string_view prefix;
	std::string_view release; // the major number in the name
};

// Reads a name of the form <prefix>_r<digits>.tsv or .tsv.gz, for one of `prefixes`.
std::optional<NamedFile> ParseReleaseFileName(std::string_view name, std::initializer_list<std::string_view> prefixes) {
	for (const auto prefix : prefixes) {
		if (!name.starts_with(prefix) || name.substr(prefix.size(), 2) != "_r") {
			continue;
		}
		auto rest = name.substr(prefix.size() + 2);
		for (const std::string_view extension : {".tsv", ".tsv.gz"}) {
			if (rest.ends_with(extension) && IsDigits(rest.substr(0, rest.size() - extension.size()))) {
				return NamedFile {prefix, rest.substr(0, rest.size() - extension.size())};
			}
		}
	}
	return std::nullopt;
}

std::string JoinNames(const std::vector<std::string> &names) {
	std::string out;
	for (const auto &name : names) {
		out += (out.empty() ? "" : ", ") + name;
	}
	return out;
}

// At most one file per name group: two would read the same genomes twice.
void RequireAtMostOne(const std::vector<std::string> &found, const char *what) {
	if (found.size() > 1) {
		throw std::runtime_error(std::string("more than one ") + what + " file: " + JoinNames(found));
	}
}

} // namespace

Release ParseRelease(std::string_view release) {
	const auto dot = release.find('.');
	const auto major = release.substr(0, dot);
	const auto minor = dot == std::string_view::npos ? std::string_view("0") : release.substr(dot + 1);
	if (!IsReleaseNumber(major) || !IsReleaseNumber(minor)) {
		throw std::invalid_argument("release must be a GTDB release number such as '232' or '214.1', not '" +
		                            std::string(release) + "'");
	}
	return Release {std::string(major), std::string(minor)};
}

std::string ReleaseFileName(ReleaseFile file, const Release &release) {
	switch (file) {
	case ReleaseFile::BAC120_METADATA:
		return "bac120_metadata_r" + release.major + ".tsv.gz";
	case ReleaseFile::AR53_METADATA:
		return "ar53_metadata_r" + release.major + ".tsv.gz";
	case ReleaseFile::QC_FAILED:
		return "qc_failed_r" + release.major + ".tsv";
	}
	throw std::logic_error("unknown GTDB release file");
}

std::string ServedPath(ReleaseFile file, const Release &release) {
	const auto dir = "release" + release.major + "/" + release.ToString() + "/";
	return dir + (file == ReleaseFile::QC_FAILED ? "auxillary_files/" : "") + ReleaseFileName(file, release);
}

std::string CachedPath(ReleaseFile file, const Release &release) {
	return release.ToString() + "/" + ReleaseFileName(file, release);
}

std::string_view FileName(std::string_view path) {
	if (RemoteFileHelper::IsRemotePath(std::string(path))) {
		path = path.substr(0, path.find_first_of("?#"));
	}
	const auto slash = path.find_last_of("/\\");
	return slash == std::string_view::npos ? path : path.substr(slash + 1);
}

std::optional<std::string> MetadataDomain(std::string_view file_name) {
	if (file_name.starts_with("bac120_")) {
		return "Bacteria";
	}
	if (file_name.starts_with("ar53_") || file_name.starts_with("ar122_")) {
		return "Archaea";
	}
	return std::nullopt;
}

std::vector<std::string> MatchMetadataFiles(const std::vector<std::string> &names) {
	std::vector<std::string> bacterial;
	std::vector<std::string> archaeal;
	std::vector<std::string_view> releases;
	for (const auto &name : names) {
		const auto parsed = ParseReleaseFileName(name, {"bac120_metadata", "ar53_metadata", "ar122_metadata"});
		if (!parsed) {
			continue;
		}
		(MetadataDomain(name) == "Bacteria" ? bacterial : archaeal).push_back(name);
		releases.push_back(parsed->release);
	}
	if (bacterial.empty() && archaeal.empty()) {
		throw std::runtime_error("no GTDB metadata file: expected bac120_metadata_r<N>.tsv[.gz] or "
		                         "ar53_metadata_r<N>.tsv[.gz]");
	}
	std::sort(bacterial.begin(), bacterial.end());
	std::sort(archaeal.begin(), archaeal.end());
	RequireAtMostOne(bacterial, "bacterial metadata");
	RequireAtMostOne(archaeal, "archaeal metadata");
	std::vector<std::string> matched = bacterial;
	matched.insert(matched.end(), archaeal.begin(), archaeal.end());
	if (releases.size() == 2 && releases[0] != releases[1]) {
		throw std::runtime_error("metadata files of different releases: " + JoinNames(matched));
	}
	return matched;
}

std::string MatchQcFailedFile(const std::vector<std::string> &names) {
	std::vector<std::string> found;
	for (const auto &name : names) {
		if (ParseReleaseFileName(name, {"qc_failed"})) {
			found.push_back(name);
		}
	}
	if (found.empty()) {
		throw std::runtime_error("no GTDB qc_failed file: expected qc_failed_r<N>.tsv[.gz]");
	}
	std::sort(found.begin(), found.end());
	RequireAtMostOne(found, "qc_failed");
	return found.front();
}

std::vector<std::string> SplitHeader(std::string_view line) {
	std::vector<std::string> names;
	std::unordered_map<std::string, size_t> seen; // lower-cased name -> its column number
	size_t start = 0;
	while (true) {
		const auto tab = line.find('\t', start);
		std::string name(line.substr(start, tab == std::string_view::npos ? std::string_view::npos : tab - start));
		const auto column = names.size() + 1;
		if (name.empty()) {
			throw std::runtime_error("column " + std::to_string(column) + " of the header has no name");
		}
		if (name.find('"') != std::string::npos) {
			throw std::runtime_error("column " + std::to_string(column) + " of the header holds a double quote, '" +
			                         name + "': GTDB's column names never do, and the header is not read as CSV");
		}
		std::string folded = name;
		std::transform(folded.begin(), folded.end(), folded.begin(), duckdb::StringUtil::CharacterToLower);
		const auto inserted = seen.emplace(folded, column);
		if (!inserted.second) {
			const auto first = inserted.first->second;
			throw std::runtime_error("column " + std::to_string(column) + " of the header, '" + name +
			                         "', repeats the name of column " + std::to_string(first) + ", '" +
			                         names[first - 1] + "' (column names ignore case)");
		}
		names.push_back(std::move(name));
		if (tab == std::string_view::npos) {
			return names;
		}
		start = tab + 1;
	}
}

} // namespace gtdb
} // namespace miint
