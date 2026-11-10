#include "GtdbFiles.hpp"

#include <zlib.h>

#include <algorithm>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <unordered_map>

namespace miint {
namespace gtdb {

namespace {

bool IsDigits(std::string_view s) {
	return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) { return c >= '0' && c <= '9'; });
}

bool StartsWith(std::string_view s, std::string_view prefix) {
	return s.substr(0, prefix.size()) == prefix;
}

bool EndsWith(std::string_view s, std::string_view suffix) {
	return s.size() >= suffix.size() && s.substr(s.size() - suffix.size()) == suffix;
}

std::string JoinPath(std::string_view dir, std::string_view name) {
	while (!dir.empty() && dir.back() == '/') {
		dir.remove_suffix(1);
	}
	return std::string(dir) + "/" + std::string(name);
}

// GTDB's name for a file of a release: <prefix>_r<major>.tsv[.gz].
struct NamedFile {
	std::string_view prefix;
	std::string_view release; // the major number in the name
};

// Reads a name of the form <prefix>_r<digits>.tsv or .tsv.gz, for one of `prefixes`.
std::optional<NamedFile> ParseReleaseFileName(std::string_view name, std::initializer_list<std::string_view> prefixes) {
	for (const auto prefix : prefixes) {
		if (!StartsWith(name, prefix) || name.substr(prefix.size(), 2) != "_r") {
			continue;
		}
		auto rest = name.substr(prefix.size() + 2);
		for (const std::string_view extension : {".tsv", ".tsv.gz"}) {
			if (EndsWith(rest, extension) && IsDigits(rest.substr(0, rest.size() - extension.size()))) {
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
	if (!IsDigits(major) || !IsDigits(minor)) {
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

std::string ReleaseFileUrl(std::string_view base_url, ReleaseFile file, const Release &release) {
	auto dir = JoinPath(JoinPath(base_url, "release" + release.major), release.ToString());
	if (file == ReleaseFile::QC_FAILED) {
		dir = JoinPath(dir, "auxillary_files");
	}
	return JoinPath(dir, ReleaseFileName(file, release));
}

std::string ReleaseCachePath(std::string_view cache_dir, ReleaseFile file, const Release &release) {
	return JoinPath(JoinPath(cache_dir, release.ToString()), ReleaseFileName(file, release));
}

std::string_view FileName(std::string_view path) {
	if (StartsWith(path, "http://") || StartsWith(path, "https://")) {
		path = path.substr(0, path.find_first_of("?#"));
	}
	const auto slash = path.find_last_of("/\\");
	return slash == std::string_view::npos ? path : path.substr(slash + 1);
}

std::optional<std::string> MetadataDomain(std::string_view file_name) {
	if (StartsWith(file_name, "bac120_")) {
		return "Bacteria";
	}
	if (StartsWith(file_name, "ar53_") || StartsWith(file_name, "ar122_")) {
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
		(parsed->prefix == "bac120_metadata" ? bacterial : archaeal).push_back(name);
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
		std::transform(folded.begin(), folded.end(), folded.begin(),
		               [](char c) { return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c; });
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

// ---------------------------------------------------------------------------
// HeaderLineReader
// ---------------------------------------------------------------------------

struct HeaderLineReader::Inflater {
	z_stream zs;
	bool ended = false;

	Inflater() {
		std::memset(&zs, 0, sizeof(zs));
		// 16 + MAX_WBITS: gzip framing.
		if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
			throw std::runtime_error("cannot start a gzip decompressor");
		}
	}
	~Inflater() {
		inflateEnd(&zs);
	}
	Inflater(const Inflater &) = delete;
	Inflater &operator=(const Inflater &) = delete;
};

namespace {
constexpr size_t GZIP_MAGIC_SIZE = 2;
constexpr char GZIP_MAGIC[GZIP_MAGIC_SIZE] = {'\x1f', '\x8b'};
constexpr size_t INFLATE_WINDOW = 65536;

std::runtime_error NoFirstLine(const char *why) {
	return std::runtime_error(std::string("no first line: ") + why);
}
} // namespace

HeaderLineReader::HeaderLineReader(size_t max_line) : max_line_(max_line) {
}
HeaderLineReader::~HeaderLineReader() = default;
HeaderLineReader::HeaderLineReader(HeaderLineReader &&) noexcept = default;
HeaderLineReader &HeaderLineReader::operator=(HeaderLineReader &&) noexcept = default;

bool HeaderLineReader::Feed(const char *data, size_t size) {
	if (done_) {
		return true;
	}
	if (!decided_) {
		const auto take = std::min(size, GZIP_MAGIC_SIZE - pending_.size());
		pending_.append(data, take);
		data += take;
		size -= take;
		if (pending_.size() < GZIP_MAGIC_SIZE) {
			return false;
		}
		if (Decide()) {
			return true;
		}
	}
	return gzip_ ? Inflate(data, size) : Append(data, size);
}

// Settles plain or gzip from the bytes held so far, and hands them on.
bool HeaderLineReader::Decide() {
	decided_ = true;
	gzip_ = pending_.size() == GZIP_MAGIC_SIZE && std::memcmp(pending_.data(), GZIP_MAGIC, GZIP_MAGIC_SIZE) == 0;
	if (gzip_) {
		inflater_ = std::make_unique<Inflater>();
	}
	std::string held;
	held.swap(pending_);
	return gzip_ ? Inflate(held.data(), held.size()) : Append(held.data(), held.size());
}

// Adds decoded bytes to the line. Returns true once it holds the whole line.
bool HeaderLineReader::Append(const char *data, size_t size) {
	const auto newline = static_cast<const char *>(std::memchr(data, '\n', size));
	line_.append(data, newline ? static_cast<size_t>(newline - data) : size);
	if (newline) {
		if (!line_.empty() && line_.back() == '\r') {
			line_.pop_back();
		}
		done_ = true;
	}
	// A '\r' that may yet turn out to end the line does not count against the cap.
	const auto length = !done_ && !line_.empty() && line_.back() == '\r' ? line_.size() - 1 : line_.size();
	if (length > max_line_) {
		throw std::runtime_error("the first line is longer than " + std::to_string(max_line_) + " bytes");
	}
	return done_;
}

// A gzip file may hold several members, one after another (bgzip writes them): when one
// ends, the bytes that follow start the next.
bool HeaderLineReader::Inflate(const char *data, size_t size) {
	auto &zs = inflater_->zs;
	char window[INFLATE_WINDOW];
	while (size > 0) {
		if (inflater_->ended) {
			if (inflateReset(&zs) != Z_OK) {
				throw std::runtime_error("cannot restart the gzip decompressor");
			}
			inflater_->ended = false;
		}
		// avail_in is a uInt: hand a large buffer over in slices.
		const auto slice = static_cast<uInt>(std::min<size_t>(size, 1u << 30));
		zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(data));
		zs.avail_in = slice;
		do {
			zs.next_out = reinterpret_cast<Bytef *>(window);
			zs.avail_out = sizeof(window);
			const int ret = inflate(&zs, Z_NO_FLUSH);
			if (ret == Z_STREAM_END) {
				inflater_->ended = true;
			} else if (ret != Z_OK && ret != Z_BUF_ERROR) {
				throw std::runtime_error(std::string("corrupt gzip stream: ") + (zs.msg ? zs.msg : "inflate failed"));
			}
			if (Append(window, sizeof(window) - zs.avail_out)) {
				return true;
			}
		} while (!inflater_->ended && (zs.avail_in > 0 || zs.avail_out == 0));
		data += slice - zs.avail_in;
		size -= slice - zs.avail_in;
	}
	return false;
}

std::string HeaderLineReader::Finish() {
	if (!decided_) {
		Decide();
	}
	if (!done_) {
		if (!gzip_) {
			throw NoFirstLine("the file ends before its first newline");
		}
		throw NoFirstLine(inflater_ && inflater_->ended ? "the gzip stream ends before its first newline"
		                                                : "the gzip stream is cut short before its first newline");
	}
	return line_;
}

} // namespace gtdb
} // namespace miint
