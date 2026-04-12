#pragma once

#include "duckdb/main/database.hpp"
#include <string>
#include <vector>

#if !defined(_WIN32) && !defined(EMSCRIPTEN)
#define MIINT_ASPERA_SUPPORTED 1
#else
#define MIINT_ASPERA_SUPPORTED 0
#endif

namespace miint {

// Return the temp directory: $TMPDIR if set and non-empty, otherwise /tmp.
std::string GetTempDir();

struct AsperaPath {
	std::string host;
	std::string remote_path;
};

struct AsperaConfig {
	std::string ascp_path;
	std::string key_path;
	std::string host;     // e.g. "fasp.sra.ebi.ac.uk"
	std::string user;     // e.g. "era-fasp"
	int port;             // e.g. 33001
	std::string max_rate; // e.g. "300m"
};

class AsperaUtils {
public:
	// Find ascp binary in PATH. Returns empty string if not found.
	static std::string FindAscp();

	// Find the openssh key at known filesystem locations.
	// If ascp_path is provided, also checks relative to the binary.
	// Returns empty string if not found.
	static std::string FindKey(const std::string &ascp_path = "");

	// Download the well-known Aspera bypass key via HTTP and cache locally.
	// Returns cached path, or empty string on failure.
	static std::string DownloadAndCacheKey(duckdb::DatabaseInstance &db);

	// FindKey() first, then DownloadAndCacheKey() if not found.
	// Throws IOException if both fail and required=true.
	static std::string ResolveKey(duckdb::DatabaseInstance &db, const std::string &ascp_path = "",
	                              bool required = false);

	// Parse ENA fastq_aspera field (semicolon-delimited paths like
	// "fasp.sra.ebi.ac.uk:/vol1/fastq/...;fasp.sra.ebi.ac.uk:/vol1/fastq/...")
	// into {host, remote_path} pairs.
	static std::vector<AsperaPath> ParseAsperaPaths(const std::string &aspera_field);

	// Build a default AsperaConfig from resolved binary and key paths.
	static AsperaConfig BuildConfig(const std::string &ascp_path, const std::string &key_path);

	// Check if Aspera is available (ascp in PATH). Caches result for process lifetime.
	static bool IsAvailable();

	// URL to download the well-known Aspera bypass key
	static constexpr const char *KEY_DOWNLOAD_URL =
	    "https://raw.githubusercontent.com/wwood/kingfisher-download/main/kingfisher/data/asperaweb_id_dsa.openssh";

	// Default cache path for downloaded key
	static constexpr const char *KEY_CACHE_FILENAME = "asperaweb_id_dsa.openssh";
};

} // namespace miint
