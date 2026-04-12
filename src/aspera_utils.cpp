#include "aspera_utils.hpp"

namespace miint {

std::string GetTempDir() {
	const char *tmpdir = getenv("TMPDIR");
	if (tmpdir && tmpdir[0] != '\0') {
		return std::string(tmpdir);
	}
	return "/tmp";
}

} // namespace miint

#if MIINT_ASPERA_SUPPORTED

#include "duckdb/common/exception.hpp"
#include "duckdb/common/http_util.hpp"
#include <array>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <mutex>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

namespace miint {

// ---------------------------------------------------------------------------
// Binary discovery (fork/exec pattern from Bowtie2Aligner)
// ---------------------------------------------------------------------------

std::string AsperaUtils::FindAscp() {
	int pipefd[2];
	if (pipe(pipefd) == -1) {
		return "";
	}

	pid_t pid = fork();
	if (pid == -1) {
		close(pipefd[0]);
		close(pipefd[1]);
		return "";
	}

	if (pid == 0) {
		// Child
		close(pipefd[0]);
		dup2(pipefd[1], STDOUT_FILENO);
		close(pipefd[1]);
		int devnull = open("/dev/null", O_WRONLY);
		if (devnull != -1) {
			dup2(devnull, STDERR_FILENO);
			close(devnull);
		}
		execlp("which", "which", "ascp", nullptr);
		_exit(127);
	}

	// Parent
	close(pipefd[1]);
	std::array<char, 512> buffer;
	std::string result;
	ssize_t n;
	while ((n = read(pipefd[0], buffer.data(), buffer.size())) > 0) {
		result.append(buffer.data(), static_cast<size_t>(n));
	}
	close(pipefd[0]);

	int status;
	waitpid(pid, &status, 0);
	if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		return "";
	}

	while (!result.empty() && (result.back() == '\n' || result.back() == '\r')) {
		result.pop_back();
	}
	return result;
}

// ---------------------------------------------------------------------------
// Key discovery
// ---------------------------------------------------------------------------

static bool FileExists(const std::string &path) {
	struct stat st;
	return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

static std::string HomeDir() {
	const char *home = getenv("HOME");
	return home ? std::string(home) : "";
}

static std::string ParentDir(const std::string &path) {
	auto pos = path.rfind('/');
	if (pos == std::string::npos || pos == 0) {
		return "";
	}
	auto parent = path.substr(0, pos);
	pos = parent.rfind('/');
	if (pos == std::string::npos) {
		return "";
	}
	return parent.substr(0, pos);
}

std::string AsperaUtils::FindKey(const std::string &ascp_path) {
	std::vector<std::string> candidates;

	// 1. CONDA_PREFIX
	const char *conda = getenv("CONDA_PREFIX");
	if (conda) {
		candidates.push_back(std::string(conda) + "/etc/asperaweb_id_dsa.openssh");
	}

	// 2-4. Well-known home paths
	std::string home = HomeDir();
	if (!home.empty()) {
		candidates.push_back(home + "/.aspera/connect/etc/asperaweb_id_dsa.openssh");
		candidates.push_back(home + "/.aspera/cli/etc/asperaweb_id_dsa.openssh");
		candidates.push_back(home + "/.aspera/sdk/aspera_bypass_dsa.pem");
	}

	// 5. Relative to ascp binary: $(dirname $(dirname ascp))/etc/asperaweb_id_dsa.openssh
	if (!ascp_path.empty()) {
		std::string parent = ParentDir(ascp_path);
		if (!parent.empty()) {
			candidates.push_back(parent + "/etc/asperaweb_id_dsa.openssh");
			candidates.push_back(parent + "/etc/aspera_bypass_dsa.pem");
		}
	}

	for (const auto &path : candidates) {
		if (FileExists(path)) {
			return path;
		}
	}
	return "";
}

// ---------------------------------------------------------------------------
// Key download + cache (requires DuckDB HTTP API — only in static build)
// ---------------------------------------------------------------------------

#ifdef MIINT_STATIC_BUILD

static bool MkdirP(const std::string &path) {
	std::string current;
	for (size_t i = 0; i < path.size(); i++) {
		current += path[i];
		if (path[i] == '/' && i > 0) {
			mkdir(current.c_str(), 0755);
		}
	}
	return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
}

std::string AsperaUtils::DownloadAndCacheKey(duckdb::DatabaseInstance &db) {
	std::string home = HomeDir();
	if (home.empty()) {
		return "";
	}

	std::string cache_dir = home + "/.aspera/connect/etc";
	std::string cache_path = cache_dir + "/" + KEY_CACHE_FILENAME;

	if (FileExists(cache_path)) {
		return cache_path;
	}

	try {
		std::string url(KEY_DOWNLOAD_URL);
		duckdb::HTTPHeaders headers(db);
		auto &http_util = duckdb::HTTPUtil::Get(db);
		auto params = http_util.InitializeParameters(db, url);
		duckdb::GetRequestInfo get_request(url, headers, *params, nullptr, nullptr);
		auto response = http_util.Request(get_request);

		if (!response->Success() || response->body.empty()) {
			return "";
		}

		if (response->body.find("-----BEGIN DSA PRIVATE KEY-----") == std::string::npos ||
		    response->body.find("-----END DSA PRIVATE KEY-----") == std::string::npos) {
			return "";
		}

		MkdirP(cache_dir);
		std::ofstream out(cache_path, std::ios::binary);
		if (!out) {
			return "";
		}
		out.write(response->body.data(), static_cast<std::streamsize>(response->body.size()));
		out.close();

		chmod(cache_path.c_str(), 0600);
		return cache_path;
	} catch (...) {
		return "";
	}
}

std::string AsperaUtils::ResolveKey(duckdb::DatabaseInstance &db, const std::string &ascp_path, bool required) {
	std::string key = FindKey(ascp_path);
	if (!key.empty()) {
		return key;
	}

	key = DownloadAndCacheKey(db);
	if (!key.empty()) {
		return key;
	}

	if (required) {
		throw duckdb::IOException("read_ena_fastx: Aspera SSH key not found. Searched: ~/.aspera/connect/etc/, "
		                          "~/.aspera/cli/etc/, $CONDA_PREFIX/etc/, and attempted download from %s",
		                          KEY_DOWNLOAD_URL);
	}
	return "";
}

#else // !MIINT_STATIC_BUILD

std::string AsperaUtils::DownloadAndCacheKey(duckdb::DatabaseInstance &) {
	return "";
}

std::string AsperaUtils::ResolveKey(duckdb::DatabaseInstance &, const std::string &ascp_path, bool required) {
	std::string key = FindKey(ascp_path);
	if (!key.empty()) {
		return key;
	}
	if (required) {
		throw std::runtime_error("Aspera SSH key not found");
	}
	return "";
}

#endif // MIINT_STATIC_BUILD

AsperaConfig AsperaUtils::BuildConfig(const std::string &ascp_path, const std::string &key_path) {
	return AsperaConfig {ascp_path, key_path, "fasp.sra.ebi.ac.uk", "era-fasp", 33001, "300m"};
}

bool AsperaUtils::IsAvailable() {
	static std::once_flag flag;
	static bool result = false;
	std::call_once(flag, []() { result = !FindAscp().empty(); });
	return result;
}

} // namespace miint

#else // !MIINT_ASPERA_SUPPORTED

namespace miint {

std::string AsperaUtils::FindAscp() {
	return "";
}
std::string AsperaUtils::FindKey(const std::string &) {
	return "";
}
std::string AsperaUtils::DownloadAndCacheKey(duckdb::DatabaseInstance &) {
	return "";
}
std::string AsperaUtils::ResolveKey(duckdb::DatabaseInstance &, const std::string &, bool required) {
	if (required) {
		throw std::runtime_error("read_ena_fastx: Aspera is not supported on this platform");
	}
	return "";
}
AsperaConfig AsperaUtils::BuildConfig(const std::string &ascp_path, const std::string &key_path) {
	return AsperaConfig {ascp_path, key_path, "fasp.sra.ebi.ac.uk", "era-fasp", 33001, "300m"};
}
bool AsperaUtils::IsAvailable() {
	return false;
}

} // namespace miint

#endif // MIINT_ASPERA_SUPPORTED

// Platform-independent functions (outside #if guards)
namespace miint {

std::vector<AsperaPath> AsperaUtils::ParseAsperaPaths(const std::string &aspera_field) {
	std::vector<AsperaPath> paths;
	if (aspera_field.empty()) {
		return paths;
	}

	std::string::size_type start = 0;
	while (true) {
		auto semi = aspera_field.find(';', start);
		std::string entry;
		if (semi == std::string::npos) {
			entry = aspera_field.substr(start);
		} else {
			entry = aspera_field.substr(start, semi - start);
		}

		if (!entry.empty()) {
			auto colon = entry.find(':');
			if (colon != std::string::npos) {
				paths.push_back({entry.substr(0, colon), entry.substr(colon + 1)});
			} else {
				paths.push_back({"fasp.sra.ebi.ac.uk", "/" + entry});
			}
		}

		if (semi == std::string::npos) {
			break;
		}
		start = semi + 1;
	}
	return paths;
}

} // namespace miint
