#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace miint {

// The taxdump member files consumed by the taxonomy reader. merged.dmp and
// delnodes.dmp are optional (empty string if absent from the archive).
struct TaxdumpFiles {
	std::string nodes;    // nodes.dmp    (required)
	std::string names;    // names.dmp    (required)
	std::string merged;   // merged.dmp   (optional)
	std::string delnodes; // delnodes.dmp (optional)
};

// Decompression + tar extraction for NCBI taxonomy dumps. Pure C++ (no DuckDB):
// the gzip layer is zlib, the tar layer is the vendored microtar. All failures
// throw std::runtime_error.
class TaxdumpArchive {
public:
	// Inflate a gzip (.gz) byte buffer to a string.
	static std::string Gunzip(const std::string &gz_bytes);

	// Extract the named members from an uncompressed tar byte buffer. Members that
	// are not present in the archive are simply absent from the returned map.
	static std::unordered_map<std::string, std::string> ExtractTarMembers(const std::string &tar_bytes,
	                                                                      const std::vector<std::string> &names);

	// Gunzip a taxdump.tar.gz buffer and pull out the taxdump members. Requires
	// nodes.dmp and names.dmp to be present; merged.dmp and delnodes.dmp are
	// optional. Throws if a required member is missing.
	static TaxdumpFiles ExtractTaxdump(const std::string &targz_bytes);
};

} // namespace miint
