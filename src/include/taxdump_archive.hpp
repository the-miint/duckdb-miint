#pragma once

#include <string>
#include <unordered_map>
#include <vector>

namespace miint {

// The taxdump member files consumed by the taxonomy reader. A member that was not
// requested (see TaxdumpMemberSet) is left empty.
struct TaxdumpFiles {
	std::string nodes;    // nodes.dmp
	std::string names;    // names.dmp
	std::string merged;   // merged.dmp
	std::string delnodes; // delnodes.dmp
};

// Which taxdump members a caller needs.
//
// Each reader requests only the members it parses, for two reasons. Cost: the full dump
// is ~510 MB extracted, so loading all four to answer from delnodes.dmp (~7 MB) read
// two orders of magnitude more than the query needed, once per function in a query that
// joined several. Correctness: a member is then required exactly when some caller
// actually reads it, so "absent" can be reported as an error instead of being flattened
// into "empty" -- an empty merge map read as authoritative would assert that every
// retired taxid is live again.
struct TaxdumpMemberSet {
	bool nodes = false;
	bool names = false;
	bool merged = false;
	bool delnodes = false;

	// Every member -- used when populating the on-disk cache, which must be complete
	// even though no single reader needs all four.
	static TaxdumpMemberSet All() {
		return TaxdumpMemberSet {true, true, true, true};
	}
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

	// Gunzip a taxdump.tar.gz buffer and pull out the requested taxdump members.
	// Every requested member must be present: throws std::runtime_error naming the
	// missing one otherwise. Members that were not requested are left empty.
	static TaxdumpFiles ExtractTaxdump(const std::string &targz_bytes, const TaxdumpMemberSet &members);
};

} // namespace miint
