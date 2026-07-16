#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace miint {

// One live taxonomy node. Field names mirror the tree-table schema emitted by
// read_newick (taxid -> node_index, parent_taxid -> parent_index) so the taxonomy
// drops straight into the existing tree tooling.
struct TaxdumpNode {
	int64_t taxid;                       // node_index
	std::optional<int64_t> parent_taxid; // parent_index; nullopt for the self-parented root
	std::string name;                    // scientific name (empty if none present)
	std::string rank;                    // raw NCBI rank string (e.g. "superkingdom", "no rank")
	bool is_tip;                         // true iff no node lists this taxid as its parent
};

// A retired -> current taxid remap from merged.dmp.
struct TaxdumpMerge {
	int64_t old_taxid;
	int64_t new_taxid;
};

// Cleanroom parser for NCBI taxonomy dump (taxdump) files. The dump is a set of
// pipe-delimited text files: fields joined by "\t|\t" and each line terminated by
// "\t|". Nothing here touches DuckDB — it is pure string parsing so it can be unit
// tested directly.
class TaxdumpParser {
public:
	// Parse nodes.dmp joined with names.dmp into live tree nodes.
	//   - only name_class == "scientific name" entries populate `name`
	//   - the root (taxid == parent) gets parent_taxid = nullopt
	//   - is_tip is derived from the set of referenced parents
	// Nodes are returned in nodes.dmp order.
	static std::vector<TaxdumpNode> ParseNodes(const std::string &nodes_dmp, const std::string &names_dmp);

	// Parse merged.dmp into (old_taxid, new_taxid) pairs, in file order.
	static std::vector<TaxdumpMerge> ParseMerged(const std::string &merged_dmp);

	// Parse delnodes.dmp into deleted taxids, in file order.
	static std::vector<int64_t> ParseDeleted(const std::string &delnodes_dmp);

	// Split a single dmp line into its field values. Strips a trailing newline and
	// the "\t|" line terminator, then splits on the "\t|\t" field separator.
	// Trailing empty fields are preserved.
	static std::vector<std::string> SplitFields(std::string_view line);
};

} // namespace miint
