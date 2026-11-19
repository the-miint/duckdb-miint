#pragma once

#include "NewickTree.hpp"
#include "TreeIndex.hpp"
#include "read_newick.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace duckdb {

// Validate that a table or view has the required columns for tree data
// Required: node_index (integer), parent_index (integer, nullable)
// Optional: name (VARCHAR), branch_length (DOUBLE/FLOAT), edge_id (integer)
void ValidateTreeTableSchema(ClientContext &context, const std::string &table_name);

// Read tree node data from a table or view and return as NodeInput structs.
// Uses a separate Connection to avoid deadlocking the current context.
std::vector<miint::NodeInput> ReadTreeTable(ClientContext &context, const std::string &table_name);

// Shared by the functions that return a tree (shear_tree, tree_resolve_multifurcations,
// tree_resolve_placement). A tree's `name` keeps its id type (VARCHAR, BIGINT or UUID) from input
// to output, and text labels live in an optional `label` VARCHAR column, so a tree can hold typed
// tip ids and text labels on its internal nodes together.

// A tree table's `name` and `label` columns, as bound.
struct TreeNameColumns {
	bool has_name = false;
	// VARCHAR when there is no name column
	LogicalType name_type = LogicalType::VARCHAR;
	bool has_label = false;
};

// Bind a tree table for a function that returns a tree: check it (ValidateTreeTableSchema; `name`
// an id type; `label` VARCHAR), throwing BinderException naming `fn` otherwise, and set the output
// schema: read_newick's (without filepath), with `name` in the tree's type, and `label VARCHAR`
// last when the tree has one.
TreeNameColumns BindTreeOutput(ClientContext &context, const char *fn, const std::string &table_name,
                               vector<LogicalType> &return_types, vector<Identifier> &names);

// Read and build a tree bound by BindTreeOutput. Names are read in the bound type and held as text,
// as ReadTreeTable holds them. When the tree has a `label` column, `labels` receives each node's
// label (nullopt for NULL), indexed by the built tree's node index. A column whose type has changed
// since bind is refused (InvalidInputException).
miint::NewickTree BuildTree(ClientContext &context, const std::string &table_name, const TreeNameColumns &columns,
                            std::vector<std::optional<std::string>> &labels);

// The bind data the functions that return a tree share; each function's Data derives from it.
struct TreeOutputData : public TableFunctionData {
	TreeOutputData(std::string tree_table_name_p, TreeNameColumns tree_columns_p)
	    : tree_table_name(std::move(tree_table_name_p)), tree_columns(std::move(tree_columns_p)) {
	}
	std::string tree_table_name;
	TreeNameColumns tree_columns;
};

// The rows a function that returns a tree emits, with their labels.
struct TreeOutputState : public GlobalTableFunctionState {
	std::vector<ReadNewickTableFunction::NodeRow> rows;
	// Each row's label; empty when the tree has no label column
	std::vector<std::optional<std::string>> labels;
	size_t current_row_idx = 0;

	idx_t MaxThreads() const override {
		return 1;
	}

	// Take the output tree's rows, with `labels` (BuildTree's, for the input tree) carried onto them.
	// With `source` (NewickTree::shear's), output node j takes the label of input node source[j];
	// without it, output node i takes input node i's, which suits an operation that keeps every
	// input node at its index and appends the nodes it creates, and those have no label.
	void SetRows(const miint::NewickTree &tree, const std::vector<uint32_t> *source = nullptr);
};

// The table function's execute for the functions that return a tree: emits TreeOutputState's rows,
// with names in TreeOutputData's type.
void TreeOutputExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

// A built tree plus the caller's original node_index for each dense tree index. Because
// ReadTreeTable sorts by node_index and NewickTree::build preserves input order, dense
// index i corresponds to node_ids[i] — a stable join key for reporting results against the
// caller's node_index rather than the internal dense index.
struct TreeWithNodeIds {
	miint::NewickTree tree;
	std::vector<int64_t> node_ids;
};

// Read a tree table/view, capture its node_index mapping, and build the NewickTree. Shared
// by the comparative-methods table functions (independent contrasts, ancestral states,
// parsimony, ML). Throws InvalidInputException if the tree fails to build.
TreeWithNodeIds BuildTreeAndNodeIds(ClientContext &context, const std::string &table_name);

// A tree table read as compact arrays, for functions that need the tree's shape and branch
// lengths but only a few of its names.
struct TreeTableArrays {
	miint::DenseTree tree;
	// (dense index, name) for every node whose name is in the requested set, tips and
	// internal nodes alike, in no particular order. A name on several nodes appears once per node.
	std::vector<std::pair<uint32_t, std::string>> named;
};

// Read a tree table or view once, streaming, keeping each node's name only when it is in
// `keep_names`. ReadTreeTable and BuildTreeAndNodeIds hold every row, name included, twice
// over: loading the Greengenes2 2024.09 tree (24.5M nodes) from Parquet that way was
// killed at an 18 GB memory cap, and this way peaked at 2.6 GB resident. Names are matched
// as text: VARCHAR names as they are, others cast with name::VARCHAR, as ReadTreeTable reads
// them, and not read at all when `keep_names` is empty; a missing name or branch_length
// column reads as NULL.
//
// `name_type` is the type Bind saw for the name column; it is unused when no names are
// read. A view can change type without a catalog change (one over a file that is
// rewritten, under a prepared statement), and names matched as text would then match the
// wrong rows or none, so a name column read with another type is refused. The check uses
// the type of the read itself.
//
// Throws std::invalid_argument when the name column's type has changed, or the rows are not
// a numbered tree (no rows, a NULL or repeated node_index, a parent_index that is no row's
// node_index), so the caller can say which function failed; a failed read throws
// InvalidInputException. The tree's shape (one root, no cycles) and branch lengths are left
// to TreeIndex::Build.
TreeTableArrays ReadTreeArrays(ClientContext &context, const std::string &table_name,
                               const std::unordered_set<std::string_view> &keep_names, const LogicalType &name_type);

} // namespace duckdb
