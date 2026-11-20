#include "tree_table_reader.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "read_newick.hpp"
#include "tree_function_utils.hpp"
#include "miint_streaming_query.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace duckdb {

static void CheckTreeColumns(const TableOrViewColumns &info, const std::string &table_name) {
	// node_index and parent_index are required; name/branch_length/edge_id are
	// optional (ReadTreeTable substitutes NULL when they are absent).
	for (const char *required : {"node_index", "parent_index"}) {
		if (!HasColumn(info, required)) {
			throw BinderException("Tree table '%s' missing required column '%s'", table_name, required);
		}
	}
}

void ValidateTreeTableSchema(ClientContext &context, const std::string &table_name) {
	CheckTreeColumns(GetTableOrViewColumns(context, table_name, "Tree table"), table_name);
}

// The rows of a tree table, in node_index order. Without `columns`, names are cast to VARCHAR, as
// ReadTreeTable reads them; with them (BindTreeOutput's, for the function `fn`), names are read in
// the bound type and, with `labels`, the VARCHAR `label` column is read too, each checked against
// bind.
static std::vector<miint::NodeInput> ReadTreeRows(ClientContext &context, const std::string &table_name, const char *fn,
                                                  const TreeNameColumns *columns,
                                                  std::vector<std::optional<std::string>> *labels) {
	std::vector<miint::NodeInput> result;

	auto conn = MakeReadOnlyHelperConnection(context);

	// Optional columns (name, branch_length, edge_id) default to NULL when the
	// tree table lacks them, so a minimal node_index+parent_index table is
	// accepted rather than failing mid-read. node_index/parent_index are
	// guaranteed present by the bind-time check (ValidateTreeTableSchema or
	// BindTreeOutput).
	//
	// The numeric columns are cast to canonical types so GetData<> always
	// matches the backing store (DuckDB handles INTEGER->BIGINT, FLOAT->DOUBLE
	// etc. transparently). Names are cast to VARCHAR, or read in their bound
	// type, which is checked below.
	auto info = GetTableOrViewColumns(context, table_name, "Tree table");
	const bool has_name = columns ? columns->has_name : HasColumn(info, "name");
	const LogicalType name_type = columns ? columns->name_type : LogicalType::VARCHAR;
	std::string name_proj = "CAST(NULL AS VARCHAR)";
	if (has_name) {
		name_proj = columns ? "name" : "name::VARCHAR";
	}
	std::string bl_proj = HasColumn(info, "branch_length") ? "branch_length::DOUBLE" : "CAST(NULL AS DOUBLE)";
	std::string edge_proj = HasColumn(info, "edge_id") ? "edge_id::BIGINT" : "CAST(NULL AS BIGINT)";

	// ORDER BY node_index: build() assigns node-array positions (and therefore
	// each node's slot in its parent's children list) by the row order this query
	// returns. Without an explicit order, that follows scan order, which is not
	// guaranteed to match node_index for views/joins, Parquet-backed tables,
	// parallel scans, or preserve_insertion_order=false. node_index is the
	// canonical child order (it encodes read_newick's original left-to-right
	// order), so ordering by it makes order-sensitive consumers
	// (tree_resolve_multifurcations) reproducible regardless of storage.
	std::string query = "SELECT node_index::BIGINT, parent_index::BIGINT, " + name_proj + ", " + bl_proj + ", " +
	                    edge_proj + (labels ? ", label" : "") + " FROM " + SQLIdentifier::ToString(table_name) +
	                    " ORDER BY node_index";

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from tree table '%s': %s", table_name, query_result->GetError());
	}

	auto &materialized = *query_result;

	// Bind fixed the output types; a name read uncast, or a label, must still have them.
	const std::string relation = "tree table '" + table_name + "'";
	if (columns && columns->has_name) {
		CheckBoundType(fn, relation, "name", materialized.GetTypes()[2], name_type);
	}
	if (labels) {
		CheckBoundType(fn, relation, "label", materialized.GetTypes()[5], LogicalType::VARCHAR);
	}

	result.reserve(materialized.RowCount());
	if (labels) {
		labels->reserve(materialized.RowCount());
	}

	std::vector<std::string> names;
	std::vector<bool> name_nulls;
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &node_index_vec = chunk->data[0];
		auto &parent_index_vec = chunk->data[1];
		auto &branch_length_vec = chunk->data[3];
		auto &edge_id_vec = chunk->data[4];

		UnifiedVectorFormat node_data, parent_data, bl_data, edge_data, label_data;
		node_index_vec.ToUnifiedFormat(node_data);
		parent_index_vec.ToUnifiedFormat(parent_data);
		branch_length_vec.ToUnifiedFormat(bl_data);
		edge_id_vec.ToUnifiedFormat(edge_data);
		ExtractIdColumnAsStrings(*chunk, 2, name_type, names, name_nulls);
		if (labels) {
			chunk->data[5].ToUnifiedFormat(label_data);
		}

		auto node_indices = UnifiedVectorFormat::GetData<int64_t>(node_data);
		auto parent_indices = UnifiedVectorFormat::GetData<int64_t>(parent_data);
		auto branch_lengths = UnifiedVectorFormat::GetData<double>(bl_data);
		auto edge_ids = UnifiedVectorFormat::GetData<int64_t>(edge_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			miint::NodeInput node;

			auto ni_idx = node_data.sel->get_index(i);
			if (!node_data.validity.RowIsValid(ni_idx)) {
				throw InvalidInputException("NULL node_index in tree table '%s'", table_name);
			}
			node.node_id = node_indices[ni_idx];

			auto pi_idx = parent_data.sel->get_index(i);
			if (parent_data.validity.RowIsValid(pi_idx)) {
				node.parent_id = parent_indices[pi_idx];
			} else {
				node.parent_id = std::nullopt; // Root
			}

			// A NULL name is held as '', as an unnamed node is
			node.name = std::move(names[i]);

			node.branch_length = std::numeric_limits<double>::quiet_NaN();
			auto bl_idx = bl_data.sel->get_index(i);
			if (bl_data.validity.RowIsValid(bl_idx)) {
				node.branch_length = branch_lengths[bl_idx];
			}

			auto ei_idx = edge_data.sel->get_index(i);
			if (edge_data.validity.RowIsValid(ei_idx)) {
				node.edge_id = edge_ids[ei_idx];
			}

			if (labels) {
				auto lb_idx = label_data.sel->get_index(i);
				if (label_data.validity.RowIsValid(lb_idx)) {
					labels->emplace_back(UnifiedVectorFormat::GetData<string_t>(label_data)[lb_idx].GetString());
				} else {
					labels->emplace_back(std::nullopt);
				}
			}

			result.push_back(std::move(node));
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Tree table '%s' is empty", table_name);
	}

	return result;
}

std::vector<miint::NodeInput> ReadTreeTable(ClientContext &context, const std::string &table_name) {
	return ReadTreeRows(context, table_name, nullptr, nullptr, nullptr);
}

// Build a tree from rows read from `table_name`, naming it in the error.
static miint::NewickTree BuildNewickTree(const std::vector<miint::NodeInput> &rows, const std::string &table_name) {
	try {
		return miint::NewickTree::build(rows);
	} catch (const std::exception &e) {
		throw InvalidInputException("Failed to build tree from '%s': %s", table_name, e.what());
	}
}

TreeNameColumns BindTreeOutput(ClientContext &context, const char *fn, const std::string &table_name,
                               vector<LogicalType> &return_types, vector<Identifier> &names) {
	auto info = GetTableOrViewColumns(context, table_name, "Tree table");
	CheckTreeColumns(info, table_name);
	TreeNameColumns columns;
	columns.has_name = HasColumn(info, "name");
	if (columns.has_name) {
		columns.name_type = TipNameType(fn, table_name, info);
	}
	LogicalType label_type;
	if (FindColumnType(info, "label", label_type)) {
		if (label_type.id() != LogicalTypeId::VARCHAR) {
			throw BinderException("%s: tree table '%s' column 'label' must be VARCHAR, got %s", fn, table_name,
			                      label_type.ToString());
		}
		columns.has_label = true;
	}

	std::vector<std::string> schema_names;
	std::vector<LogicalType> schema_types;
	ReadNewickTableFunction::GetSchema(schema_names, schema_types, false);
	schema_types[1] = columns.name_type;
	if (columns.has_label) {
		schema_names.emplace_back("label");
		schema_types.emplace_back(LogicalType::VARCHAR);
	}
	for (const auto &name : schema_names) {
		names.emplace_back(name);
	}
	return_types.insert(return_types.end(), schema_types.begin(), schema_types.end());
	return columns;
}

miint::NewickTree BuildTree(ClientContext &context, const char *fn, const std::string &table_name,
                            const TreeNameColumns &columns, std::vector<std::optional<std::string>> &labels) {
	labels.clear();
	// The rows are released on return, before the caller's operation, so they don't add to its peak
	// memory: an operation that builds a new tree (shear, resolution) holds two trees at once.
	auto node_inputs = ReadTreeRows(context, table_name, fn, &columns, columns.has_label ? &labels : nullptr);
	// build() keeps the rows' order, so node i is row i, the row labels[i] came from
	return BuildNewickTree(node_inputs, table_name);
}

void TreeOutputState::SetRows(const miint::NewickTree &tree, const std::vector<uint32_t> *source) {
	rows = ReadNewickTableFunction::TreeToRows(tree);
	if (labels.empty()) {
		return;
	}
	if (source) {
		std::vector<std::optional<std::string>> carried;
		carried.reserve(source->size());
		for (auto orig : *source) {
			carried.push_back(std::move(labels[orig]));
		}
		labels = std::move(carried);
	} else {
		labels.resize(tree.num_nodes());
	}
}

void TreeOutputExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<TreeOutputData>();
	auto &state = data_p.global_state->Cast<TreeOutputState>();

	if (state.current_row_idx >= state.rows.size()) {
		output.SetChildCardinality(0);
		return;
	}

	size_t rows_to_output = std::min<size_t>(STANDARD_VECTOR_SIZE, state.rows.size() - state.current_row_idx);

	ReadNewickTableFunction::EmitNodeRows(state.rows, state.current_row_idx, rows_to_output, output, 0, false, "",
	                                      bind_data.tree_columns.name_type,
	                                      bind_data.tree_columns.has_label ? &state.labels : nullptr);

	state.current_row_idx += rows_to_output;
	output.SetChildCardinality(rows_to_output);
}

TreeWithNodeIds BuildTreeAndNodeIds(ClientContext &context, const std::string &table_name) {
	auto node_inputs = ReadTreeTable(context, table_name);

	// build() preserves input order for node indices, and ReadTreeTable sorts by
	// node_index, so tree dense index i corresponds to node_inputs[i]. Capture the
	// original node_index up front so results are reported against the caller's
	// node_index (a stable join key), not the internal dense index.
	TreeWithNodeIds out;
	out.node_ids.reserve(node_inputs.size());
	for (const auto &ni : node_inputs) {
		out.node_ids.push_back(ni.node_id);
	}

	out.tree = BuildNewickTree(node_inputs, table_name);
	return out;
}

TreeTableArrays ReadTreeArrays(ClientContext &context, const std::string &table_name,
                               const std::unordered_set<std::string_view> &keep_names, const LogicalType &name_type) {
	auto info = GetTableOrViewColumns(context, table_name, "Tree table");
	// Names are read only when some are kept.
	const bool has_names = HasColumn(info, "name") && !keep_names.empty();
	// The names are matched as text. VARCHAR names, the usual case, are read as they are, so
	// the type of the read shows whether they still are VARCHAR; it is checked before any
	// row is read as text. Names bound with another type are cast to text, and the column as
	// it is read is projected too, last, only so that its type can be checked: a second copy
	// of every VARCHAR name would cost more.
	const bool names_uncast = name_type.id() == LogicalTypeId::VARCHAR;
	const idx_t name_type_col = names_uncast ? 2 : 4;
	std::string name_proj = "CAST(NULL AS VARCHAR)";
	if (has_names) {
		name_proj = names_uncast ? "name" : "name::VARCHAR";
	}
	const std::string bl_proj = HasColumn(info, "branch_length") ? "branch_length::DOUBLE" : "CAST(NULL AS DOUBLE)";
	// No ORDER BY: DenseTree numbers the rows itself, and sorting here would carry every
	// name through the sort.
	const std::string query = "SELECT node_index::BIGINT, parent_index::BIGINT, " + name_proj + ", " + bl_proj +
	                          (has_names && !names_uncast ? ", name" : "") + " FROM " +
	                          SQLIdentifier::ToString(table_name);

	auto conn = MakeReadOnlyHelperConnection(context);
	auto stream = SubmitStream(conn, query);
	if (stream->HasError()) {
		throw InvalidInputException("Failed to read from tree table '%s': %s", table_name, stream->GetError());
	}
	if (has_names && stream->GetTypes()[name_type_col] != name_type) {
		throw std::invalid_argument("column 'name' is " + stream->GetTypes()[name_type_col].ToString() +
		                            " now but was " + name_type.ToString() +
		                            " when the statement was bound; run or prepare it again");
	}

	std::vector<int64_t> node_ids;
	std::vector<int64_t> parent_ids;
	std::vector<bool> has_parent;
	std::vector<double> branch_lengths;
	std::vector<std::pair<int64_t, std::string>> named_ids;
	while (auto chunk = stream->Fetch()) {
		const idx_t rows = chunk->size();
		if (rows == 0) {
			break;
		}
		context.InterruptCheck();
		UnifiedVectorFormat node_data, parent_data, name_data, bl_data;
		chunk->data[0].ToUnifiedFormat(node_data);
		chunk->data[1].ToUnifiedFormat(parent_data);
		chunk->data[2].ToUnifiedFormat(name_data);
		chunk->data[3].ToUnifiedFormat(bl_data);
		auto node_values = UnifiedVectorFormat::GetData<int64_t>(node_data);
		auto parent_values = UnifiedVectorFormat::GetData<int64_t>(parent_data);
		auto name_values = UnifiedVectorFormat::GetData<string_t>(name_data);
		auto bl_values = UnifiedVectorFormat::GetData<double>(bl_data);
		for (idx_t i = 0; i < rows; ++i) {
			const auto ni = node_data.sel->get_index(i);
			if (!node_data.validity.RowIsValid(ni)) {
				throw std::invalid_argument("a row's node_index is NULL");
			}
			node_ids.push_back(node_values[ni]);

			const auto pi = parent_data.sel->get_index(i);
			const bool parent_valid = parent_data.validity.RowIsValid(pi);
			has_parent.push_back(parent_valid);
			parent_ids.push_back(parent_valid ? parent_values[pi] : 0);

			const auto bi = bl_data.sel->get_index(i);
			branch_lengths.push_back(bl_data.validity.RowIsValid(bi) ? bl_values[bi]
			                                                         : std::numeric_limits<double>::quiet_NaN());

			const auto mi = name_data.sel->get_index(i);
			if (name_data.validity.RowIsValid(mi)) {
				const std::string_view name(name_values[mi].GetData(), name_values[mi].GetSize());
				if (keep_names.count(name)) {
					named_ids.emplace_back(node_values[ni], std::string(name));
				}
			}
		}
	}
	if (stream->HasError()) {
		throw InvalidInputException("Failed to read from tree table '%s': %s", table_name, stream->GetError());
	}

	TreeTableArrays out;
	out.tree = miint::DenseTree::FromRows(std::move(node_ids), std::move(parent_ids), std::move(has_parent),
	                                      std::move(branch_lengths));
	out.named.reserve(named_ids.size());
	for (auto &named : named_ids) {
		out.named.emplace_back(out.tree.find(named.first), std::move(named.second));
	}
	return out;
}

} // namespace duckdb
