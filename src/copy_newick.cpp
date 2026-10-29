#include "copy_newick.hpp"
#include "copy_format_common.hpp"
#include "bgzf_write_mode.hpp"
#include "NewickTree.hpp"
#include "placement_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/execution/operator/persistent/physical_copy_to_file.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/planner/expression/bound_cast_expression.hpp"

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct NewickCopyBindData : public FunctionData {
	std::string file_path;
	bool include_edge_ids = false;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	int compression_level = ::miint::DEFAULT_GZIP_LEVEL; // COMPRESSION_LEVEL, for local gzip output

	// The PLACEMENTS relation. Its rows are read when the COPY runs (InitializeGlobal), so
	// planning reads none and a prepared COPY writes the placements as they are when executed.
	std::optional<std::string> placements_table;

	// Column indices (DConstants::INVALID_INDEX if not present)
	idx_t node_index_idx = DConstants::INVALID_INDEX;
	idx_t parent_index_idx = DConstants::INVALID_INDEX;
	idx_t name_idx = DConstants::INVALID_INDEX;
	idx_t branch_length_idx = DConstants::INVALID_INDEX;
	idx_t edge_id_idx = DConstants::INVALID_INDEX;

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<NewickCopyBindData>();
		result->file_path = file_path;
		result->include_edge_ids = include_edge_ids;
		result->compression = compression;
		result->compression_level = compression_level;
		result->placements_table = placements_table;
		result->node_index_idx = node_index_idx;
		result->parent_index_idx = parent_index_idx;
		result->name_idx = name_idx;
		result->branch_length_idx = branch_length_idx;
		result->edge_id_idx = edge_id_idx;
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<NewickCopyBindData>();
		return file_path == other.file_path && include_edge_ids == other.include_edge_ids &&
		       compression == other.compression && compression_level == other.compression_level &&
		       placements_table == other.placements_table && node_index_idx == other.node_index_idx &&
		       parent_index_idx == other.parent_index_idx && name_idx == other.name_idx &&
		       branch_length_idx == other.branch_length_idx && edge_id_idx == other.edge_id_idx;
	}
};

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct NewickCopyGlobalState : public GlobalFunctionData {
	std::mutex lock;
	std::string file_path;
	FileCompressionType compression;
	int compression_level;
	bool include_edge_ids;

	// Read from the PLACEMENTS relation when the COPY started; empty without one
	std::vector<miint::Placement> placements;

	// File handle (opened in finalize, not at init)
	unique_ptr<CopyFileHandle> file;

	// Collected node data from all threads
	std::vector<miint::NodeInput> nodes;

	NewickCopyGlobalState(const std::string &path, FileCompressionType comp, int level, bool edge_ids,
	                      std::vector<miint::Placement> placements_p)
	    : file_path(path), compression(comp), compression_level(level), include_edge_ids(edge_ids),
	      placements(std::move(placements_p)) {
	}
};

//===--------------------------------------------------------------------===//
// Local State
//===--------------------------------------------------------------------===//
struct NewickCopyLocalState : public LocalFunctionData {
	// Local buffer to collect nodes before merging into global state
	std::vector<miint::NodeInput> local_nodes;
};

//===--------------------------------------------------------------------===//
// Select: cast the columns the sink reads to the types it reads them as
//===--------------------------------------------------------------------===//
// The columns the sink reads, where bind records each, the one type the sink reads it as, and the
// types NewickCopySelect casts to it. An untyped NULL column is cast too, so `NULL AS name` writes
// no names.
struct NewickColumnType {
	const char *name;
	idx_t NewickCopyBindData::*index;
	LogicalTypeId type;
	bool (*castable)(const LogicalType &type);
	const char *description;
};

static bool IsIntegralType(const LogicalType &type) {
	return type.IsIntegral();
}

static bool IsNumericType(const LogicalType &type) {
	return type.IsNumeric();
}

static const NewickColumnType NEWICK_COLUMN_TYPES[] = {
    {"node_index", &NewickCopyBindData::node_index_idx, LogicalTypeId::BIGINT, IsIntegralType, "an integer type"},
    {"parent_index", &NewickCopyBindData::parent_index_idx, LogicalTypeId::BIGINT, IsIntegralType, "an integer type"},
    {"edge_id", &NewickCopyBindData::edge_id_idx, LogicalTypeId::BIGINT, IsIntegralType, "an integer type"},
    {"branch_length", &NewickCopyBindData::branch_length_idx, LogicalTypeId::DOUBLE, IsNumericType, "a numeric type"},
    {"name", &NewickCopyBindData::name_idx, LogicalTypeId::VARCHAR, nullptr, "VARCHAR"},
};

static const NewickColumnType *FindNewickColumnType(const Identifier &name) {
	for (auto &column : NEWICK_COLUMN_TYPES) {
		if (name == column.name) {
			return &column;
		}
	}
	return nullptr;
}

static vector<unique_ptr<Expression>> NewickCopySelect(CopyToSelectInput &input) {
	vector<unique_ptr<Expression>> result;
	bool any_cast = false;
	for (auto &expr : input.select_list) {
		const auto &type = expr->GetReturnType();
		const auto name = expr->GetAlias();
		auto column = FindNewickColumnType(name);
		const bool cast = column && type.id() != column->type &&
		                  (type.id() == LogicalTypeId::SQLNULL || (column->castable && column->castable(type)));
		if (!cast) {
			result.push_back(std::move(expr));
			continue;
		}
		auto cast_expr =
		    BoundCastExpression::AddCastToType(input.context, std::move(expr), LogicalType(column->type), false);
		cast_expr->SetAlias(name);
		result.push_back(std::move(cast_expr));
		any_cast = true;
	}
	// An empty list leaves the input as it is
	if (!any_cast) {
		result.clear();
	}
	return result;
}

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> NewickCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                               const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<NewickCopyBindData>();
	result->file_path = input.info.file_path;

	// Find the columns the sink reads. One still not of the sink's type is one NewickCopySelect could
	// not cast; other columns (like is_tip, filepath) are ignored.
	for (idx_t i = 0; i < names.size(); i++) {
		auto column = FindNewickColumnType(names[i]);
		if (!column) {
			continue;
		}
		if (sql_types[i].id() != column->type) {
			throw BinderException("COPY FORMAT NEWICK: column '%s' must be %s, got %s", names[i].GetIdentifierName(),
			                      column->description, sql_types[i].ToString());
		}
		(*result).*(column->index) = i;
	}

	// Validate required columns
	if (result->node_index_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT NEWICK requires 'node_index' column");
	}
	if (result->parent_index_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT NEWICK requires 'parent_index' column");
	}
	// Parse parameters
	Value compression_param;
	const vector<Value> *compression_level_values = nullptr;
	bool edge_ids_specified = false;
	std::optional<std::string> placements_table;
	for (auto &option : input.info.options) {
		if (option.first == "edge_ids") {
			result->include_edge_ids = CopyOptionFlag(option.second);
			edge_ids_specified = true;
		} else if (option.first == "compression") {
			compression_param = option.second[0];
		} else if (option.first == "compression_level") {
			compression_level_values = &option.second;
		} else if (option.first == "placements") {
			placements_table = option.second[0].ToString();
			// Validate placements table/view exists and has correct schema
			ValidatePlacementTableSchema(context, placements_table.value());
		} else {
			ThrowUnreadCopyOption("NEWICK", option.first);
		}
	}

	// Use common compression detection (handles auto-detection from .gz extension)
	result->compression = DetectCompressionType(result->file_path, compression_param);
	if (compression_level_values) {
		result->compression_level =
		    GzipCompressionLevel(*compression_level_values, result->compression, result->file_path);
	}

	// Default: include edge_ids if edge_id column exists
	if (!edge_ids_specified) {
		result->include_edge_ids = (result->edge_id_idx != DConstants::INVALID_INDEX);
	}

	// Validate edge_id column if EDGE_IDS=true is explicitly requested
	if (edge_ids_specified && result->include_edge_ids && result->edge_id_idx == DConstants::INVALID_INDEX) {
		throw BinderException("EDGE_IDS option requires 'edge_id' column in input");
	}

	// Validate edge_id column if PLACEMENTS is specified (required for placement lookup)
	if (placements_table.has_value() && result->edge_id_idx == DConstants::INVALID_INDEX) {
		throw BinderException("PLACEMENTS option requires 'edge_id' column in input tree data");
	}

	result->placements_table = placements_table;

	return result;
}

//===--------------------------------------------------------------------===//
// Initialize Global
//===--------------------------------------------------------------------===//
static unique_ptr<GlobalFunctionData> NewickCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                                 const string &file_path) {
	auto &fdata = bind_data.Cast<NewickCopyBindData>();
	std::vector<miint::Placement> placements;
	if (fdata.placements_table) {
		placements = ReadPlacementTable(context, *fdata.placements_table);
	}
	return make_uniq<NewickCopyGlobalState>(file_path, fdata.compression, fdata.compression_level,
	                                        fdata.include_edge_ids, std::move(placements));
}

//===--------------------------------------------------------------------===//
// Initialize Local
//===--------------------------------------------------------------------===//
static unique_ptr<LocalFunctionData> NewickCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	return make_uniq<NewickCopyLocalState>();
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void NewickCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                           LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<NewickCopyBindData>();
	auto &lstate = lstate_p.Cast<NewickCopyLocalState>();

	// Get unified formats for all columns
	UnifiedVectorFormat node_index_data, parent_index_data;
	UnifiedVectorFormat name_data, branch_length_data, edge_id_data;

	input.data[fdata.node_index_idx].ToUnifiedFormat(node_index_data);
	input.data[fdata.parent_index_idx].ToUnifiedFormat(parent_index_data);

	auto node_indices = UnifiedVectorFormat::GetData<int64_t>(node_index_data);
	auto parent_indices = UnifiedVectorFormat::GetData<int64_t>(parent_index_data);

	bool has_name = fdata.name_idx != DConstants::INVALID_INDEX;
	bool has_branch_length = fdata.branch_length_idx != DConstants::INVALID_INDEX;
	bool has_edge_id = fdata.edge_id_idx != DConstants::INVALID_INDEX;

	const string_t *names_ptr = nullptr;
	const double *branch_lengths_ptr = nullptr;
	const int64_t *edge_ids_ptr = nullptr;

	if (has_name) {
		input.data[fdata.name_idx].ToUnifiedFormat(name_data);
		names_ptr = UnifiedVectorFormat::GetData<string_t>(name_data);
	}
	if (has_branch_length) {
		input.data[fdata.branch_length_idx].ToUnifiedFormat(branch_length_data);
		branch_lengths_ptr = UnifiedVectorFormat::GetData<double>(branch_length_data);
	}
	if (has_edge_id) {
		input.data[fdata.edge_id_idx].ToUnifiedFormat(edge_id_data);
		edge_ids_ptr = UnifiedVectorFormat::GetData<int64_t>(edge_id_data);
	}

	// Process each row
	lstate.local_nodes.reserve(lstate.local_nodes.size() + input.size());

	for (idx_t row = 0; row < input.size(); row++) {
		miint::NodeInput node;

		// node_index (required, must not be NULL)
		auto node_idx_row = node_index_data.sel->get_index(row);
		if (!node_index_data.validity.RowIsValid(node_idx_row)) {
			throw InvalidInputException("NULL value in node_index column (row " + std::to_string(row) + ")");
		}
		node.node_id = node_indices[node_idx_row];

		// parent_index (NULL means root)
		auto parent_idx_row = parent_index_data.sel->get_index(row);
		if (parent_index_data.validity.RowIsValid(parent_idx_row)) {
			node.parent_id = parent_indices[parent_idx_row];
		} else {
			node.parent_id = std::nullopt; // Root node
		}

		// name (optional)
		if (has_name) {
			auto name_row = name_data.sel->get_index(row);
			if (name_data.validity.RowIsValid(name_row)) {
				node.name = names_ptr[name_row].GetString();
			}
		}

		// branch_length (optional, default to NaN)
		node.branch_length = std::numeric_limits<double>::quiet_NaN();
		if (has_branch_length) {
			auto bl_row = branch_length_data.sel->get_index(row);
			if (branch_length_data.validity.RowIsValid(bl_row)) {
				double bl = branch_lengths_ptr[bl_row];
				if (bl < 0) {
					throw InvalidInputException("Negative branch_length " + std::to_string(bl) + " for node_index " +
					                            std::to_string(node.node_id));
				}
				node.branch_length = bl;
			}
		}

		// edge_id (optional)
		if (has_edge_id) {
			auto edge_row = edge_id_data.sel->get_index(row);
			if (edge_id_data.validity.RowIsValid(edge_row)) {
				node.edge_id = edge_ids_ptr[edge_row];
			}
		}

		lstate.local_nodes.push_back(std::move(node));
	}
}

//===--------------------------------------------------------------------===//
// Combine
//===--------------------------------------------------------------------===//
static void NewickCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                              LocalFunctionData &lstate_p) {
	auto &gstate = gstate_p.Cast<NewickCopyGlobalState>();
	auto &lstate = lstate_p.Cast<NewickCopyLocalState>();

	if (lstate.local_nodes.empty()) {
		return;
	}

	std::lock_guard<std::mutex> guard(gstate.lock);
	gstate.nodes.reserve(gstate.nodes.size() + lstate.local_nodes.size());
	for (auto &node : lstate.local_nodes) {
		gstate.nodes.push_back(std::move(node));
	}
	lstate.local_nodes.clear();
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void NewickCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<NewickCopyGlobalState>();

	if (gstate.nodes.empty()) {
		throw InvalidInputException("COPY FORMAT NEWICK: No data to write");
	}

	// Build the tree from collected nodes
	miint::NewickTree tree;
	try {
		tree = miint::NewickTree::build(gstate.nodes);
	} catch (const std::exception &e) {
		throw InvalidInputException("COPY FORMAT NEWICK: Failed to build tree: %s", e.what());
	}

	// Insert placements if any were provided (read when the COPY started)
	if (!gstate.placements.empty()) {
		// Validate tree has at least one non-NULL edge_id
		bool has_edge_id = false;
		for (uint32_t i = 0; i < tree.num_nodes(); i++) {
			if (tree.edge_id(i).has_value()) {
				has_edge_id = true;
				break;
			}
		}
		if (!has_edge_id) {
			throw InvalidInputException("COPY FORMAT NEWICK with PLACEMENTS: Tree has no edge_id values (all NULL). "
			                            "Placements require a tree with edge identifiers.");
		}

		// Insert placements into tree
		try {
			tree.insert_fully_resolved(gstate.placements);
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("COPY FORMAT NEWICK: Failed to insert placements: %s", e.what());
		}
	}

	// Clear edge_ids if not requested
	if (!gstate.include_edge_ids) {
		for (uint32_t i = 0; i < tree.num_nodes(); i++) {
			tree.set_edge_id(i, std::nullopt);
		}
	}

	// Serialize to Newick
	std::string newick = tree.to_newick();

	// Open file and write (CopyFileHandle handles compression)
	auto &fs = FileSystem::GetFileSystem(context);
	gstate.file = make_uniq<CopyFileHandle>(fs, gstate.file_path, gstate.compression, 1, gstate.compression_level);
	gstate.file->WriteString(newick);
	gstate.file->Close();
}

//===--------------------------------------------------------------------===//
// Initialize Operator: one tree, one file
//===--------------------------------------------------------------------===//
static void NewickCopyInitializeOperator(GlobalFunctionData &gstate, const PhysicalOperator &op) {
	auto &copy = op.Cast<PhysicalCopyToFile>();
	// Each file would get only the rows sent to it, part of the tree
	if (copy.per_thread_output || copy.partition_output) {
		throw InvalidInputException("COPY FORMAT NEWICK writes one tree to one file, so it takes no PER_THREAD_OUTPUT, "
		                            "PARTITION_BY or ORDER_BY");
	}
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
static void NewickCopyOptions(ClientContext &context, CopyOptionsInput &input) {
	DeclareCopyOption(input, "edge_ids", LogicalType::BOOLEAN);
	DeclareCopyOption(input, "compression", LogicalType::VARCHAR);
	DeclareCompressionLevelCopyOption(input);
	DeclareCopyOption(input, "placements", LogicalType::VARCHAR);
}

CopyFunction CopyNewickFunction::GetFunction() {
	CopyFunction func("newick");
	func.copy_options = NewickCopyOptions;
	func.copy_to_select = NewickCopySelect;
	func.copy_to_bind = NewickCopyBind;
	func.copy_to_initialize_local = NewickCopyInitializeLocal;
	func.copy_to_initialize_global = NewickCopyInitializeGlobal;
	func.copy_to_sink = NewickCopySink;
	func.copy_to_combine = NewickCopyCombine;
	func.copy_to_finalize = NewickCopyFinalize;
	func.initialize_operator = NewickCopyInitializeOperator;
	return func;
}

void CopyNewickFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
