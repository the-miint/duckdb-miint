#include "copy_newick.hpp"
#include "copy_format_common.hpp"
#include "NewickTree.hpp"
#include "placement_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct NewickCopyBindData : public FunctionData {
	std::string file_path;
	bool include_edge_ids = false;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;

	// Placements are read at bind time (not during finalize) to avoid deadlocks
	// when executing queries during an active COPY operation.
	// This is consistent with SQL snapshot semantics - placement data is captured
	// when the query is planned.
	std::vector<miint::Placement> placements;

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
		result->placements = placements;
		result->node_index_idx = node_index_idx;
		result->parent_index_idx = parent_index_idx;
		result->name_idx = name_idx;
		result->branch_length_idx = branch_length_idx;
		result->edge_id_idx = edge_id_idx;
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<NewickCopyBindData>();
		// Note: we don't compare placements for equality - if all other fields match,
		// we consider the bind data equal (placements are derived from the table reference)
		return file_path == other.file_path && include_edge_ids == other.include_edge_ids &&
		       compression == other.compression && placements.size() == other.placements.size() &&
		       node_index_idx == other.node_index_idx && parent_index_idx == other.parent_index_idx &&
		       name_idx == other.name_idx && branch_length_idx == other.branch_length_idx &&
		       edge_id_idx == other.edge_id_idx;
	}
};

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct NewickCopyGlobalState : public GlobalFunctionData {
	std::mutex lock;
	std::string file_path;
	FileCompressionType compression;
	bool include_edge_ids;

	// Placements are copied from bind data (read at bind time)
	std::vector<miint::Placement> placements;

	// File handle (opened in finalize, not at init)
	unique_ptr<CopyFileHandle> file;

	// Collected node data from all threads
	std::vector<miint::NodeInput> nodes;

	NewickCopyGlobalState(const std::string &path, FileCompressionType comp, bool edge_ids,
	                      std::vector<miint::Placement> placements_p)
	    : file_path(path), compression(comp), include_edge_ids(edge_ids), placements(std::move(placements_p)) {
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
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> NewickCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                               const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<NewickCopyBindData>();
	result->file_path = input.info.file_path;

	// Find column indices
	for (idx_t i = 0; i < names.size(); i++) {
		const auto &name = names[i];
		if (StringUtil::CIEquals(name, "node_index")) {
			result->node_index_idx = i;
		} else if (StringUtil::CIEquals(name, "parent_index")) {
			result->parent_index_idx = i;
		} else if (StringUtil::CIEquals(name, "name")) {
			result->name_idx = i;
		} else if (StringUtil::CIEquals(name, "branch_length")) {
			result->branch_length_idx = i;
		} else if (StringUtil::CIEquals(name, "edge_id")) {
			result->edge_id_idx = i;
		}
		// Ignore unknown columns (like is_tip, filepath)
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
	bool edge_ids_specified = false;
	std::optional<std::string> placements_table;
	for (auto &option : input.info.options) {
		auto loption = StringUtil::Lower(option.first);
		if (loption == "edge_ids") {
			if (option.second.size() != 1 || option.second[0].type().id() != LogicalTypeId::BOOLEAN) {
				throw BinderException("EDGE_IDS option requires a boolean value");
			}
			result->include_edge_ids = BooleanValue::Get(option.second[0]);
			edge_ids_specified = true;
		} else if (loption == "compression") {
			if (option.second.size() != 1 || option.second[0].type().id() != LogicalTypeId::VARCHAR) {
				throw BinderException("COMPRESSION option requires a string value");
			}
			compression_param = option.second[0];
		} else if (loption == "placements") {
			if (option.second.size() != 1 || option.second[0].type().id() != LogicalTypeId::VARCHAR) {
				throw BinderException("PLACEMENTS option requires a string value (table or view name)");
			}
			placements_table = option.second[0].ToString();
			// Validate placements table/view exists and has correct schema
			ValidatePlacementTableSchema(context, placements_table.value());
		} else {
			throw BinderException("Unknown option for COPY FORMAT NEWICK: %s", option.first);
		}
	}

	// Use common compression detection (handles auto-detection from .gz extension)
	result->compression = DetectCompressionType(result->file_path, compression_param);

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

	// Read placements at bind time to avoid deadlocks during COPY execution.
	// This is safe because bind happens before query execution starts.
	if (placements_table.has_value()) {
		result->placements = ReadPlacementTable(context, placements_table.value());
	}

	return result;
}

//===--------------------------------------------------------------------===//
// Initialize Global
//===--------------------------------------------------------------------===//
static unique_ptr<GlobalFunctionData> NewickCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                                 const string &file_path) {
	auto &fdata = bind_data.Cast<NewickCopyBindData>();
	return make_uniq<NewickCopyGlobalState>(file_path, fdata.compression, fdata.include_edge_ids,
	                                        fdata.placements);
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

	input.data[fdata.node_index_idx].ToUnifiedFormat(input.size(), node_index_data);
	input.data[fdata.parent_index_idx].ToUnifiedFormat(input.size(), parent_index_data);

	auto node_indices = UnifiedVectorFormat::GetData<int64_t>(node_index_data);
	auto parent_indices = UnifiedVectorFormat::GetData<int64_t>(parent_index_data);

	bool has_name = fdata.name_idx != DConstants::INVALID_INDEX;
	bool has_branch_length = fdata.branch_length_idx != DConstants::INVALID_INDEX;
	bool has_edge_id = fdata.edge_id_idx != DConstants::INVALID_INDEX;

	const string_t *names_ptr = nullptr;
	const double *branch_lengths_ptr = nullptr;
	const int64_t *edge_ids_ptr = nullptr;

	if (has_name) {
		input.data[fdata.name_idx].ToUnifiedFormat(input.size(), name_data);
		names_ptr = UnifiedVectorFormat::GetData<string_t>(name_data);
	}
	if (has_branch_length) {
		input.data[fdata.branch_length_idx].ToUnifiedFormat(input.size(), branch_length_data);
		branch_lengths_ptr = UnifiedVectorFormat::GetData<double>(branch_length_data);
	}
	if (has_edge_id) {
		input.data[fdata.edge_id_idx].ToUnifiedFormat(input.size(), edge_id_data);
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
					throw InvalidInputException("Negative branch_length " + std::to_string(bl) +
					                            " for node_index " + std::to_string(node.node_id));
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

	// Insert placements if any were provided (read at bind time)
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
			throw InvalidInputException(
			    "COPY FORMAT NEWICK with PLACEMENTS: Tree has no edge_id values (all NULL). "
			    "Placements require a tree with edge identifiers.");
		}

		// Insert placements into tree (placements were already read at bind time)
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
	gstate.file = make_uniq<CopyFileHandle>(fs, gstate.file_path, gstate.compression);
	gstate.file->WriteString(newick);
	gstate.file->Close();
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
CopyFunction CopyNewickFunction::GetFunction() {
	CopyFunction func("newick");
	func.copy_to_bind = NewickCopyBind;
	func.copy_to_initialize_local = NewickCopyInitializeLocal;
	func.copy_to_initialize_global = NewickCopyInitializeGlobal;
	func.copy_to_sink = NewickCopySink;
	func.copy_to_combine = NewickCopyCombine;
	func.copy_to_finalize = NewickCopyFinalize;
	return func;
}

void CopyNewickFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
