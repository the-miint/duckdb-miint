#include "copy_gff.hpp"

#include "copy_option_parsing.hpp"
#include "predict_genes.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include "duckdb/parser/parser.hpp"
#include "duckdb/parser/query_node.hpp"
#include "duckdb/parser/statement/copy_statement.hpp"
#include "duckdb/parser/statement/select_statement.hpp"
#include "duckdb/planner/binder.hpp"

namespace duckdb {

namespace {

constexpr const char *kDialectGff3 = "gff3";
constexpr const char *kDialectProdigal = "prodigal";

// The CTE the user's query becomes, which the line macros read by name (query_table resolves a CTE).
constexpr const char *kInput = "__miint_gff_input";
// The same rows with a NULL genome_id, for Prodigal's GFF of rows without one.
constexpr const char *kInputWithGenome = "__miint_gff_input_genome";

// read_gff's columns, the input of DIALECT 'gff3'.
const vector<string> kGffColumns = {"seqid", "source", "type",  "position",  "stop_position",
                                    "score", "strand", "phase", "attributes"};

// predict_genes' columns that Prodigal's GFF writes (genome_id aside: without it, the rows are one genome).
const vector<string> kProdigalColumns = {"read_id",    "gene_index",      "position",        "stop_position",
                                         "strand",     "partial_left",    "partial_right",   "start_type",
                                         "rbs_motif",  "rbs_spacer",      "score",           "cscore",
                                         "sscore",     "rscore",          "uscore",          "tscore",
                                         "confidence", "gene_gc_content", "contig_length",   "translation_table",
                                         "uses_sd",    "model_desc",      "model_gc_content"};

// The COPY options FORMAT GFF passes on to the CSV writer as given. The others the CSV writer takes would break a
// GFF file: a second version line (APPEND, FILE_SIZE_BYTES, PER_THREAD_OUTPUT), columns that are not the line
// (WRITE_PARTITION_COLUMNS), or lines out of the relation's order (PRESERVE_ORDER).
const vector<string> kPassedOptions = {"compression",         "file_extension",   "use_tmp_file", "overwrite",
                                       "overwrite_or_ignore", "filename_pattern", "return_files"};

bool IsNullType(const LogicalType &t) {
	return t.id() == LogicalTypeId::SQLNULL;
}

bool IsPlainVarchar(const LogicalType &t) {
	return t.id() == LogicalTypeId::VARCHAR && !t.HasAlias();
}

bool IsTextMap(const LogicalType &t) {
	return t.id() == LogicalTypeId::MAP && IsPlainVarchar(MapType::KeyType(t)) && IsPlainVarchar(MapType::ValueType(t));
}

// The index of each wanted column among names, or a BinderException naming the first missing one.
vector<idx_t> FindColumns(const vector<string> &names, const vector<string> &wanted, const char *whose) {
	vector<idx_t> result;
	for (const auto &w : wanted) {
		idx_t found = DConstants::INVALID_INDEX;
		for (idx_t i = 0; i < names.size(); i++) {
			if (StringUtil::CIEquals(names[i], w)) {
				found = i;
				break;
			}
		}
		if (found == DConstants::INVALID_INDEX) {
			throw BinderException("COPY FORMAT GFF: the query has no %s column; it needs %s", w, whose);
		}
		result.push_back(found);
	}
	return result;
}

void RequireType(bool ok, const string &column, const LogicalType &got, const char *wanted) {
	if (!ok) {
		throw BinderException("COPY FORMAT GFF: %s must be %s (got %s)", column, wanted, got.ToString());
	}
}

// read_gff's column types, checked as the query is bound; the line macros then check the values.
void CheckGff3Types(const vector<string> &names, const vector<LogicalType> &types) {
	const auto cols = FindColumns(names, kGffColumns, "read_gff's columns");
	const auto type = [&](idx_t i) -> const LogicalType & {
		return types[cols[i]];
	};
	const auto text = [](const LogicalType &t) {
		return IsPlainVarchar(t) || IsNullType(t);
	};
	const auto integral = [](const LogicalType &t) {
		return t.IsIntegral() || IsNullType(t);
	};
	RequireType(text(type(0)), "seqid", type(0), "VARCHAR");
	RequireType(text(type(1)), "source", type(1), "VARCHAR");
	RequireType(text(type(2)), "type", type(2), "VARCHAR");
	RequireType(integral(type(3)), "position", type(3), "an integer");
	RequireType(integral(type(4)), "stop_position", type(4), "an integer");
	RequireType(type(5).IsNumeric() || IsNullType(type(5)), "score", type(5), "a number");
	RequireType(text(type(6)), "strand", type(6), "VARCHAR");
	RequireType(integral(type(7)), "phase", type(7), "an integer");
	RequireType(text(type(8)) || IsTextMap(type(8)), "attributes", type(8), "MAP(VARCHAR, VARCHAR) or VARCHAR");
}

bool HasColumn(const vector<string> &names, const char *name) {
	for (const auto &n : names) {
		if (StringUtil::CIEquals(n, name)) {
			return true;
		}
	}
	return false;
}

BoundStatement GffCopyPlan(Binder &binder, CopyStatement &stmt) {
	auto &info = *stmt.info;

	string dialect = kDialectGff3;
	bool dialect_set = false;
	string source = string("Prodigal_v") + kProdigalVersion;
	bool source_set = false;
	vector<Value> partition_by;
	identifier_map_t<vector<Value>> csv_options {{"file_extension", {Value("gff")}}};
	for (const auto &option : info.options) {
		const auto name = StringUtil::Lower(option.first.GetIdentifierName());
		if (name == "dialect") {
			dialect = StringUtil::Lower(CopyOptionValue(option.first, option.second).ToString());
			if (dialect != kDialectGff3 && dialect != kDialectProdigal) {
				throw BinderException("COPY FORMAT GFF: DIALECT must be '%s' or '%s' (got '%s')", kDialectGff3,
				                      kDialectProdigal, CopyOptionValue(option.first, option.second).ToString());
			}
			dialect_set = true;
		} else if (name == "source") {
			source = CopyOptionValue(option.first, option.second).ToString();
			source_set = true;
		} else if (name == "partition_by") {
			partition_by = option.second;
		} else if (std::find(kPassedOptions.begin(), kPassedOptions.end(), name) != kPassedOptions.end()) {
			csv_options[option.first] = option.second;
		} else {
			throw BinderException("COPY FORMAT GFF does not take %s: it would break the GFF file (a second version "
			                      "line, or columns or lines out of place)",
			                      StringUtil::Upper(name));
		}
	}

	// The query's columns, to check them here and name what is missing, rather than fail inside the macro.
	vector<string> names;
	vector<LogicalType> types;
	{
		auto node_copy = info.select_statement->Copy();
		auto child_binder = Binder::CreateBinder(binder.context, &binder);
		auto bound = child_binder->Bind(*node_copy);
		for (const auto &n : bound.names) {
			names.push_back(n.GetIdentifierName());
		}
		types = bound.types;
	}

	const bool prodigal = dialect == kDialectProdigal;
	if (!dialect_set && HasColumn(names, "gene_index")) {
		// predict_genes' rows are not GFF3 features; rather than fail on read_gff's missing columns, say what to add.
		throw BinderException("COPY FORMAT GFF: these look like predict_genes' columns; add DIALECT '%s' to write "
		                      "Prodigal's GFF",
		                      kDialectProdigal);
	}

	bool per_genome = false;
	if (!partition_by.empty()) {
		if (!prodigal || partition_by.size() != 1 || !StringUtil::CIEquals(partition_by[0].ToString(), "genome_id")) {
			throw BinderException("COPY FORMAT GFF: PARTITION_BY is (genome_id) with DIALECT '%s', a file per genome; "
			                      "a file per anything else would need its own version line",
			                      kDialectProdigal);
		}
		per_genome = true;
		csv_options["partition_by"] = partition_by;
	}

	string sql;
	if (!prodigal) {
		if (source_set) {
			throw BinderException("COPY FORMAT GFF: SOURCE is for DIALECT '%s'; read_gff's columns carry their own "
			                      "source",
			                      kDialectProdigal);
		}
		CheckGff3Types(names, types);
		sql = StringUtil::Format("WITH %s AS (SELECT 1) SELECT line FROM _miint_gff3_lines('%s')", kInput, kInput);
	} else {
		if (source.empty()) {
			throw BinderException("COPY FORMAT GFF: SOURCE is empty");
		}
		for (const auto &ch : source) {
			if (static_cast<unsigned char>(ch) < 0x20 || ch == 0x7f) {
				throw BinderException("COPY FORMAT GFF: SOURCE '%s' has a control character, which a GFF line cannot "
				                      "hold",
				                      source);
			}
		}
		FindColumns(names, kProdigalColumns, "predict_genes' columns");
		const bool has_genome = HasColumn(names, "genome_id");
		if (per_genome && !has_genome) {
			throw BinderException("COPY FORMAT GFF: PARTITION_BY (genome_id), but the query has no genome_id column");
		}
		// Without genome_id, the rows are one genome's.
		const string with_genome = has_genome
		                               ? ""
		                               : StringUtil::Format(", %s AS (SELECT NULL::VARCHAR AS genome_id, * FROM %s)",
		                                                    kInputWithGenome, kInput);
		sql = StringUtil::Format("WITH %s AS (SELECT 1)%s SELECT %sline FROM _miint_prodigal_gff_lines('%s', source := "
		                         "%s, per_genome := %s)",
		                         kInput, with_genome, per_genome ? "genome_id, " : "",
		                         has_genome ? kInput : kInputWithGenome, KeywordHelper::WriteQuoted(source, '\''),
		                         per_genome ? "true" : "false");
	}

	// The user's query becomes the input CTE.
	Parser parser(binder.context);
	parser.ParseQuery(sql);
	auto &select = parser.statements[0]->Cast<SelectStatement>();
	select.node->cte_map.map[Identifier(kInput)]->query_node = std::move(info.select_statement);
	info.select_statement = std::move(select.node);

	// Now the CSV writer writes the lines, as they are.
	info.format = "csv";
	info.options = std::move(csv_options);
	info.options["header"] = {Value::BOOLEAN(false)};
	info.options["quote"] = {Value("")};
	info.options["escape"] = {Value("")};
	return binder.Bind(stmt);
}

void GffCopyOptions(ClientContext &, CopyOptionsInput &input) {
	DeclareCopyOption(input, "compression", LogicalType::VARCHAR);
	DeclareCopyOption(input, "source", LogicalType::VARCHAR);
	DeclareCopyOption(input, "dialect", LogicalType::VARCHAR);
}

} // namespace

CopyFunction CopyGffFunction::GetFunction() {
	CopyFunction func("gff");
	func.copy_options = GffCopyOptions;
	func.plan = GffCopyPlan;
	func.extension = "gff";
	return func;
}

void CopyGffFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
