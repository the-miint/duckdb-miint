#include "copy_gff.hpp"
#include "copy_format_common.hpp"
#include "copy_option_parsing.hpp"
#include "predict_genes.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

#include <cstdio>
#include <mutex>
#include <string>
#include <vector>

namespace duckdb {

namespace {

// The program in column 2 of Prodigal's GFF, and in its "# Model Data" version, as the CLI names them.
const string kProdigalSource = string("Prodigal_v") + kProdigalVersion;
const string kProdigalModelVersion = string("Prodigal.v") + kProdigalVersion;

// read_gff's columns, the input of a GFF3 COPY, in order.
const vector<string> kGffColumns = {"seqid", "source", "type",  "position",  "stop_position",
                                    "score", "strand", "phase", "attributes"};
enum GffColumn : idx_t { SEQID, SOURCE, TYPE, POSITION, STOP_POSITION, SCORE, STRAND, PHASE, ATTRIBUTES };

// predict_genes' columns that Prodigal's GFF writes. genome_id may be missing (PARTITION_BY genome_id drops it from
// each file's columns); a contig is then told from the next by read_id alone.
const vector<string> kProdigalColumns = {"genome_id",       "read_id",       "gene_index",
                                         "position",        "stop_position", "strand",
                                         "partial_left",    "partial_right", "start_type",
                                         "rbs_motif",       "rbs_spacer",    "score",
                                         "cscore",          "sscore",        "rscore",
                                         "uscore",          "tscore",        "confidence",
                                         "gene_gc_content", "contig_length", "translation_table",
                                         "uses_sd",         "model_desc",    "model_gc_content"};
enum ProdigalColumn : idx_t {
	P_GENOME_ID,
	P_READ_ID,
	P_GENE_INDEX,
	P_POSITION,
	P_STOP_POSITION,
	P_STRAND,
	P_PARTIAL_LEFT,
	P_PARTIAL_RIGHT,
	P_START_TYPE,
	P_RBS_MOTIF,
	P_RBS_SPACER,
	P_SCORE,
	P_CSCORE,
	P_SSCORE,
	P_RSCORE,
	P_USCORE,
	P_TSCORE,
	P_CONFIDENCE,
	P_GENE_GC_CONTENT,
	P_CONTIG_LENGTH,
	P_TRANSLATION_TABLE,
	P_USES_SD,
	P_MODEL_DESC,
	P_MODEL_GC_CONTENT
};

bool IsIntegral(const LogicalType &t) {
	return t.IsIntegral() || t.id() == LogicalTypeId::SQLNULL;
}

bool IsNumber(const LogicalType &t) {
	return t.IsNumeric() || t.id() == LogicalTypeId::SQLNULL;
}

bool IsText(const LogicalType &t) {
	return t.id() == LogicalTypeId::VARCHAR || t.id() == LogicalTypeId::SQLNULL;
}

bool IsTextMap(const LogicalType &t) {
	return t.id() == LogicalTypeId::MAP && MapType::KeyType(t).id() == LogicalTypeId::VARCHAR &&
	       MapType::ValueType(t).id() == LogicalTypeId::VARCHAR;
}

struct GffCopyBindData : public FunctionData {
	string file_path;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	int compression_level = SEQUENCE_DEFAULT_GZIP_LEVEL;
	bool prodigal = false;
	string source = kProdigalSource;
	// Input column of each of kGffColumns / kProdigalColumns (INVALID_INDEX for a missing genome_id), and for GFF3 the
	// type the sink reads it as.
	vector<idx_t> columns;
	vector<LogicalType> read_types;

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<GffCopyBindData>(*this);
	}
	bool Equals(const FunctionData &other_p) const override {
		const auto &o = other_p.Cast<GffCopyBindData>();
		return file_path == o.file_path && compression == o.compression && compression_level == o.compression_level &&
		       prodigal == o.prodigal && source == o.source && columns == o.columns && read_types == o.read_types;
	}
};

// Each of `wanted` by name, or a BinderException naming the first missing one.
vector<idx_t> FindColumns(const vector<Identifier> &names, const vector<string> &wanted, const char *what) {
	vector<idx_t> out;
	for (const auto &w : wanted) {
		idx_t found = DConstants::INVALID_INDEX;
		for (idx_t i = 0; i < names.size(); i++) {
			if (StringUtil::CIEquals(names[i].GetIdentifierName(), w)) {
				found = i;
			}
		}
		if (found == DConstants::INVALID_INDEX) {
			throw BinderException("COPY FORMAT GFF: missing column '%s' (%s: %s)", w, what,
			                      StringUtil::Join(wanted, ", "));
		}
		out.push_back(found);
	}
	return out;
}

void RequireType(bool ok, const string &column, const LogicalType &got, const char *expected) {
	if (!ok) {
		throw BinderException("COPY FORMAT GFF: column '%s' must be %s (got %s)", column, expected, got.ToString());
	}
}

unique_ptr<FunctionData> GffCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                     const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<GffCopyBindData>();
	result->file_path = input.info.file_path;

	Value compression_param;
	const vector<Value> *compression_level_values = nullptr;
	bool source_set = false;
	bool dialect_set = false;
	for (auto &option : input.info.options) {
		if (option.first == "compression") {
			compression_param = option.second[0];
		} else if (option.first == "compression_level") {
			compression_level_values = &option.second;
		} else if (option.first == "source") {
			result->source = CopyOptionValue(option.first, option.second).ToString();
			source_set = true;
		} else if (option.first == "dialect") {
			const auto dialect = StringUtil::Lower(CopyOptionValue(option.first, option.second).ToString());
			if (dialect != "gff3" && dialect != "prodigal") {
				throw BinderException("COPY FORMAT GFF: DIALECT must be 'gff3' or 'prodigal' (got '%s')", dialect);
			}
			result->prodigal = dialect == "prodigal";
			dialect_set = true;
		} else {
			ThrowUnreadCopyOption("GFF", option.first);
		}
	}
	result->compression = DetectCompressionType(result->file_path, compression_param);
	if (compression_level_values) {
		result->compression_level =
		    GzipCompressionLevel(*compression_level_values, result->compression, result->file_path);
	}

	if (!dialect_set) {
		// predict_genes' rows are not GFF3 features; rather than fail on read_gff's missing columns, say what to add.
		for (const auto &name : names) {
			if (StringUtil::CIEquals(name.GetIdentifierName(), "gene_index")) {
				throw BinderException("COPY FORMAT GFF: these look like predict_genes' columns; add DIALECT 'prodigal' "
				                      "to write Prodigal's GFF");
			}
		}
	}

	if (!result->prodigal) {
		if (source_set) {
			throw BinderException("COPY FORMAT GFF: SOURCE is for DIALECT 'prodigal'; read_gff's columns carry their "
			                      "own source");
		}
		result->columns = FindColumns(names, kGffColumns, "read_gff's columns");
		const auto type = [&](GffColumn c) -> const LogicalType & {
			return sql_types[result->columns[c]];
		};
		RequireType(type(SEQID).id() == LogicalTypeId::VARCHAR, "seqid", type(SEQID), "VARCHAR");
		RequireType(type(SOURCE).id() == LogicalTypeId::VARCHAR, "source", type(SOURCE), "VARCHAR");
		RequireType(type(TYPE).id() == LogicalTypeId::VARCHAR, "type", type(TYPE), "VARCHAR");
		RequireType(type(POSITION).IsIntegral(), "position", type(POSITION), "an integer");
		RequireType(type(STOP_POSITION).IsIntegral(), "stop_position", type(STOP_POSITION), "an integer");
		RequireType(IsNumber(type(SCORE)), "score", type(SCORE), "a number");
		RequireType(IsText(type(STRAND)), "strand", type(STRAND), "VARCHAR");
		RequireType(IsIntegral(type(PHASE)), "phase", type(PHASE), "an integer");
		RequireType(IsText(type(ATTRIBUTES)) || IsTextMap(type(ATTRIBUTES)), "attributes", type(ATTRIBUTES),
		            "MAP(VARCHAR, VARCHAR) or VARCHAR");
		result->read_types = {LogicalType::VARCHAR,
		                      LogicalType::VARCHAR,
		                      LogicalType::VARCHAR,
		                      LogicalType::BIGINT,
		                      LogicalType::BIGINT,
		                      LogicalType::DOUBLE,
		                      LogicalType::VARCHAR,
		                      LogicalType::BIGINT,
		                      IsTextMap(type(ATTRIBUTES)) ? type(ATTRIBUTES) : LogicalType::VARCHAR};
	} else {
		for (const auto &ch : result->source) {
			if (static_cast<unsigned char>(ch) < 0x20 || ch == 0x7f) {
				throw BinderException("COPY FORMAT GFF: SOURCE '%s' has a control character, which a GFF line cannot "
				                      "hold",
				                      result->source);
			}
		}
		if (result->source.empty()) {
			throw BinderException("COPY FORMAT GFF: SOURCE is empty");
		}
		bool has_genome = false;
		for (const auto &name : names) {
			has_genome = has_genome || StringUtil::CIEquals(name.GetIdentifierName(), kProdigalColumns[P_GENOME_ID]);
		}
		vector<string> wanted = kProdigalColumns;
		if (!has_genome) {
			wanted[P_GENOME_ID] = kProdigalColumns[P_READ_ID]; // Stands in for it, so the indices stay.
		}
		result->columns = FindColumns(names, wanted, "predict_genes' columns");
		if (!has_genome) {
			result->columns[P_GENOME_ID] = DConstants::INVALID_INDEX;
		}
	}
	return std::move(result);
}

struct GffCopyGlobalState : public GlobalFunctionData {
	std::mutex lock;
	unique_ptr<CopyFileHandle> file;
	// Prodigal's GFF: the contig whose header lines were written last, and its number in the file.
	bool have_contig = false;
	Value genome;
	Value read_id;
	int64_t seqnum = 0;
};

unique_ptr<GlobalFunctionData> GffCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                       const string &file_path) {
	auto &data = bind_data.Cast<GffCopyBindData>();
	auto state = make_uniq<GffCopyGlobalState>();
	auto &fs = FileSystem::GetFileSystem(context);
	state->file = make_uniq<CopyFileHandle>(fs, file_path, data.compression, 1, data.compression_level);
	// Prodigal writes the directive with two spaces.
	state->file->WriteString(data.prodigal ? "##gff-version  3\n" : "##gff-version 3\n");
	return std::move(state);
}

unique_ptr<LocalFunctionData> GffCopyInitializeLocal(ExecutionContext &, FunctionData &) {
	return make_uniq<LocalFunctionData>();
}

// GFF3 reserves these in column 9 (and tab, newline, carriage return, '%' and control characters everywhere);
// read_gff decodes them.
void AppendEncoded(string &out, const string &text) {
	static const char *hex = "0123456789ABCDEF";
	for (unsigned char ch : text) {
		if (ch < 0x20 || ch == 0x7f || ch == '%' || ch == ';' || ch == '=' || ch == '&' || ch == ',') {
			out += '%';
			out += hex[ch >> 4];
			out += hex[ch & 0xf];
		} else {
			out += static_cast<char>(ch);
		}
	}
}

// A column 1-3 value: written as it is, since read_gff does not decode those columns; so a value that would break
// the line, or (a seqid starting with '#') turn it into a comment, is refused.
string FieldText(const Value &v, const char *column, idx_t row, bool first_column) {
	if (v.IsNull()) {
		throw InvalidInputException("COPY FORMAT GFF: NULL %s (row %llu)", column,
		                            static_cast<unsigned long long>(row));
	}
	const auto text = v.ToString();
	if (text.empty()) {
		throw InvalidInputException("COPY FORMAT GFF: empty %s (row %llu)", column,
		                            static_cast<unsigned long long>(row));
	}
	for (unsigned char ch : text) {
		if (ch < 0x20 || ch == 0x7f) {
			throw InvalidInputException("COPY FORMAT GFF: %s '%s' (row %llu) has a control character, which a GFF "
			                            "line cannot hold",
			                            column, text, static_cast<unsigned long long>(row));
		}
	}
	if (first_column && text[0] == '#') {
		throw InvalidInputException("COPY FORMAT GFF: %s '%s' (row %llu) starts with '#', which makes the line "
		                            "a comment",
		                            column, text, static_cast<unsigned long long>(row));
	}
	return text;
}

void AppendField(string &out, const Value &v, const char *column, idx_t row) {
	out += FieldText(v, column, row, string(column) == "seqid");
}

void WriteGffRows(ClientContext &context, const GffCopyBindData &data, GffCopyGlobalState &state, DataChunk &input) {
	DataChunk cast;
	cast.Initialize(context, data.read_types);
	cast.SetChildCardinality(input.size());
	for (idx_t c = 0; c < data.columns.size(); c++) {
		VectorOperations::Cast(context, input.data[data.columns[c]], cast.data[c], input.size());
	}
	string out;
	for (idx_t row = 0; row < input.size(); row++) {
		AppendField(out, cast.GetValue(SEQID, row), "seqid", row);
		out += '\t';
		AppendField(out, cast.GetValue(SOURCE, row), "source", row);
		out += '\t';
		AppendField(out, cast.GetValue(TYPE, row), "type", row);

		const Value position = cast.GetValue(POSITION, row);
		const Value stop = cast.GetValue(STOP_POSITION, row);
		if (position.IsNull() || BigIntValue::Get(position) < 1) {
			throw InvalidInputException("COPY FORMAT GFF: position must be at least 1 (row %llu has %s)",
			                            static_cast<unsigned long long>(row), position.ToString());
		}
		// stop_position is half-open; GFF's end is the last base, and a feature has at least one.
		if (stop.IsNull() || BigIntValue::Get(stop) <= BigIntValue::Get(position)) {
			throw InvalidInputException("COPY FORMAT GFF: stop_position must be greater than position (row %llu "
			                            "has %s and %s)",
			                            static_cast<unsigned long long>(row), position.ToString(), stop.ToString());
		}
		out += '\t' + std::to_string(BigIntValue::Get(position)) + '\t' + std::to_string(BigIntValue::Get(stop) - 1) +
		       '\t';

		const Value score = cast.GetValue(SCORE, row);
		out += score.IsNull() ? "." : score.ToString();
		out += '\t';

		const Value strand = cast.GetValue(STRAND, row);
		if (strand.IsNull()) {
			out += '.';
		} else {
			const auto &s = StringValue::Get(strand);
			if (s != "+" && s != "-" && s != "." && s != "?") {
				throw InvalidInputException("COPY FORMAT GFF: strand must be '+', '-', '.', '?' or NULL (row %llu has "
				                            "'%s')",
				                            static_cast<unsigned long long>(row), s);
			}
			out += s;
		}
		out += '\t';

		const Value phase = cast.GetValue(PHASE, row);
		if (phase.IsNull()) {
			out += '.';
		} else {
			const auto p = BigIntValue::Get(phase);
			if (p < 0 || p > 2) {
				throw InvalidInputException("COPY FORMAT GFF: phase must be 0, 1, 2 or NULL (row %llu has %lld)",
				                            static_cast<unsigned long long>(row), static_cast<long long>(p));
			}
			out += std::to_string(p);
		}
		out += '\t';

		const Value attributes = cast.GetValue(ATTRIBUTES, row);
		if (attributes.IsNull()) {
			out += '.';
		} else if (attributes.type().id() == LogicalTypeId::VARCHAR) {
			const auto &text = StringValue::Get(attributes);
			if (text.find_first_of("\t\n\r") != string::npos) {
				throw InvalidInputException("COPY FORMAT GFF: attributes text (row %llu) has a tab or line break; "
				                            "give attributes as a MAP to have them encoded",
				                            static_cast<unsigned long long>(row));
			}
			out += text.empty() ? "." : text;
		} else {
			const auto &entries = MapValue::GetChildren(attributes);
			if (entries.empty()) {
				out += '.';
			}
			for (idx_t i = 0; i < entries.size(); i++) {
				const auto &kv = StructValue::GetChildren(entries[i]);
				const auto key = kv[0].ToString();
				if (kv[1].IsNull()) {
					throw InvalidInputException("COPY FORMAT GFF: attribute '%s' (row %llu) is NULL; GFF has no "
					                            "attribute without a value",
					                            key, static_cast<unsigned long long>(row));
				}
				if (i > 0) {
					out += ';';
				}
				AppendEncoded(out, key);
				out += '=';
				AppendEncoded(out, StringValue::Get(kv[1]));
			}
		}
		out += '\n';
	}
	state.file->WriteString(out);
}

// printf, as Prodigal formats its numbers.
template <typename... ARGS>
string Format(const char *fmt, ARGS... args) {
	char buf[512];
	const int n = std::snprintf(buf, sizeof(buf), fmt, args...);
	return string(buf, static_cast<size_t>(MaxValue(n, 0)));
}

void WriteProdigalRows(const GffCopyBindData &data, GffCopyGlobalState &state, DataChunk &input) {
	const auto value = [&](ProdigalColumn c, idx_t row) {
		return data.columns[c] == DConstants::INVALID_INDEX ? Value() : input.GetValue(data.columns[c], row);
	};
	// A value the line needs: NULL is refused, rather than written as "NULL" or read as 0.
	const auto required = [&](ProdigalColumn c, idx_t row) {
		Value v = value(c, row);
		if (v.IsNull()) {
			throw InvalidInputException("COPY FORMAT GFF: NULL %s (row %llu)", kProdigalColumns[c],
			                            static_cast<unsigned long long>(row));
		}
		return v;
	};
	// Text written inside a line, NULL as `none`.
	const auto text = [&](ProdigalColumn c, idx_t row, const char *none) {
		const Value v = value(c, row);
		return v.IsNull() ? string(none) : FieldText(v, kProdigalColumns[c].c_str(), row, false);
	};
	string out;
	for (idx_t row = 0; row < input.size(); row++) {
		const Value gene_index = value(P_GENE_INDEX, row);
		const Value read_id = value(P_READ_ID, row);
		if (read_id.IsNull() && gene_index.IsNull()) {
			continue; // A genome's status row (training failed, or no genes) has no contig.
		}
		const string seqid = FieldText(read_id, "read_id", row, true);
		const Value genome = value(P_GENOME_ID, row);
		if (!state.have_contig || !Value::NotDistinctFrom(genome, state.genome) ||
		    !Value::NotDistinctFrom(read_id, state.read_id)) {
			state.have_contig = true;
			state.genome = genome;
			state.read_id = read_id;
			state.seqnum++;
			const string model_desc = text(P_MODEL_DESC, row, "");
			out +=
			    Format("# Sequence Data: seqnum=%lld;seqlen=%lld;seqhdr=\"%s\"\n", static_cast<long long>(state.seqnum),
			           static_cast<long long>(required(P_CONTIG_LENGTH, row).GetValue<int64_t>()), seqid.c_str());
			out += Format("# Model Data: version=%s;run_type=%s;model=\"%s\";gc_cont=%.2f;transl_table=%d;uses_sd=%d\n",
			              kProdigalModelVersion.c_str(), model_desc.empty() ? "Single" : "Metagenomic",
			              model_desc.empty() ? "Ab initio" : model_desc.c_str(),
			              required(P_MODEL_GC_CONTENT, row).GetValue<double>() * 100.0,
			              required(P_TRANSLATION_TABLE, row).GetValue<int32_t>(),
			              required(P_USES_SD, row).GetValue<bool>() ? 1 : 0);
		}
		if (gene_index.IsNull()) {
			continue; // A contig without genes: its header lines only, as the CLI writes it.
		}
		const string strand = required(P_STRAND, row).ToString();
		if (strand != "+" && strand != "-") {
			throw InvalidInputException("COPY FORMAT GFF: strand must be '+' or '-' (row %llu has '%s')",
			                            static_cast<unsigned long long>(row), strand);
		}
		const double score = required(P_SCORE, row).GetValue<double>();
		out += seqid + '\t' + data.source + "\tCDS\t" + std::to_string(required(P_POSITION, row).GetValue<int64_t>()) +
		       '\t' + std::to_string(required(P_STOP_POSITION, row).GetValue<int64_t>() - 1) + '\t' +
		       Format("%.1f", score) + '\t' + strand + "\t0\t";
		out += Format("ID=%lld_%d;partial=%d%d;", static_cast<long long>(state.seqnum), gene_index.GetValue<int32_t>(),
		              required(P_PARTIAL_LEFT, row).GetValue<bool>() ? 1 : 0,
		              required(P_PARTIAL_RIGHT, row).GetValue<bool>() ? 1 : 0);
		out += "start_type=" + FieldText(required(P_START_TYPE, row), "start_type", row, false) +
		       ";rbs_motif=" + text(P_RBS_MOTIF, row, "None") + ";rbs_spacer=" + text(P_RBS_SPACER, row, "None");
		out +=
		    Format(";gc_cont=%.3f;conf=%.2f;score=%.2f;cscore=%.2f;sscore=%.2f;rscore=%.2f;uscore=%.2f;tscore=%.2f;\n",
		           required(P_GENE_GC_CONTENT, row).GetValue<double>(), required(P_CONFIDENCE, row).GetValue<double>(),
		           score, required(P_CSCORE, row).GetValue<double>(), required(P_SSCORE, row).GetValue<double>(),
		           required(P_RSCORE, row).GetValue<double>(), required(P_USCORE, row).GetValue<double>(),
		           required(P_TSCORE, row).GetValue<double>());
	}
	state.file->WriteString(out);
}

void GffCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p, LocalFunctionData &,
                 DataChunk &input) {
	auto &data = bind_data.Cast<GffCopyBindData>();
	auto &state = gstate_p.Cast<GffCopyGlobalState>();
	// No execution_mode, so COPY sinks on one thread, in the query's order: Prodigal's header lines rely on it.
	// PER_THREAD_OUTPUT (which COPY handles before bind, so it cannot be refused here) gives each thread a file of its
	// own; docs/writing.md says not to use it with predict_genes' columns.
	std::lock_guard<std::mutex> guard(state.lock);
	if (data.prodigal) {
		WriteProdigalRows(data, state, input);
	} else {
		WriteGffRows(context.client, data, state, input);
	}
}

void GffCopyCombine(ExecutionContext &, FunctionData &, GlobalFunctionData &, LocalFunctionData &) {
}

void GffCopyFinalize(ClientContext &, FunctionData &, GlobalFunctionData &gstate_p) {
	auto &state = gstate_p.Cast<GffCopyGlobalState>();
	state.file->Close();
}

void GffCopyOptions(ClientContext &, CopyOptionsInput &input) {
	DeclareCopyOption(input, "compression", LogicalType::VARCHAR);
	DeclareCompressionLevelCopyOption(input);
	DeclareCopyOption(input, "source", LogicalType::VARCHAR);
	DeclareCopyOption(input, "dialect", LogicalType::VARCHAR);
}

} // namespace

CopyFunction CopyGffFunction::GetFunction() {
	CopyFunction func("gff");
	func.copy_options = GffCopyOptions;
	func.copy_to_bind = GffCopyBind;
	func.copy_to_initialize_local = GffCopyInitializeLocal;
	func.copy_to_initialize_global = GffCopyInitializeGlobal;
	func.copy_to_sink = GffCopySink;
	func.copy_to_combine = GffCopyCombine;
	func.copy_to_finalize = GffCopyFinalize;
	func.extension = "gff";
	return func;
}

void CopyGffFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
