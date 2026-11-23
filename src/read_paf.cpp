#include "read_paf.hpp"
#include "miint_function_reference.hpp"

#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

#include <stdexcept>

namespace duckdb {

namespace {

using miint::TagType;

// The tags minimap2 lists for PAF (minimap2(1), "OUTPUT FORMAT"), in its order, as typed
// columns named the way read_alignments names its tag columns, except SA, which minimap2
// writes only in SAM. Any other tag stays in the `tags` map.
struct TypedTag {
	const char *key;
	TagType type;
};

constexpr TypedTag MINIMAP2_TAGS[] = {
    {"tp", TagType::PrintableChar}, {"cm", TagType::Integer},       {"s1", TagType::Integer}, {"s2", TagType::Integer},
    {"NM", TagType::Integer},       {"MD", TagType::String},        {"AS", TagType::Integer}, {"ms", TagType::Integer},
    {"nn", TagType::Integer},       {"ts", TagType::PrintableChar}, {"cg", TagType::String},  {"cs", TagType::String},
    {"dv", TagType::Float},         {"de", TagType::Float},         {"rl", TagType::Integer}, {"zd", TagType::Integer},
};
constexpr idx_t TYPED_TAG_COUNT = sizeof(MINIMAP2_TAGS) / sizeof(MINIMAP2_TAGS[0]);

enum Column : idx_t {
	READ_ID,
	QUERY_LENGTH,
	QUERY_START,
	QUERY_STOP,
	STRAND,
	REFERENCE,
	REFERENCE_LENGTH,
	POSITION,
	STOP_POSITION,
	RESIDUE_MATCHES,
	BLOCK_LENGTH,
	MAPQ,
	FIRST_TYPED_TAG
};
constexpr idx_t TAGS = FIRST_TYPED_TAG + TYPED_TAG_COUNT;

std::string ColumnName(const TypedTag &tag) {
	return "tag_" + StringUtil::Lower(tag.key);
}

LogicalType ColumnType(TagType type) {
	switch (type) {
	case TagType::Integer:
		return LogicalType::BIGINT;
	case TagType::Float:
		return LogicalType::DOUBLE;
	default:
		return LogicalType::VARCHAR;
	}
}

// The index of `key` in MINIMAP2_TAGS, or TYPED_TAG_COUNT for any other tag. Every minimap2
// key is two characters, so a longer one, like each of myloasm's, is turned away at once.
idx_t TypedTagIndex(const std::string &key) {
	if (key.size() != 2) {
		return TYPED_TAG_COUNT;
	}
	for (idx_t i = 0; i < TYPED_TAG_COUNT; i++) {
		if (key[0] == MINIMAP2_TAGS[i].key[0] && key[1] == MINIMAP2_TAGS[i].key[1]) {
			return i;
		}
	}
	return TYPED_TAG_COUNT;
}

// Move minimap2's tags out of `tags` into `typed`, each checked against its column's type,
// and leave the rest in `tags`. Throws std::runtime_error naming the tag; ScanLines adds the
// file and line.
void TakeTypedTags(std::vector<miint::SamStyleTag> &tags, std::vector<ReadPafTableFunction::TypedTagValue> &typed) {
	for (auto &value : typed) {
		value.present = false;
	}
	size_t kept = 0;
	for (size_t i = 0; i < tags.size(); i++) {
		auto &tag = tags[i];
		const idx_t found = TypedTagIndex(tag.key);
		if (found == TYPED_TAG_COUNT) {
			if (kept != i) {
				tags[kept] = std::move(tag);
			}
			kept++;
			continue;
		}
		const TypedTag &spec = MINIMAP2_TAGS[found];
		const char letter = static_cast<char>(spec.type);
		if (tag.type) {
			if (*tag.type != spec.type) {
				throw std::runtime_error("tag '" + tag.key + "' has type '" + static_cast<char>(*tag.type) +
				                         "'; read_paf reads it into " + ColumnName(spec) + " as type '" + letter + "'");
			}
		} else if (!miint::IsValidTagValue(spec.type, tag.value)) {
			// An untyped field's value has not been checked against any type yet.
			throw std::runtime_error("'" + tag.key + ":" + tag.value + "' has a value that is not valid for type '" +
			                         letter + "', the type of " + ColumnName(spec));
		}
		auto &out = typed[found];
		const string_t text(tag.value.data(), static_cast<uint32_t>(tag.value.size()));
		switch (spec.type) {
		case TagType::Integer:
			// The value matches the integer grammar, so the cast fails only on overflow.
			if (!TryCast::Operation<string_t, int64_t>(text, out.integer, false)) {
				throw std::runtime_error("tag '" + tag.key + "' value '" + tag.value + "' does not fit " +
				                         ColumnName(spec) + ", a 64-bit integer");
			}
			break;
		case TagType::Float:
			// DuckDB's cast reads every form the float grammar accepts, one too large for a DOUBLE
			// as an infinity (as htslib's strtod does), so this fails only if that changes.
			if (!TryCast::Operation<string_t, double>(text, out.real, false)) {
				throw std::runtime_error("tag '" + tag.key + "' value '" + tag.value + "' is not a number " +
				                         ColumnName(spec) + " can hold");
			}
			break;
		default:
			out.text = std::move(tag.value);
		}
		out.present = true;
	}
	tags.resize(kept);
}

void SetString(Vector &vec, idx_t row, const std::string &value) {
	FlatVector::GetDataMutable<string_t>(vec)[row] = StringVector::AddString(vec, value);
}

void SetInt64(Vector &vec, idx_t row, int64_t value) {
	FlatVector::GetDataMutable<int64_t>(vec)[row] = value;
}

void SetTypedTag(Vector &vec, idx_t row, TagType type, const ReadPafTableFunction::TypedTagValue &value) {
	if (!value.present) {
		FlatVector::ValidityMutable(vec).SetInvalid(row);
		return;
	}
	switch (type) {
	case TagType::Integer:
		SetInt64(vec, row, value.integer);
		break;
	case TagType::Float:
		FlatVector::GetDataMutable<double>(vec)[row] = value.real;
		break;
	default:
		SetString(vec, row, value.text);
	}
}

} // namespace

unique_ptr<FunctionData> ReadPafTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<Identifier> &names) {
	names = {"read_id",          "query_length", "query_start",   "query_stop",      "strand",       "reference",
	         "reference_length", "position",     "stop_position", "residue_matches", "block_length", "mapq"};
	return_types = {LogicalType::VARCHAR, LogicalType::BIGINT,  LogicalType::BIGINT, LogicalType::BIGINT,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::BIGINT,
	                LogicalType::BIGINT,  LogicalType::BIGINT,  LogicalType::BIGINT, LogicalType::UTINYINT};
	for (const auto &tag : MINIMAP2_TAGS) {
		names.emplace_back(ColumnName(tag));
		return_types.push_back(ColumnType(tag.type));
	}
	names.push_back("tags");
	return_types.push_back(LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR));
	return BindLineScan(context, input, names, return_types);
}

unique_ptr<LocalTableFunctionState> ReadPafTableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	auto state = make_uniq<LocalState>();
	state->typed_tags.resize(TYPED_TAG_COUNT);
	return std::move(state);
}

void ReadPafTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &local = data_p.local_state->Cast<LocalState>();
	const auto &r = local.record;
	ScanLines(
	    context, data_p, output,
	    [&](std::string_view line) {
		    local.parser.Parse(line, local.record);
		    TakeTypedTags(local.record.tags, local.typed_tags);
		    return true;
	    },
	    [&](idx_t row) {
		    SetString(output.data[READ_ID], row, r.read_id);
		    SetInt64(output.data[QUERY_LENGTH], row, r.query_length);
		    SetInt64(output.data[QUERY_START], row, r.query_start);
		    SetInt64(output.data[QUERY_STOP], row, r.query_stop);
		    SetString(output.data[STRAND], row, r.strand);
		    SetString(output.data[REFERENCE], row, r.reference);
		    SetInt64(output.data[REFERENCE_LENGTH], row, r.reference_length);
		    SetInt64(output.data[POSITION], row, r.position);
		    SetInt64(output.data[STOP_POSITION], row, r.stop_position);
		    SetInt64(output.data[RESIDUE_MATCHES], row, r.residue_matches);
		    SetInt64(output.data[BLOCK_LENGTH], row, r.block_length);
		    FlatVector::GetDataMutable<uint8_t>(output.data[MAPQ])[row] = r.mapq;
		    for (idx_t i = 0; i < TYPED_TAG_COUNT; i++) {
			    SetTypedTag(output.data[FIRST_TYPED_TAG + i], row, MINIMAP2_TAGS[i].type, local.typed_tags[i]);
		    }
		    SetTagMap(output.data[TAGS], row, r.tags);
	    });
}

void ReadPafTableFunction::Register(ExtensionLoader &loader) {
	auto tf = MakeLineScanFunction("read_paf", Execute, Bind, InitLocal);
	tf.GetSignature().GetParameter(0).SetName("filename");
	CreateTableFunctionInfo info(std::move(tf));
	info.on_conflict = OnCreateConflict::ALTER_ON_CONFLICT;
	info.descriptions.push_back(::miint::reference::read_paf());
	loader.RegisterFunction(std::move(info));
}

} // namespace duckdb
