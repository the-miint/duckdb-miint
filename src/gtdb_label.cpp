#include "gtdb_label.hpp"

#include "table_function_common.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/common/vector/struct_vector.hpp"
#include "duckdb/function/scalar_function.hpp"

#include <optional>
#include <string_view>

namespace duckdb {

namespace {

// The characters of a GTDB name: ASCII letters, digits and '-'. StringUtil's character
// tests are ASCII-only; <cctype>'s follow the locale.
bool IsNameChar(char c) {
	return StringUtil::CharacterIsAlphaNumeric(c) || c == '-';
}

bool IsName(std::string_view s) {
	if (s.empty()) {
		return false;
	}
	for (const char c : s) {
		if (!IsNameChar(c)) {
			return false;
		}
	}
	return true;
}

// A species epithet: words of name characters and '.', separated by single spaces, as in
// "bar subsp. baz".
bool IsEpithet(std::string_view s) {
	if (s.empty() || s.front() == ' ' || s.back() == ' ') {
		return false;
	}
	for (size_t i = 0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			if (s[i - 1] == ' ') {
				return false;
			}
		} else if (!IsNameChar(s[i]) && s[i] != '.') {
			return false;
		}
	}
	return true;
}

// A suffix, without its leading underscore: runs of capitals and digits joined by single
// underscores, so "Y_123" is one suffix and "A_" and "_A" are none. Capitals are tested by
// range: StringUtil has no single-character test for them, and std::isupper follows the
// locale.
bool IsSuffix(std::string_view s) {
	if (s.empty() || s.front() == '_' || s.back() == '_') {
		return false;
	}
	for (size_t i = 0; i < s.size(); ++i) {
		const char c = s[i];
		if (c == '_') {
			if (s[i - 1] == '_') {
				return false;
			}
		} else if (!((c >= 'A' && c <= 'Z') || StringUtil::CharacterIsDigit(c))) {
			return false;
		}
	}
	return true;
}

struct Labelled {
	std::string_view base;
	std::optional<std::string_view> suffix;
};

// Everything before the first '_' is the base, which `is_base` must accept; everything
// after it is the suffix. No '_' means no suffix.
std::optional<Labelled> SplitSuffix(std::string_view s, bool (*is_base)(std::string_view)) {
	const auto underscore = s.find('_');
	const auto base = s.substr(0, underscore);
	if (!is_base(base)) {
		return std::nullopt;
	}
	if (underscore == std::string_view::npos) {
		return Labelled {base, std::nullopt};
	}
	const auto suffix = s.substr(underscore + 1);
	if (!IsSuffix(suffix)) {
		return std::nullopt;
	}
	return Labelled {base, suffix};
}

void SetPart(Vector &part, idx_t row, std::optional<std::string_view> value) {
	if (value) {
		FlatVector::GetDataMutable<string_t>(part)[row] = StringVector::AddString(part, value->data(), value->size());
	} else {
		FlatVector::SetNull(part, row, true);
	}
}

// Split each row's name with `split`, which writes the parts of a name that fits and
// returns false for one that does not. A NULL name, or one that does not fit, is a NULL
// row, its parts included (SetStructRowNull).
template <class SPLIT>
void ExecuteSplit(DataChunk &args, Vector &result, SPLIT split) {
	UnifiedVectorFormat name_data;
	args.data[0].ToUnifiedFormat(name_data);
	auto names = UnifiedVectorFormat::GetData<string_t>(name_data);
	auto &parts = StructVector::GetEntries(result);
	for (idx_t row = 0; row < args.size(); ++row) {
		const auto idx = name_data.sel->get_index(row);
		if (!name_data.validity.RowIsValid(idx) ||
		    !split(std::string_view(names[idx].GetData(), names[idx].GetSize()), parts, row)) {
			SetStructRowNull(result, row);
		}
	}
}

void GtdbLabelExecute(DataChunk &args, ExpressionState &, Vector &result) {
	ExecuteSplit(args, result, [](std::string_view name, vector<Vector> &parts, idx_t row) {
		const auto label = SplitSuffix(name, IsName);
		if (!label) {
			return false;
		}
		SetPart(parts[0], row, label->base);
		SetPart(parts[1], row, label->suffix);
		return true;
	});
}

// The genus label, one space, then the epithet with its own optional suffix, which can
// only sit at the epithet's end since an epithet holds no '_'.
void GtdbSpeciesLabelExecute(DataChunk &args, ExpressionState &, Vector &result) {
	ExecuteSplit(args, result, [](std::string_view name, vector<Vector> &parts, idx_t row) {
		const auto space = name.find(' ');
		if (space == std::string_view::npos) {
			return false;
		}
		const auto genus = SplitSuffix(name.substr(0, space), IsName);
		const auto epithet = SplitSuffix(name.substr(space + 1), IsEpithet);
		if (!genus || !epithet) {
			return false;
		}
		SetPart(parts[0], row, genus->base);
		SetPart(parts[1], row, genus->suffix);
		SetPart(parts[2], row, epithet->base);
		SetPart(parts[3], row, epithet->suffix);
		return true;
	});
}

LogicalType VarcharStruct(std::initializer_list<const char *> names) {
	child_list_t<LogicalType> fields;
	for (const auto *name : names) {
		fields.emplace_back(name, LogicalType::VARCHAR);
	}
	return LogicalType::STRUCT(std::move(fields));
}

} // namespace

void GtdbLabelFunctions::Register(ExtensionLoader &loader) {
	const FunctionParameter name(Identifier("name"), LogicalType::VARCHAR);
	loader.RegisterFunction(ScalarFunction("gtdb_label", {name}, VarcharStruct({"base", "suffix"}), GtdbLabelExecute));
	loader.RegisterFunction(ScalarFunction("gtdb_species_label", {name},
	                                       VarcharStruct({"genus_base", "genus_suffix", "epithet", "epithet_suffix"}),
	                                       GtdbSpeciesLabelExecute));
}

} // namespace duckdb
