#include "krepp_index_create.hpp"

#include "KreppPlacer.hpp"
#include "NewickTree.hpp"
#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "tree_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "miint_named_parameter.hpp"
#include "miint_streaming_query.hpp"

namespace duckdb {

namespace {

constexpr const char *kCallerName = "krepp_index_create";

// Upper bound for the `threads` parameter. See the check in the binder for
// what goes wrong without one and how this number was chosen.
constexpr uint32_t kMaxThreads = 256;

// The characters NewickTree::to_newick quotes a label for. A reference name
// carrying any of them would reach krepp quoted in the tree and bare in the
// input map, and krepp matches the two by string equality - so the reference
// would simply never be visited (build_for_subtree walks the tree and looks
// each tip up in the map). No k-mers, no message: exactly the silent shortfall
// this function exists to prevent.
// krepp decides what its input file IS from the first non-whitespace byte of
// the whole file: '>' or '@' means FASTA/FASTQ, anything else means the TSV map
// (ext/krepp/src/index.cpp, read_input_file). The map's first line starts with a
// reference name, so a reference that sorts first and begins with '>' would flip
// krepp into per-sequence mode, where a guide tree is refused outright - a
// confusing error about a tree the caller did supply.
bool NameConfusesKreppInputSniffer(const std::string &name) {
	return !name.empty() && (name.front() == '>' || name.front() == '@');
}

bool NameNeedsNewickQuoting(const std::string &name) {
	for (unsigned char c : name) {
		if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';' || c == '{' || c == '}' || c == '\'' || c == '"' ||
		    c == '[' || c == ']' || std::isspace(c)) {
			return true;
		}
	}
	return false;
}

// Scratch space for the per-reference FASTAs and the input map. Under $TMPDIR
// (default /tmp), one directory per call, removed when the build finishes or
// throws. Mirrors align_bowtie2's MakeTempIndexDir.
std::string MakeTempWorkDir() {
	const char *tmp = std::getenv("TMPDIR");
	if (!tmp || !*tmp) {
		tmp = "/tmp";
	}
	std::string tmpl = std::string(tmp) + "/miint-krepp-XXXXXX";
	std::vector<char> buf(tmpl.begin(), tmpl.end());
	buf.push_back('\0');
	if (::mkdtemp(buf.data()) == nullptr) {
		throw IOException("%s: failed to create a temp directory under '%s' (errno=%d)", kCallerName, std::string(tmp),
		                  errno);
	}
	return std::string(buf.data());
}

struct TempWorkDir {
	std::string path;

	explicit TempWorkDir(std::string path_p) : path(std::move(path_p)) {
	}
	~TempWorkDir() {
		std::error_code ec;
		std::filesystem::remove_all(path, ec);
		// `ec` ignored: a destructor cannot propagate it, and the only remedy
		// beyond letting the OS reclaim the space would be to throw here.
	}
	TempWorkDir(const TempWorkDir &) = delete;
	TempWorkDir &operator=(const TempWorkDir &) = delete;
};

// Streams the sequence relation out to one FASTA per reference.
//
// Exactly one file is open at a time. Rows sharing a read_id are almost always
// adjacent (a genome's contigs arrive together), so that is one open per
// reference; an interleaved relation still produces the right files, just with
// more reopens. The alternative - a handle per reference - is what makes 196k
// genomes a descriptor problem, and holding the corpus in memory to sort it is
// what makes it a memory problem.
//
// Files are named by ordinal, not by reference name: the name is the caller's
// data and may contain anything a filesystem does not allow.
class ReferenceFastaWriter {
public:
	explicit ReferenceFastaWriter(std::string dir) : dir_(std::move(dir)) {
	}

	void Add(const std::string &name, const std::string &sequence) {
		size_t idx;
		auto it = index_of_.find(name);
		const bool is_new = it == index_of_.end();
		if (is_new) {
			idx = names_.size();
			index_of_.emplace(name, idx);
			names_.push_back(name);
			paths_.push_back(dir_ + "/ref_" + std::to_string(idx) + ".fa");
			records_.push_back(0);
		} else {
			idx = it->second;
		}

		if (!out_.is_open() || open_idx_ != idx) {
			Close();
			// A reference seen for the first time truncates; one revisited
			// appends, so an interleaved relation does not lose its earlier
			// records.
			out_.open(paths_[idx], is_new ? std::ios::out : std::ios::app);
			if (!out_) {
				throw IOException("%s: failed to open '%s' for writing", kCallerName, paths_[idx]);
			}
			open_idx_ = idx;
		}

		records_[idx]++;
		out_ << '>' << names_[idx] << '_' << records_[idx] << '\n' << sequence << '\n';
		if (!out_) {
			throw IOException("%s: failed writing sequences for reference '%s' to '%s'", kCallerName, names_[idx],
			                  paths_[idx]);
		}
	}

	void Close() {
		if (!out_.is_open()) {
			return;
		}
		out_.close();
		// Checked, unlike a bare close(): the last buffered bytes of a
		// reference are flushed here, and the next open() calls clear() and
		// would wipe the failbit before anything looked at it.
		if (!out_) {
			throw IOException("%s: failed to finish writing '%s'", kCallerName, paths_[open_idx_]);
		}
	}

	// name<TAB>path, one line per reference, in first-seen order. krepp reads
	// this with getline + two tab-delimited fields and rejects a duplicate name;
	// there are none here because the names come from a map.
	std::string WriteInputMap() {
		Close();
		const std::string map_path = dir_ + "/input_map.tsv";
		std::ofstream map(map_path);
		for (size_t i = 0; i < names_.size(); ++i) {
			map << names_[i] << '\t' << paths_[i] << '\n';
		}
		map.close();
		if (!map) {
			throw IOException("%s: failed to write the reference map '%s'", kCallerName, map_path);
		}
		return map_path;
	}

	const std::vector<std::string> &names() const {
		return names_;
	}

private:
	std::string dir_;
	std::vector<std::string> names_;
	std::vector<std::string> paths_;
	std::vector<uint64_t> records_;
	std::unordered_map<std::string, size_t> index_of_;
	std::ofstream out_;
	size_t open_idx_ = 0;
};

// Reads an optional named parameter, leaving `target` alone when absent.
// Same shape as place_krepp's, kept local for the same reason.
template <typename T>
void ReadOptional(const named_argument_map_t &params, const char *key, T &target, T (*convert)(const Value &)) {
	auto it = params.find(key);
	if (it != params.end() && !it->second.IsNull()) {
		target = convert(it->second);
	}
}

// k, w and h are uint8_t in krepp's config. Taking them as INTEGER and
// narrowing silently would turn `k := 285` into k = 29 - a valid-looking build
// with the wrong k - so the range is checked before the cast rather than after.
uint8_t ReadByteParam(const named_argument_map_t &params, const char *key, uint8_t fallback, bool *was_set = nullptr) {
	auto it = params.find(key);
	if (it == params.end() || it->second.IsNull()) {
		return fallback;
	}
	const int64_t value = it->second.GetValue<int64_t>();
	if (value < 0 || value > 255) {
		throw BinderException("%s: %s must be between 0 and 255 (got %lld)", kCallerName, std::string(key),
		                      static_cast<long long>(value));
	}
	if (was_set) {
		*was_set = true;
	}
	return static_cast<uint8_t>(value);
}

uint32_t ReadUIntParam(const named_argument_map_t &params, const char *key, uint32_t fallback) {
	auto it = params.find(key);
	if (it == params.end() || it->second.IsNull()) {
		return fallback;
	}
	const int64_t value = it->second.GetValue<int64_t>();
	// Every caller declares its parameter as LogicalType::INTEGER, so DuckDB has
	// already cast to INT32 and rejected anything outside that before Bind runs:
	// `threads := 2147483648` arrives as an Invalid Input Error from the cast and
	// never reaches here. A negative value does reach here, though, because INT32
	// holds it. The old message named 0..UINT32_MAX, advertising a ceiling this
	// function cannot actually be handed.
	if (value < 0) {
		throw BinderException("%s: %s must not be negative (got %lld)", kCallerName, std::string(key),
		                      static_cast<long long>(value));
	}
	if (value > static_cast<int64_t>(UINT32_MAX)) {
		// Unreachable while every caller is INTEGER. Kept so that declaring one
		// BIGINT later truncates loudly instead of silently wrapping.
		throw BinderException("%s: %s must be at most %u (got %lld)", kCallerName, std::string(key), UINT32_MAX,
		                      static_cast<long long>(value));
	}
	return static_cast<uint32_t>(value);
}

// k, w and h of one partial, read back from disk rather than recomputed so the
// row reports the index that exists, including the w and h krepp derived from k
// when they were not given.
//
// Reads the BINARY `metadata<suffix>`, not the `.txt` sidecar beside it. The
// sidecar is advisory - krepp writes it second, synthesises it when absent, and
// never compares it - so parsing it here meant reporting values krepp might not
// be using. ValidateIndexLayout needs the same three fields, so the read lives
// in krepp_detail and both callers get the same answer.
void ReadIndexMetadata(const std::string &index_dir, const std::string &suffix, int32_t &k, int32_t &w, int32_t &h) {
	try {
		miint::krepp_detail::ReadPartialConfig(index_dir + "/metadata" + suffix, suffix, k, w, h);
	} catch (const std::exception &e) {
		// Fail rather than report zeros next to status='ok'. These three columns
		// are the only way a caller learns what w and h krepp derived, so a
		// silent 0 would be a wrong answer wearing a success badge.
		throw IOException("%s: %s", kCallerName, std::string(e.what()));
	}
}

} // namespace

unique_ptr<FunctionData> KreppIndexCreateTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<Identifier> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 2) {
		throw BinderException("%s requires sequence_table and output_path parameters", kCallerName);
	}
	data->sequence_table = input.inputs[0].ToString();
	RejectCTERelationName(input, data->sequence_table);
	data->output_path = input.inputs[1].ToString();

	ReadOptional<std::string>(input.named_parameters, "tree_table", data->tree_table,
	                          [](const Value &v) { return v.ToString(); });
	if (!data->tree_table.empty()) {
		RejectCTERelationName(input, data->tree_table);
	}
	ReadOptional<std::string>(input.named_parameters, "newick_path", data->newick_path,
	                          [](const Value &v) { return v.ToString(); });

	// A backbone is required, not optional. krepp will happily build without
	// one - it generates a tree from the reference names and then skips writing
	// any tree file at all (index.cpp: "Skipped saving a backbone for the
	// index!"). place_krepp cannot use the result without being handed a
	// separate Newick, so an index built that way here would be a trap.
	if (data->tree_table.empty() == data->newick_path.empty()) {
		throw BinderException("%s requires exactly one of tree_table or newick_path; the index needs a backbone tree "
		                      "to be usable with place_krepp",
		                      kCallerName);
	}

	bool w_set = false;
	bool h_set = false;
	data->options.k = ReadByteParam(input.named_parameters, "k", data->options.k);
	const uint8_t w = ReadByteParam(input.named_parameters, "w", 0, &w_set);
	const uint8_t h = ReadByteParam(input.named_parameters, "h", 0, &h_set);
	if (w_set) {
		data->options.w = w;
	}
	if (h_set) {
		data->options.h = h;
	}
	data->options.m = ReadUIntParam(input.named_parameters, "m", data->options.m);
	data->options.r = ReadUIntParam(input.named_parameters, "r", data->options.r);
	// krepp's CLI bounds -m with CLI::PositiveNumber; the library path this uses
	// does not, and BaseLSH::validate_configuration checks only w, h and k. So
	// m = 0 reaches set_nrows(), which computes `hash_size % m` - division by
	// zero, undefined behaviour whose symptom depends on the CPU. Measured on
	// arm64 (AArch64 UDIV returns 0 rather than trapping) it does not crash: it
	// builds an index with zero k-mers at every node and reports status='ok',
	// which is worse than the SIGFPE the same code would raise on x86_64.
	// r is deliberately not bounded below - krepp's CLI allows r = 0
	// (CLI::NonNegativeNumber) and nothing divides by it.
	if (data->options.m < 1) {
		throw BinderException("%s: m must be at least 1 (got %u); krepp partitions the LSH space modulo m and "
		                      "divides by it",
		                      kCallerName, data->options.m);
	}
	data->options.threads = ReadUIntParam(input.named_parameters, "threads", data->options.threads);
	if (data->options.threads < 1) {
		throw BinderException("%s: threads must be at least 1 (got %u)", kCallerName, data->options.threads);
	}
	// Refused rather than clamped: a build with krepp's OpenMP regions compiled
	// out would take the number, run on one core anyway, and report success -
	// and the only way to tell would be the wall clock.
	// Nothing else bounds `threads` from above: set_num_threads does not clamp
	// (ext/krepp/src/common.cpp:28) and neither does krepp's CLI, so the value
	// reaches omp_set_num_threads unchanged. libomp does not fail that call when
	// it cannot create the threads - it calls abort(), taking the whole DuckDB
	// process down with no SQL error and no chance to catch it:
	//   OMP: Error #34: System unable to allocate necessary resources for OMP thread
	// The ceiling is not a fixed property of the machine. Measured here: a bare
	// process built a team of 8192 and aborted at 16384, while inside DuckDB -
	// whose own pool has already taken threads - 1000 succeeded and 10000
	// aborted. So this is a guard rail set far below the lowest abort seen
	// anywhere, not a tuned value; no index build benefits from oversubscribing
	// this hard, and a query must not be able to kill the server.
	if (data->options.threads > kMaxThreads) {
		throw BinderException("%s: threads must be at most %u (got %u)", kCallerName, kMaxThreads,
		                      data->options.threads);
	}
	if (data->options.threads > 1 && !miint::KreppIndexThreadsSupported()) {
		throw BinderException("%s: threads := %u needs an OpenMP runtime and this build has none, so krepp's "
		                      "index regions were compiled out. Rebuild with libomp available (brew install "
		                      "libomp, or -DMIINT_LIBOMP_PREFIX=/prefix), or pass threads := 1",
		                      kCallerName, data->options.threads);
	}
	data->options.sdust_t = ReadUIntParam(input.named_parameters, "sdust_t", data->options.sdust_t);
	data->options.sdust_w = ReadUIntParam(input.named_parameters, "sdust_w", data->options.sdust_w);
	ReadOptional<bool>(input.named_parameters, "frac", data->options.frac,
	                   [](const Value &v) { return BooleanValue::Get(v); });

	// Same contract as the aligners and save_bowtie2_index: read_id plus
	// sequence1, BIGINT ids allowed. The id is stringified in InitGlobal, since
	// a krepp reference name is text either way.
	data->schema = ValidateSequenceTableSchema(context, data->sequence_table, /*allow_bigint=*/true);
	if (!data->tree_table.empty()) {
		ValidateTreeTableSchema(context, data->tree_table);
	}

	return_types = data->types;
	names = StringsToIdentifiers(data->names);
	return std::move(data);
}

// The filename suffix krepp will give every file this build writes. The shape
// itself lives in krepp_detail::PartialSuffix, which mirrors the IndexMultiple
// constructor (ext/krepp/src/index.cpp:249-251); this only supplies the options.
std::string PartialSuffixFor(const miint::KreppIndexOptions &options) {
	return miint::krepp_detail::PartialSuffix(options.m, options.r, options.frac);
}

// Decide whether `output_path` is somewhere this build may write.
//
// Absent or empty: yes. Already holding partials of the SAME index at a
// different residue: yes - that is the multi-partial build, and adding to it is
// the point. It exists only under frac := false (see PartialHashConfig).
// Holding a different index, this exact residue, a frac := true partial, or
// debris: no, and we neither delete nor overwrite any of it.
void CheckOutputPathAcceptsPartial(const KreppIndexCreateTableFunction::Data &data) {
	std::error_code ec;
	if (!std::filesystem::exists(data.output_path, ec)) {
		return;
	}
	// is_directory before is_empty: is_empty() answers "size == 0" for a regular
	// file, so an empty file here would return early as though it were an empty
	// directory, and the mistake would only surface much later inside krepp's
	// create_directory - after the temporary FASTAs had already been written.
	if (!std::filesystem::is_directory(data.output_path, ec)) {
		throw IOException("%s: output_path '%s' exists and is not a directory", kCallerName, data.output_path);
	}
	if (std::filesystem::is_empty(std::filesystem::path(data.output_path), ec)) {
		return;
	}

	std::map<std::string, std::set<std::string>> existing;
	try {
		existing = miint::krepp_detail::ValidateIndexLayout(data.output_path);
	} catch (const std::exception &e) {
		// Debris from a build that died partway, files beside no complete index,
		// or two indexes already mixed together. None of those is a directory to
		// add a partial to, and none of them is ours to clean up. Files that are
		// not index files (DiscoverPartials skips them) do not cause this beside a
		// complete index: they are left alone and the build goes ahead.
		throw IOException("%s: output_path '%s' is not empty and does not hold one complete krepp index for this "
		                  "build to add a partial to (%s). krepp writes alongside whatever is already there, so "
		                  "building here would leave a mixed directory - remove it or choose another path",
		                  kCallerName, data.output_path, std::string(e.what()));
	}

	const std::string ours = PartialSuffixFor(data.options);
	if (existing.find(ours) != existing.end()) {
		// A frac := true index is a single partial (see PartialHashConfig), so a
		// different r is no way out for it.
		throw IOException("%s: output_path '%s' already holds the partial this build would write ('%s'). krepp does "
		                  "not clear what is there, so this would write over it. %s",
		                  kCallerName, data.output_path, ours,
		                  data.options.frac ? "Use another path" : "Use a different r, or another path");
	}
	const std::string our_config = miint::krepp_detail::PartialHashConfig(ours);
	const std::string their_config = miint::krepp_detail::PartialHashConfig(existing.begin()->first);
	if (our_config != their_config) {
		throw IOException("%s: output_path '%s' holds an index built with a different hash configuration "
		                  "('%s' vs this build's '%s'); krepp would load them as one index and reject them. "
		                  "Partials of one index must share m and frac, and differ only in r, which frac := true "
		                  "does not allow: a frac := true index is a single partial",
		                  kCallerName, data.output_path, their_config.substr(1), our_config.substr(1));
	}

	// m and frac are in the filename; k, w and h are not, and they decide which
	// k-mers get indexed and where the LSH positions land. Compare them against
	// a partial that is already there.
	//
	// Compare the EFFECTIVE values, not the ones the caller typed. Unset w and h
	// are derived from k (ext/krepp/src/index.cpp:237-238), so what has to agree
	// is what krepp will actually use - the partial on disk may have been built
	// with an explicit h that a later job leaves unset. Gating these on
	// has_value() let exactly that through: same k, one partial at h := 10 and
	// the next with h unset (k - 16), accepted here and then rejected by krepp
	// at query time, long after the build. w is worse - LSHF::check_compatible
	// (ext/krepp/src/lshf.cpp:163-170) compares m, h, k, frac, r under
	// frac := true, and the positions, but never w, so a w mismatch is caught
	// nowhere downstream and just leaves one index holding two different sets
	// of minimizers.
	//
	// The positions themselves are drawn from `gen`, which BuildKreppIndex
	// reseeds to a fixed state before every build - that is what makes two
	// separately built partials share a hash function at all.
	int32_t their_k = 0, their_w = 0, their_h = 0;
	ReadIndexMetadata(data.output_path, existing.begin()->first, their_k, their_w, their_h);
	// Narrowing to uint8_t before widening mirrors krepp: its k, w and h are all
	// uint8_t (ext/krepp/src/index.hpp:107-109), so a k below 16 wraps there too
	// and the prediction has to wrap with it.
	const auto our_w = static_cast<int32_t>(data.options.w.value_or(static_cast<uint8_t>(data.options.k + 6)));
	const auto our_h = static_cast<int32_t>(data.options.h.value_or(static_cast<uint8_t>(data.options.k - 16)));
	const char *which = nullptr;
	int32_t theirs = 0, mine = 0;
	if (their_k != static_cast<int32_t>(data.options.k)) {
		which = "k";
		theirs = their_k;
		mine = static_cast<int32_t>(data.options.k);
	} else if (their_w != our_w) {
		which = "w";
		theirs = their_w;
		mine = our_w;
	} else if (their_h != our_h) {
		which = "h";
		theirs = their_h;
		mine = our_h;
	}
	if (which != nullptr) {
		throw IOException("%s: output_path '%s' holds partials built with %s := %d, but this build uses %s := %d. "
		                  "Every partial of one index must agree on k, w, h, m and frac, and differ only in r, "
		                  "which frac := true does not allow",
		                  kCallerName, data.output_path, which, theirs, which, mine);
	}
}

unique_ptr<GlobalTableFunctionState> KreppIndexCreateTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// krepp creates its output directory and writes files whose names encode the
	// resolved config; it never clears what is already there. That is a problem
	// for a SECOND INDEX in the same directory - krepp would load both as one -
	// but it is exactly how a large index is meant to be built: one residue per
	// job, every job pointed at the same directory, `-r` varying and everything
	// else fixed. Refusing every non-empty directory blocked the second case to
	// prevent the first. So check which one it is, and never delete anything.
	CheckOutputPathAcceptsPartial(data);

	TempWorkDir work(MakeTempWorkDir());

	// ---- The backbone -----------------------------------------------------
	//
	// Both paths end with the same two things: the Newick text krepp will read,
	// and the tip names, which are the only thing that decides whether a
	// reference is indexed at all. The tree object itself is not kept - a
	// backbone at Greengenes2 scale is hundreds of megabytes, and the text has
	// to be held anyway to be written out.
	std::string newick_text;
	std::string newick_origin;
	std::unordered_set<std::string> tip_names;
	if (!data.tree_table.empty()) {
		newick_origin = "tree_table '" + data.tree_table + "'";
		// ReadTreeTable throws on an empty relation before returning, so there
		// is no empty case to handle here (src/tree_table_reader.cpp).
		const auto nodes = ReadTreeTable(context, data.tree_table);
		const auto tree = miint::NewickTree::build(nodes);
		newick_text = tree.to_newick();
		const auto names_v = tree.tip_names();
		tip_names.insert(names_v.begin(), names_v.end());
	} else {
		newick_origin = "newick_path '" + data.newick_path + "'";
		std::ifstream in(data.newick_path, std::ios::binary);
		if (!in) {
			throw IOException("%s: cannot read %s", kCallerName, newick_origin);
		}
		newick_text.assign((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		miint::krepp_detail::ValidateNewickLexically(newick_text, newick_origin);
		// Parsed only here: the tree_table path already has the tree.
		const auto tree = miint::NewickTree::parse(newick_text);
		const auto names_v = tree.tip_names();
		tip_names.insert(names_v.begin(), names_v.end());
	}

	// The generated path gets the same lexical check the file path already had.
	// Belt and braces rather than a real check - to_newick writes one line and
	// quotes what has to be quoted - but a linear scan is nothing next to
	// indexing the references, and it is the only thing standing between a
	// future change in to_newick and a tree krepp would exit on.
	if (!data.tree_table.empty()) {
		miint::krepp_detail::ValidateNewickLexically(newick_text, newick_origin);
	}

	const std::string newick_file = work.path + "/backbone.nwk";
	{
		std::ofstream out(newick_file);
		out << newick_text;
		if (newick_text.empty() || newick_text.back() != '\n') {
			out << '\n';
		}
		out.close();
		if (!out) {
			throw IOException("%s: failed to write the backbone tree to '%s'", kCallerName, newick_file);
		}
	}

	// ---- The references ---------------------------------------------------
	ReferenceFastaWriter writer(work.path);
	{
		// Deliberately not QuerySequenceStream, which is the house helper for
		// this shape: every row here is validated against its own reference
		// name and written straight to disk, so there is nothing to gain from
		// the SequenceRecordBatch it accumulates in between. Revisit if a
		// second caller wants the same per-row rejections.
		auto conn = MakeReadOnlyHelperConnection(context);
		// DuckDB pauses the producer of a streaming result once the chunks
		// buffered for Fetch reach streaming_buffer_size - all of it, or a fixed
		// share of it when the plan uses the batched collector - and sizes each
		// chunk with Vector::GetAllocationSize, which is 16 bytes a row for VARCHAR
		// however long the strings are (duckdb/src/common/types/vector.cpp:895). At
		// the default 976.5 KiB, with sequence_table a view over Parquet, peak RSS
		// while the FASTA were written grew with the sequence bytes: 2.05 GiB for
		// 400 random 4 Mbp genomes, 3.73 GiB for 800 (m := 64, r := 0,
		// frac := false, threads := 8). At 1,000 bytes it was 0.79 and 0.81 GiB,
		// against 0.76 GiB for a plain scan of the same Parquet file; the index
		// files other than the timestamped metadata .txt were byte-identical, and
		// krepp's build started after 51 s and 101 s instead of 44 s and 89 s. As a
		// read_fastx view over the same 800 genomes, the feed peaked at 3.13 GiB
		// uncapped and 0.46 GiB capped, against 0.43 GiB for a plain scan of the
		// view, with the same index files other than the metadata .txt. Written into
		// this connection's ClientConfig rather than through SET, which
		// lock_configuration refuses (DBConfig::CheckLock); no other connection
		// sees it.
		ClientConfig::GetConfig(*conn.context).streaming_buffer_size = 1000;
		// Cast to VARCHAR so a BIGINT read_id and a text one reach krepp the
		// same way; a reference name is text on both sides of the map.
		// sequence2 is selected only to refuse it. A krepp reference is one
		// sequence; handed a paired-end relation this would otherwise index
		// half the data and report status='ok'. Same contract as
		// ReadSubjectTable: the column may exist, its values may not.
		const std::string sql = "SELECT read_id::VARCHAR AS read_id, sequence1" +
		                        std::string(data.schema.has_sequence2 ? ", sequence2" : "") + " FROM " +
		                        KeywordHelper::WriteOptionallyQuoted(data.sequence_table);
		auto result = SubmitStream(conn, sql);
		if (result->HasError()) {
			throw InvalidInputException("%s: failed to read '%s': %s", kCallerName, data.sequence_table,
			                            result->GetError());
		}
		while (auto chunk = result->Fetch()) {
			for (idx_t i = 0; i < chunk->size(); i++) {
				const auto id_value = chunk->GetValue(0, i);
				const auto seq_value = chunk->GetValue(1, i);
				if (id_value.IsNull() || seq_value.IsNull()) {
					throw InvalidInputException("%s: '%s' has a NULL read_id or sequence1; every reference needs a "
					                            "name and a sequence",
					                            kCallerName, data.sequence_table);
				}
				if (data.schema.has_sequence2 && !chunk->GetValue(2, i).IsNull()) {
					throw InvalidInputException("%s: '%s' has a non-NULL sequence2; a krepp reference is a single "
					                            "sequence, and indexing only sequence1 would silently drop half the "
					                            "data",
					                            kCallerName, data.sequence_table);
				}
				const auto name = id_value.GetValue<std::string>();
				const auto sequence = seq_value.GetValue<std::string>();
				if (name.empty()) {
					throw InvalidInputException("%s: '%s' has an empty read_id", kCallerName, data.sequence_table);
				}
				if (sequence.empty()) {
					throw InvalidInputException("%s: reference '%s' has an empty sequence1", kCallerName, name);
				}
				// The same alphabet place_krepp holds query sequences to, for
				// two reasons rather than one. A byte above 127 indexes past
				// the end of krepp's 128-entry nucleotide table. And a newline
				// followed by '>' would open a SECOND FASTA record inside this
				// reference's file - which krepp does not reject or even
				// notice: DynHT::fill_table loops over every record in the file
				// and folds them all into the same leaf, so the reference would
				// silently carry k-mers from content no row claimed.
				const size_t bad = sequence.find_first_not_of(miint::krepp_detail::kNucleotideAlphabet);
				if (bad != std::string::npos) {
					throw InvalidInputException(
					    "%s: sequence for reference '%s' contains a character that is not a nucleotide code: byte %d "
					    "at offset %lld",
					    kCallerName, name, static_cast<int32_t>(static_cast<unsigned char>(sequence[bad])),
					    static_cast<long long>(bad));
				}
				if (NameConfusesKreppInputSniffer(name)) {
					throw InvalidInputException("%s: reference name '%s' starts with '>' or '@'; krepp reads the "
					                            "first byte of its reference map to tell a map from a FASTA, so a "
					                            "name like this can turn the whole build into per-sequence mode",
					                            kCallerName, name);
				}
				if (NameNeedsNewickQuoting(name)) {
					throw InvalidInputException("%s: reference name '%s' contains a character Newick has to quote "
					                            "(whitespace or one of ()[]{},:;'\"), so it could never match a tree "
					                            "tip; rename it in sequence_table",
					                            kCallerName, name);
				}
				writer.Add(name, sequence);
			}
		}
	}
	if (writer.names().empty()) {
		throw InvalidInputException("%s: '%s' produced no sequences", kCallerName, data.sequence_table);
	}

	// ---- Names on both sides have to agree --------------------------------
	//
	// krepp walks the tree and looks each tip up in the map. A reference the
	// tree does not name is never visited and never reported - it just is not
	// in the index. A tip with no reference is reported, but only to stderr,
	// which nothing in a SQL session sees.
	int64_t unmatched_count = 0;
	vector<std::string> unmatched_sample;
	for (const auto &name : writer.names()) {
		if (tip_names.find(name) == tip_names.end()) {
			unmatched_count++;
			if (unmatched_sample.size() < 5) {
				unmatched_sample.push_back(name);
			}
		}
	}
	if (unmatched_count > 0) {
		throw InvalidInputException("%s: %lld reference name(s) in '%s' are not tips of %s (e.g. %s). krepp indexes "
		                            "by walking the tree, so those references would be silently absent from the "
		                            "index; shear the tree or fix the names first",
		                            kCallerName, static_cast<long long>(unmatched_count), data.sequence_table,
		                            newick_origin, StringUtil::Join(unmatched_sample, ", "));
	}
	const int64_t num_references = static_cast<int64_t>(writer.names().size());
	const int64_t tips_without_references = static_cast<int64_t>(tip_names.size()) - num_references;
	if (tips_without_references > 0) {
		miint::EmitWarning(context, std::string(kCallerName) + ": " + std::to_string(tips_without_references) + " of " +
		                                std::to_string(tip_names.size()) + " backbone tips have no sequence in '" +
		                                data.sequence_table +
		                                "' and are skipped by krepp; placements can still land on their edges");
	}

	// ---- Build ------------------------------------------------------------
	auto options = data.options;
	options.index_dir = data.output_path;
	options.input_map_path = writer.WriteInputMap();
	options.newick_path = newick_file;
	// krepp's own create_directory only makes the last component, and throws a
	// std::filesystem error rather than going through error_exit when a parent
	// is missing. Create the parents here, as save_bowtie2_index does.
	const std::filesystem::path index_path(data.output_path);
	if (index_path.has_parent_path()) {
		std::error_code parent_ec;
		std::filesystem::create_directories(index_path.parent_path(), parent_ec);
		if (parent_ec) {
			throw IOException("%s: failed to create directory '%s': %s", kCallerName, index_path.parent_path().string(),
			                  parent_ec.message());
		}
	}
	try {
		miint::BuildKreppIndex(options);
	} catch (const miint::KreppFatalError &e) {
		// Whatever krepp had already written to output_path stays there; it is
		// the caller's directory and nothing here deletes it. Say so, because
		// the retry will hit the not-empty check above.
		throw IOException("%s: krepp rejected the build: %s (any files already written to '%s' are left in place)",
		                  kCallerName, std::string(e.what()), data.output_path);
	} catch (const std::filesystem::filesystem_error &e) {
		throw IOException("%s: krepp could not write to '%s': %s", kCallerName, data.output_path,
		                  std::string(e.what()));
	}

	// ---- What actually landed on disk -------------------------------------
	std::map<std::string, std::set<std::string>> partials;
	try {
		partials = miint::krepp_detail::ValidateIndexLayout(data.output_path);
	} catch (const std::exception &e) {
		throw IOException("%s: krepp reported success but '%s' does not hold a complete index: %s", kCallerName,
		                  data.output_path, std::string(e.what()));
	}
	// The partial THIS build wrote, not partials.begin() - that is the
	// lexicographically first suffix in the directory, which in a multi-partial
	// build is some earlier residue. Reporting its k/w/h beside status='ok'
	// would describe an index this call did not write.
	ReadIndexMetadata(data.output_path, PartialSuffixFor(options), gstate->k, gstate->w, gstate->h);
	gstate->num_references = num_references;

	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState>
KreppIndexCreateTableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void KreppIndexCreateTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.done) {
		output.SetCardinality(0);
		return;
	}

	output.data[0].SetValue(0, Value(bind_data.output_path));
	output.data[1].SetValue(0, Value::INTEGER(gstate.k));
	output.data[2].SetValue(0, Value::INTEGER(gstate.w));
	output.data[3].SetValue(0, Value::INTEGER(gstate.h));
	output.data[4].SetValue(0, Value::BIGINT(gstate.num_references));
	output.data[5].SetValue(0, Value("ok"));
	output.SetCardinality(1);
	gstate.done = true;
}

TableFunction KreppIndexCreateTableFunction::GetFunction() {
	auto tf =
	    TableFunction(kCallerName, {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	AddNamedParameter(tf, "tree_table", LogicalType::VARCHAR);
	AddNamedParameter(tf, "newick_path", LogicalType::VARCHAR);
	AddNamedParameter(tf, "k", LogicalType::INTEGER);
	AddNamedParameter(tf, "w", LogicalType::INTEGER);
	AddNamedParameter(tf, "h", LogicalType::INTEGER);
	AddNamedParameter(tf, "m", LogicalType::INTEGER);
	AddNamedParameter(tf, "r", LogicalType::INTEGER);
	AddNamedParameter(tf, "frac", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "threads", LogicalType::INTEGER);
	AddNamedParameter(tf, "sdust_t", LogicalType::INTEGER);
	AddNamedParameter(tf, "sdust_w", LogicalType::INTEGER);
	return tf;
}

void KreppIndexCreateTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
