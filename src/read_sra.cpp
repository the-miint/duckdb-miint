#include "read_sra.hpp"

#include "miint_named_parameter.hpp"
#include "remote_file_helper.hpp"
#include "rype_input_stream.hpp"
#include "table_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/list_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/execution/partition_info.hpp"
#include "duckdb/main/client_context.hpp"

#include <mutex>

namespace duckdb {

using namespace miint::sra;

namespace {

const char *FUNCTION = "read_sra";

// Consecutive spots decoded by one thread and emitted under one batch index. Large enough that
// an order-keeping sink holds few batches (docs/internals/duckdb-engine-notes.md, "Output order
// costs threads and memory"), small enough to spread a run over the threads. The setting
// miint_sra_spots_per_unit overrides it.
constexpr uint64_t DEFAULT_SPOTS_PER_UNIT = 128 * 1024;

// Random access to a .sra through DuckDB's file system, so local paths and (with httpfs) remote
// URLs are read the same way. FileHandle::Read at a position is safe from several threads when
// the handle is opened for parallel access.
class DuckDBSraSource : public SraByteSource {
public:
	DuckDBSraSource(FileSystem &fs, const std::string &path)
	    : handle_(fs.OpenFile(path, FileFlags::FILE_FLAGS_READ | FileFlags::FILE_FLAGS_PARALLEL_ACCESS)),
	      size_(handle_->GetFileSize()) {
	}
	uint64_t Size() const override {
		return size_;
	}
	void ReadAt(uint64_t pos, uint8_t *buf, size_t n) const override {
		if (pos > size_ || n > size_ - pos) {
			throw SraFormatError("read past the end of the file");
		}
		handle_->Read(buf, n, pos);
	}

private:
	unique_ptr<FileHandle> handle_;
	uint64_t size_;
};

// A run's errors name the function and the run: a refusal (what the run is) is invalid input,
// damage is an I/O error.
[[noreturn]] void ThrowForRun(const std::string &accession, const std::exception &e) {
	std::string msg = std::string(FUNCTION) + ": " + accession + ": " + e.what();
	if (dynamic_cast<const SraUnsupportedError *>(&e)) {
		throw InvalidInputException(msg);
	}
	throw IOException(msg);
}

struct OpenRun {
	std::string path;
	std::string accession;
	std::shared_ptr<const SraTable> table;
};

struct WorkUnit {
	std::shared_ptr<const OpenRun> run;
	int64_t first_row = 0;
	int64_t end_row = 0;
	idx_t batch_index = 0;
};

struct GlobalState : public GlobalTableFunctionState {
	std::mutex lock;
	std::vector<std::string> paths;
	size_t next_path = 0;
	std::shared_ptr<const OpenRun> current;
	int64_t next_row = 0;
	int64_t end_row = 0;
	idx_t next_batch = 0;
	uint64_t spots_per_unit = DEFAULT_SPOTS_PER_UNIT;

	idx_t MaxThreads() const override {
		// A run is split into many units, so even one file uses several threads; DuckDB's
		// threads setting bounds this.
		return 16;
	}

	// Hands out the next unit, in run and spot order, opening the next run when the current one
	// is used up. The order of claims is the order of batch indexes, hence of the output.
	bool Claim(ClientContext &context, WorkUnit &unit) {
		std::lock_guard<std::mutex> guard(lock);
		while (!current || next_row >= end_row) {
			if (next_path >= paths.size()) {
				return false;
			}
			const std::string &path = paths[next_path++];
			auto run = std::make_shared<OpenRun>();
			run->path = path;
			run->accession = ReadSraTableFunction::AccessionFromPath(path);
			try {
				auto source = std::make_shared<DuckDBSraSource>(FileSystem::GetFileSystem(context), path);
				run->table = std::make_shared<const SraTable>(std::make_shared<const KarArchive>(source));
				ClassifyRun(GatherRunFacts(*run->table));
				SpotReader probe(run->table);
				next_row = probe.FirstRow();
				end_row = next_row + int64_t(probe.RowCount());
			} catch (const SraFormatError &e) {
				ThrowForRun(run->accession, e);
			} catch (const SraUnsupportedError &e) {
				ThrowForRun(run->accession, e);
			}
			current = std::move(run);
		}
		unit.run = current;
		unit.first_row = next_row;
		unit.end_row = std::min<int64_t>(end_row, next_row + int64_t(spots_per_unit));
		unit.batch_index = next_batch++;
		next_row = unit.end_row;
		return true;
	}
};

struct LocalState : public LocalTableFunctionState {
	bool has_unit = false;
	WorkUnit unit;
	idx_t batch_index = 0; // of the chunk last returned
	int64_t next_row = 0;
	std::shared_ptr<const OpenRun> reader_run;
	std::unique_ptr<SpotReader> reader;
	Spot spot;
	std::vector<SraRead> reads;
	std::vector<SraRow> pending;
	size_t pending_pos = 0;
};

void SetQualityList(Vector &vec, const std::vector<const std::vector<uint8_t> *> &quals) {
	idx_t total = 0;
	for (auto *q : quals) {
		total += q ? q->size() : 0;
	}
	ListVector::Reserve(vec, total);
	ListVector::SetListSize(vec, total);
	auto &child = ListVector::GetChildMutable(vec);
	auto child_data = FlatVector::GetDataMutable<uint8_t>(child);
	auto entries = FlatVector::GetDataMutable<list_entry_t>(vec);
	auto &validity = FlatVector::ValidityMutable(vec);
	idx_t offset = 0;
	for (idx_t i = 0; i < quals.size(); i++) {
		const auto *q = quals[i];
		entries[i].offset = offset;
		entries[i].length = q ? q->size() : 0;
		if (!q) {
			validity.SetInvalid(i);
			continue;
		}
		if (!q->empty()) {
			std::memcpy(child_data + offset, q->data(), q->size());
		}
		offset += q->size();
	}
}

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	auto data = make_uniq<ReadSraTableFunction::Data>();
	auto &fs = FileSystem::GetFileSystem(context);
	const Value &arg = input.inputs[0];
	if (arg.IsNull()) {
		throw InvalidInputException("read_sra: the path cannot be NULL");
	}
	if (arg.type().id() == LogicalTypeId::VARCHAR) {
		data->paths = ExpandGlobPatternWithInfo(fs, context, arg.ToString()).paths;
	} else if (arg.type().id() == LogicalTypeId::LIST) {
		data->paths = ExpandPathList(fs, arg, FUNCTION);
	} else {
		throw InvalidInputException("read_sra: the path must be VARCHAR or VARCHAR[]");
	}
	for (const auto &path : data->paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException("read_sra: stdin is not supported (an .sra file is read at random positions)");
		}
		if (!miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	auto &params = input.named_parameters;
	auto split = params.find("split");
	if (split != params.end() && !split->second.IsNull()) {
		std::string s = split->second.ToString();
		if (s == "split_3") {
			data->split = SplitMode::Split3;
		} else if (s == "spot") {
			data->split = SplitMode::Spot;
		} else {
			throw InvalidInputException("read_sra: split must be 'split_3' or 'spot', not '" + s + "'");
		}
	}
	auto tech = params.find("include_technical");
	if (tech != params.end() && !tech->second.IsNull()) {
		data->include_technical = tech->second.GetValue<bool>();
	}
	if (data->include_technical && data->split == SplitMode::Split3) {
		// fasterq-dump --split-3 --include-technical switches to split-files with a warning.
		throw InvalidInputException("read_sra: include_technical requires split => 'spot'");
	}
	data->include_filepath = ParseIncludeFilepathParameter(params);

	data->names = {"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"};
	data->types = {LogicalType::BIGINT,
	               LogicalType::VARCHAR,
	               LogicalType::VARCHAR,
	               LogicalType::VARCHAR,
	               LogicalType::VARCHAR,
	               LogicalType::LIST(LogicalType::UTINYINT),
	               LogicalType::LIST(LogicalType::UTINYINT)};
	if (data->split == SplitMode::Spot) {
		data->names.emplace_back("read_number");
		data->types.emplace_back(LogicalType::UTINYINT);
	}
	if (data->include_filepath) {
		data->names.emplace_back("filepath");
		data->types.emplace_back(LogicalType::VARCHAR);
	}
	for (idx_t i = 0; i < data->names.size(); i++) {
		names.emplace_back(data->names[i]);
		return_types.push_back(data->types[i]);
	}
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto state = make_uniq<GlobalState>();
	state->paths = input.bind_data->Cast<ReadSraTableFunction::Data>().paths;
	idx_t unit = GetNonNegativeOption(context, "miint_sra_spots_per_unit", 0, "0 uses the default");
	if (unit > 0) {
		state->spots_per_unit = unit;
	}
	return std::move(state);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                              GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// Moves up to `capacity` rows of the thread's current unit into `rows`, decoding spots as needed.
// Clears has_unit when the unit is used up.
void FillRows(LocalState &ls, const ReadSraTableFunction::Data &bind, idx_t capacity, std::vector<SraRow> &rows) {
	while (rows.size() < capacity) {
		if (ls.pending_pos < ls.pending.size()) {
			rows.push_back(std::move(ls.pending[ls.pending_pos++]));
			continue;
		}
		if (ls.next_row >= ls.unit.end_row) {
			ls.has_unit = false;
			return;
		}
		ls.pending.clear();
		ls.pending_pos = 0;
		ls.reader->Read(ls.next_row++, ls.spot);
		ls.reads.clear();
		RouteSpot(ls.spot, bind.split, bind.include_technical, ls.reads);
		RowsForSpot(ls.reads, bind.split, ls.pending);
	}
}

void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind = data_p.bind_data->Cast<ReadSraTableFunction::Data>();
	auto &gs = data_p.global_state->Cast<GlobalState>();
	auto &ls = data_p.local_state->Cast<LocalState>();

	std::vector<SraRow> rows;
	rows.reserve(STANDARD_VECTOR_SIZE);
	while (rows.empty()) {
		if (!ls.has_unit && ls.pending_pos >= ls.pending.size()) {
			context.BlockingInterruptCheck();
			if (!gs.Claim(context, ls.unit)) {
				output.SetChildCardinality(0);
				return;
			}
			ls.has_unit = true;
			ls.batch_index = ls.unit.batch_index;
			ls.next_row = ls.unit.first_row;
			if (ls.reader_run != ls.unit.run) {
				try {
					ls.reader = std::make_unique<SpotReader>(ls.unit.run->table);
				} catch (const SraFormatError &e) {
					ThrowForRun(ls.unit.run->accession, e);
				} catch (const SraUnsupportedError &e) {
					ThrowForRun(ls.unit.run->accession, e);
				}
				ls.reader_run = ls.unit.run;
			}
		}
		try {
			FillRows(ls, bind, STANDARD_VECTOR_SIZE, rows);
		} catch (const SraFormatError &e) {
			ThrowForRun(ls.unit.run->accession, e);
		} catch (const SraUnsupportedError &e) {
			ThrowForRun(ls.unit.run->accession, e);
		}
	}

	const idx_t n = rows.size();
	const std::string &acc = ls.unit.run->accession;
	auto seq_index = FlatVector::GetDataMutable<int64_t>(output.data[0]);
	std::vector<std::string> read_ids(n), seq1(n), seq2(n);
	std::vector<const std::vector<uint8_t> *> qual1(n), qual2(n);
	for (idx_t i = 0; i < n; i++) {
		const SraRow &r = rows[i];
		seq_index[i] = r.spot;
		read_ids[i] = acc + "." + std::to_string(r.spot);
		seq1[i] = r.read1.bases;
		qual1[i] = &r.read1.quality;
		if (r.paired) {
			seq2[i] = r.read2.bases;
			qual2[i] = &r.read2.quality;
		}
	}
	idx_t col = 1;
	SetResultVectorString(output.data[col++], read_ids);
	SetResultVectorNull(output.data[col++]); // comment: the spot name, reconstructed in a later phase
	SetResultVectorString(output.data[col++], seq1);
	SetResultVectorStringNullable(output.data[col++], seq2);
	SetQualityList(output.data[col++], qual1);
	SetQualityList(output.data[col++], qual2);
	if (bind.split == SplitMode::Spot) {
		auto numbers = FlatVector::GetDataMutable<uint8_t>(output.data[col++]);
		for (idx_t i = 0; i < n; i++) {
			if (rows[i].read_number > 255) {
				throw InvalidInputException("read_sra: " + acc + ": spot " + std::to_string(rows[i].spot) +
				                            " has more than 255 reads");
			}
			numbers[i] = uint8_t(rows[i].read_number);
		}
	}
	if (bind.include_filepath) {
		SetResultVectorFilepath(output.data[col++], ls.unit.run->path);
	}
	output.SetChildCardinality(n);
}

OperatorPartitionData GetPartitionData(ClientContext &context, TableFunctionGetPartitionInput &input) {
	if (input.partition_info.RequiresPartitionColumns()) {
		throw InternalException("read_sra: partition columns are not supported");
	}
	return OperatorPartitionData(input.local_state->Cast<LocalState>().batch_index);
}

// About 300 bytes of .sra per spot in the surveyed runs; a lower figure errs high, the safe side
// (EstimateFileScanCardinality).
unique_ptr<NodeStatistics> Cardinality(ClientContext &context, const FunctionData *bind_data) {
	return EstimateFileScanCardinality(bind_data->Cast<ReadSraTableFunction::Data>().paths, 100);
}

} // namespace

std::string ReadSraTableFunction::AccessionFromPath(const std::string &path) {
	std::string name = path.substr(path.find_last_of("/\\") == std::string::npos ? 0 : path.find_last_of("/\\") + 1);
	for (const char *ext : {".sralite", ".sra"}) {
		std::string e = ext;
		if (name.size() > e.size() && name.compare(name.size() - e.size(), e.size(), e) == 0) {
			return name.substr(0, name.size() - e.size());
		}
	}
	return name;
}

TableFunction ReadSraTableFunction::GetFunction() {
	TableFunction tf(FUNCTION, {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	AddNamedParameter(tf, "split", LogicalType::VARCHAR);
	AddNamedParameter(tf, "include_technical", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "include_filepath", LogicalType::BOOLEAN);
	tf.get_partition_data = GetPartitionData;
	tf.cardinality = Cardinality;
	return tf;
}

void ReadSraTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
