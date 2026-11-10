#include "read_sra.hpp"

#include "miint_named_parameter.hpp"
#include "remote_file_helper.hpp"
#include "rype_input_stream.hpp"
#include "table_function_common.hpp"

#include "duckdb/common/error_data.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/list_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/execution/partition_info.hpp"
#include "duckdb/main/client_context.hpp"

#include <algorithm>
#include <cstring>
#include <mutex>

namespace duckdb {

using namespace miint::sra;

namespace {

// Consecutive spots decoded by one thread and emitted under one batch index. Large enough that
// an order-keeping sink holds few batches (docs/internals/duckdb-engine-notes.md, "Output order
// costs threads and memory"), small enough to spread a run over the threads. The setting
// miint_sra_spots_per_unit overrides it.
constexpr uint64_t DEFAULT_SPOTS_PER_UNIT = 128 * 1024;

// The most threads one scan decodes on (GlobalState::MaxThreads).
constexpr idx_t MAX_SCAN_THREADS = 16;

// Random access to a .sra through DuckDB's file system, so local paths and (with httpfs) remote
// URLs are read the same way. FileHandle::Read at a position is safe from several threads when
// the handle is opened for parallel access.
class DuckDBSraSource : public SraByteSource {
public:
	DuckDBSraSource(FileSystem &fs, const std::string &path, bool direct_io)
	    : handle_(fs.OpenFile(path, FileFlags::FILE_FLAGS_READ | FileFlags::FILE_FLAGS_PARALLEL_ACCESS |
	                                    (direct_io ? FileFlags::FILE_FLAGS_DIRECT_IO : FileOpenFlags()))),
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
[[noreturn]] void ThrowForRun(const SraScanData &bind, const std::string &accession, const std::exception &e) {
	std::string msg = bind.function + ": " + accession + ": " + e.what();
	if (dynamic_cast<const SraUnsupportedError *>(&e)) {
		throw InvalidInputException(msg);
	}
	throw IOException(msg);
}

// A run that could not be fetched: the I/O and HTTP errors DuckDB's file systems raise (a 404 is
// HTTP, a refused connection IO). Anything else (an interrupt, a deadline, an exhausted resource,
// a refused permission) is the query's failure, never a run's, and is not reworded as one.
bool IsFetchFailure(const std::exception &e) {
	auto type = ErrorData(e).Type();
	return type == ExceptionType::IO || type == ExceptionType::HTTP;
}

struct OpenRun {
	std::string location;
	std::string accession;
	std::shared_ptr<const SraTable> table;
	// Hints each claimed unit's blobs to the byte source, under the global lock. Execute also
	// passes it, outside the lock, as share_from to the readers of this run: that reads only what
	// its constructor wrote (the table and the name index), which PrefetchRows does not touch.
	mutable std::unique_ptr<SpotReader> prefetcher;
};

struct WorkUnit {
	std::shared_ptr<const OpenRun> run;
	int64_t first_row = 0;
	int64_t end_row = 0;
	idx_t batch_index = 0;
};

// Opens a run: a local file, or a remote one read by direct range requests.
std::shared_ptr<const SraTable> OpenTable(ClientContext &context, const std::string &location) {
	auto &fs = FileSystem::GetFileSystem(context);
	bool remote = miint::RemoteFileHelper::IsRemotePath(location);
	auto source = std::make_shared<DuckDBSraSource>(fs, location, remote);
	return std::make_shared<const SraTable>(std::make_shared<const KarArchive>(source));
}

struct GlobalState : public GlobalTableFunctionState {
	std::mutex lock;
	size_t next_run = 0;
	std::shared_ptr<const OpenRun> current;
	int64_t next_row = 0;
	int64_t end_row = 0;
	idx_t next_batch = 0;
	uint64_t spots_per_unit = DEFAULT_SPOTS_PER_UNIT;

	idx_t MaxThreads() const override {
		// A run is split into many units, so even one run uses several threads, up to DuckDB's
		// threads setting.
		return MAX_SCAN_THREADS;
	}

	// Hands out the next unit, in run and spot order, opening the next run when the current one
	// is used up. The order of claims is the order of batch indexes, hence of the output.
	bool Claim(ClientContext &context, const SraScanData &bind, WorkUnit &unit) {
		std::lock_guard<std::mutex> guard(lock);
		while (!current || next_row >= end_row) {
			current.reset();
			if (next_run >= bind.runs.size()) {
				return false;
			}
			const SraScanData::Run &input = bind.runs[next_run++];
			auto run = std::make_shared<OpenRun>();
			run->location = input.location;
			run->accession = input.accession;
			try {
				run->table = OpenTable(context, input.location);
				ClassifyRun(GatherRunFacts(*run->table));
				run->prefetcher = std::make_unique<SpotReader>(run->table);
				next_row = run->prefetcher->FirstRow();
				end_row = next_row + int64_t(run->prefetcher->RowCount());
			} catch (const SraFormatError &e) {
				ThrowForRun(bind, run->accession, e);
			} catch (const SraUnsupportedError &e) {
				ThrowForRun(bind, run->accession, e);
			} catch (const std::exception &e) {
				if (!IsFetchFailure(e)) {
					throw;
				}
				// The run could not be opened: the error names it.
				throw IOException(bind.function + ": " + run->accession + ": " + ErrorData(e).RawMessage());
			}
			current = std::move(run);
		}
		unit.run = current;
		unit.first_row = next_row;
		unit.end_row = std::min<int64_t>(end_row, next_row + int64_t(spots_per_unit));
		unit.batch_index = next_batch++;
		next_row = unit.end_row;
		// A thread's read-ahead follows its position, but each claim jumps it ahead of the other
		// threads: fetch the unit's blobs, and the next unit's, before a thread arrives at them.
		current->prefetcher->PrefetchRows(unit.first_row,
		                                  std::min<int64_t>(end_row, unit.end_row + int64_t(spots_per_unit)));
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
	// A NULL left in a mask buffer from an earlier chunk must not survive onto a valid row; a mask
	// without a buffer is already all-valid (as in SetResultVectorListUInt8).
	auto &validity = FlatVector::ValidityMutable(vec);
	if (validity.CanHaveNull()) {
		validity.SetAllValid(quals.size());
	}
	auto &child_validity = FlatVector::ValidityMutable(child);
	if (child_validity.CanHaveNull()) {
		child_validity.SetAllValid(total);
	}
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

// The options read_sra takes, and the schema they imply.
void BindOptions(SraScanData &data, named_argument_map_t &params, vector<LogicalType> &return_types,
                 vector<Identifier> &names) {
	const std::string &fn = data.function;
	auto split = params.find("split");
	if (split != params.end() && !split->second.IsNull()) {
		std::string s = split->second.ToString();
		if (s == "split_3") {
			data.split = SplitMode::Split3;
		} else if (s == "spot") {
			data.split = SplitMode::Spot;
		} else {
			throw InvalidInputException(fn + ": split must be 'split_3' or 'spot', not '" + s + "'");
		}
	}
	auto tech = params.find("include_technical");
	if (tech != params.end() && !tech->second.IsNull()) {
		data.include_technical = tech->second.GetValue<bool>();
	}
	if (data.include_technical && data.split == SplitMode::Split3) {
		// fasterq-dump --split-3 --include-technical switches to split-files with a warning.
		throw InvalidInputException(fn + ": include_technical requires split => 'spot'");
	}
	data.include_filepath = ParseIncludeFilepathParameter(params);

	data.names = {"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"};
	data.types = {LogicalType::BIGINT,
	              LogicalType::VARCHAR,
	              LogicalType::VARCHAR,
	              LogicalType::VARCHAR,
	              LogicalType::VARCHAR,
	              LogicalType::LIST(LogicalType::UTINYINT),
	              LogicalType::LIST(LogicalType::UTINYINT)};
	if (data.split == SplitMode::Spot) {
		data.names.emplace_back("read_number");
		data.types.emplace_back(LogicalType::UTINYINT);
	}
	if (data.include_filepath) {
		data.names.emplace_back("filepath");
		data.types.emplace_back(LogicalType::VARCHAR);
	}
	for (idx_t i = 0; i < data.names.size(); i++) {
		names.emplace_back(data.names[i]);
		return_types.push_back(data.types[i]);
	}
}

unique_ptr<FunctionData> ReadSraBind(ClientContext &context, TableFunctionBindInput &input,
                                     vector<LogicalType> &return_types, vector<Identifier> &names) {
	auto data = make_uniq<SraScanData>();
	data->function = "read_sra";
	auto &fs = FileSystem::GetFileSystem(context);
	const Value &arg = input.inputs[0];
	if (arg.IsNull()) {
		throw InvalidInputException("read_sra: the path cannot be NULL");
	}
	std::vector<std::string> paths;
	if (arg.type().id() == LogicalTypeId::VARCHAR) {
		paths = ExpandGlobPatternWithInfo(fs, context, arg.ToString()).paths;
	} else if (arg.type().id() == LogicalTypeId::LIST) {
		paths = ExpandPathList(fs, arg, "read_sra");
	} else {
		throw InvalidInputException("read_sra: the path must be VARCHAR or VARCHAR[]");
	}
	for (const auto &path : paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException("read_sra: stdin is not supported (an .sra file is read at random positions)");
		}
		if (!miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
		data->runs.push_back({path, ReadSraTableFunction::AccessionFromPath(path)});
	}
	BindOptions(*data, input.named_parameters, return_types, names);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto state = make_uniq<GlobalState>();
	idx_t unit = GetNonNegativeOption(context, "miint_sra_spots_per_unit", 0, "0 uses the default");
	if (unit > 0) {
		// Capped far above any run's spot count, so row arithmetic cannot overflow.
		state->spots_per_unit = std::min<idx_t>(unit, idx_t(1) << 40);
	}
	return std::move(state);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                              GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// Moves up to `capacity` rows of the thread's current unit into `rows`, decoding spots as needed.
// Clears has_unit when the unit is used up.
void FillRows(LocalState &ls, const SraScanData &bind, idx_t capacity, std::vector<SraRow> &rows) {
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
	auto &bind = data_p.bind_data->Cast<SraScanData>();
	auto &gs = data_p.global_state->Cast<GlobalState>();
	auto &ls = data_p.local_state->Cast<LocalState>();

	std::vector<SraRow> rows;
	rows.reserve(STANDARD_VECTOR_SIZE);
	while (rows.empty()) {
		if (!ls.has_unit && ls.pending_pos >= ls.pending.size()) {
			context.BlockingInterruptCheck();
			if (!gs.Claim(context, bind, ls.unit)) {
				output.SetChildCardinality(0);
				return;
			}
			ls.has_unit = true;
			ls.batch_index = ls.unit.batch_index;
			ls.next_row = ls.unit.first_row;
			if (ls.reader_run != ls.unit.run) {
				try {
					// The run's name index is shared with the reader Claim made, not read again.
					ls.reader = std::make_unique<SpotReader>(ls.unit.run->table, ls.unit.run->prefetcher.get());
				} catch (const SraFormatError &e) {
					ThrowForRun(bind, ls.unit.run->accession, e);
				} catch (const SraUnsupportedError &e) {
					ThrowForRun(bind, ls.unit.run->accession, e);
				}
				ls.reader_run = ls.unit.run;
			}
		}
		try {
			FillRows(ls, bind, STANDARD_VECTOR_SIZE, rows);
		} catch (const SraFormatError &e) {
			ThrowForRun(bind, ls.unit.run->accession, e);
		} catch (const SraUnsupportedError &e) {
			ThrowForRun(bind, ls.unit.run->accession, e);
		} catch (const std::exception &e) {
			// A run that fails partway (a fetch the network refused) fails the query: its rows so
			// far are already out. The error names the run.
			if (!IsFetchFailure(e)) {
				throw;
			}
			throw IOException(bind.function + ": " + ls.unit.run->accession + ": " + ErrorData(e).RawMessage());
		}
	}

	const idx_t n = rows.size();
	const std::string &acc = ls.unit.run->accession;
	auto seq_index = FlatVector::GetDataMutable<int64_t>(output.data[0]);
	std::vector<std::string> read_ids(n), comments(n), seq1(n), seq2(n);
	std::vector<const std::vector<uint8_t> *> qual1(n), qual2(n);
	for (idx_t i = 0; i < n; i++) {
		const SraRow &r = rows[i];
		seq_index[i] = r.spot;
		read_ids[i] = acc + "." + std::to_string(r.spot);
		// fasterq-dump's defline continues "$sn length=$rl"; read_fastx keeps read 1's.
		comments[i] = r.read1.has_name ? r.read1.name + " " : std::string();
		comments[i] += "length=" + std::to_string(r.read1.bases.size());
		seq1[i] = r.read1.bases;
		qual1[i] = &r.read1.quality;
		if (r.paired) {
			seq2[i] = r.read2.bases;
			qual2[i] = &r.read2.quality;
		}
	}
	idx_t col = 1;
	SetResultVectorString(output.data[col++], read_ids);
	SetResultVectorString(output.data[col++], comments);
	SetResultVectorString(output.data[col++], seq1);
	SetResultVectorStringNullable(output.data[col++], seq2);
	SetQualityList(output.data[col++], qual1);
	SetQualityList(output.data[col++], qual2);
	if (bind.split == SplitMode::Spot) {
		auto numbers = FlatVector::GetDataMutable<uint8_t>(output.data[col++]);
		for (idx_t i = 0; i < n; i++) {
			if (rows[i].read_number > 255) {
				throw InvalidInputException(bind.function + ": " + acc + ": spot " + std::to_string(rows[i].spot) +
				                            " has more than 255 reads");
			}
			numbers[i] = uint8_t(rows[i].read_number);
		}
	}
	if (bind.include_filepath) {
		SetResultVectorFilepath(output.data[col++], ls.unit.run->location);
	}
	output.SetChildCardinality(n);
}

OperatorPartitionData GetPartitionData(ClientContext &context, TableFunctionGetPartitionInput &input) {
	if (input.partition_info.RequiresPartitionColumns()) {
		throw InternalException("SRA scan: partition columns are not supported");
	}
	return OperatorPartitionData(input.local_state->Cast<LocalState>().batch_index);
}

// About 300 bytes of .sra per spot in the surveyed runs; a lower figure errs high, the safe side
// (EstimateFileScanCardinality). A remote run counts as large.
unique_ptr<NodeStatistics> Cardinality(ClientContext &context, const FunctionData *bind_data) {
	std::vector<std::string> locations;
	for (const auto &r : bind_data->Cast<SraScanData>().runs) {
		locations.push_back(r.location);
	}
	return EstimateFileScanCardinality(locations, 100);
}

TableFunction MakeScan(const char *name, table_function_bind_t bind) {
	TableFunction tf(name, {LogicalType::ANY}, Execute, bind, InitGlobal, InitLocal);
	AddNamedParameter(tf, "split", LogicalType::VARCHAR);
	AddNamedParameter(tf, "include_technical", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "include_filepath", LogicalType::BOOLEAN);
	tf.get_partition_data = GetPartitionData;
	tf.cardinality = Cardinality;
	return tf;
}

} // namespace

std::string ReadSraTableFunction::AccessionFromPath(const std::string &path) {
	auto slash = path.find_last_of("/\\");
	std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
	for (const std::string ext : {".sralite", ".sra"}) {
		if (name.size() > ext.size() && name.compare(name.size() - ext.size(), ext.size(), ext) == 0) {
			return name.substr(0, name.size() - ext.size());
		}
	}
	return name;
}

TableFunction ReadSraTableFunction::GetFunction() {
	return MakeScan("read_sra", ReadSraBind);
}

void ReadSraTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
