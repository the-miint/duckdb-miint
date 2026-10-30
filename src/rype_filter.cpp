#include "rype_filter.hpp"
#include "rype_filter_schedule.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"
#include "miint_named_parameter.hpp"
#include "rype_common.hpp"

#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/function_set.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <chrono>
#include <exception>
#include <condition_variable>
#include <deque>
#include <mutex>

namespace duckdb {

namespace {

//! A mean sequence1 length at or above this, sampled from the head of the scan,
//! is long reads, which minimap2's short-read preset is not for: rype_filter then
//! runs RYpe alone and never opens the minimap2 index.
constexpr size_t LONG_READ_MEAN_LENGTH = 300;

//! minimap2 as align_minimap2 runs it by default (the sr preset, five secondary
//! alignments), less what only shapes SAM rows. =/X CIGAR operations are written
//! after an alignment is scored, so they cannot change whether a read maps.
miint::Minimap2Config FilterMinimap2Config() {
	miint::Minimap2Config config;
	config.eqx = false;
	return config;
}

//! Read the session-scoped `miint_rype_filter_group_rows` setting: rows per
//! group, or 0 to size groups from the memory budget.
idx_t GetRypeFilterGroupRows(ClientContext &context) {
	return GetNonNegativeOption(context, "miint_rype_filter_group_rows", 0, "0 sizes groups from the memory budget");
}

//! The most rows one group should carry within `budget` bytes; see
//! Schedule::GroupRows.
idx_t AutoGroupRows(ClientContext &context, const RypeIndex *index, const RypeInputStream::Sizing &sizing,
                    size_t budget) {
	const auto config = CalculateRypeBatchConfig(context, "rype_filter", index, sizing.avg_read_length,
	                                             sizing.is_paired ? 1 : 0, budget);
	if (!config) {
		return STANDARD_VECTOR_SIZE;
	}
	return RypeFilterTableFunction::Schedule::GroupRows(config->per_batch_memory, config->batch_size,
	                                                    sizing.held_bytes_per_row);
}

//! A thread that cannot park waits this long at most before looking again, so
//! that it notices an interrupted query.
constexpr auto WAIT_FOR_PROGRESS = std::chrono::milliseconds(10);

} // namespace

// ============================================================================
// Schedule
// ============================================================================

void RypeFilterTableFunction::Schedule::Release() {
	// Every RYpe output stream was released inside the call that made it, so
	// nothing of RYpe's is left reading the input stream or the index.
	units.clear();
	aligned.clear();
	minimap2_index.reset();
	if (index) {
		rype_index_free(index);
		index = nullptr;
	}
	// The stream holds the streaming result, which points into the connection's
	// context: stream first, connection last.
	input_stream.reset();
	input_connection.reset();
}

RypeFilterTableFunction::Schedule::~Schedule() {
	Release();
}

idx_t RypeFilterTableFunction::Schedule::GroupRows(idx_t per_batch_memory, idx_t batch_size, idx_t held_bytes_per_row) {
	// RYpe's estimate is linear in its batch: c = per_batch_memory / batch_size a
	// read. A group of G rows costs G * c while RYpe classifies it, and its rows
	// are held, at h bytes each (every column, not just the sequences), until
	// emitted; so are the rows of the group before it, and of the one after,
	// read ahead. G * (c + 3h) must fit the per_batch_memory RYpe would have
	// spent alone.
	const double per_read = static_cast<double>(per_batch_memory) / static_cast<double>(batch_size);
	const double groups_in_memory = static_cast<double>(MAX_HELD_GROUPS + 1);
	const double rows =
	    static_cast<double>(per_batch_memory) / (per_read + groups_in_memory * static_cast<double>(held_bytes_per_row));
	return MaxValue<idx_t>(1, MinValue<idx_t>(static_cast<idx_t>(rows), batch_size));
}

bool RypeFilterTableFunction::Schedule::TakeClassifyLocked() {
	if (error.HasError() || classifying || input_done || held_groups >= MAX_HELD_GROUPS) {
		return false;
	}
	classifying = true;
	held_groups++;
	return true;
}

bool RypeFilterTableFunction::Schedule::TakeIndexLoadLocked() {
	if (error.HasError() || !use_minimap2 || minimap2_index || index_loading) {
		return false;
	}
	index_loading = true;
	return true;
}

idx_t RypeFilterTableFunction::Schedule::TakeAlignLocked() {
	if (error.HasError() || !minimap2_index || aligners >= max_aligners) {
		return 0;
	}
	const idx_t start = MinValue<idx_t>(max_aligners - aligners, units.size());
	aligners += start;
	return start;
}

RypeFilterJobs RypeFilterTableFunction::Schedule::TakeJobsLocked(bool inline_jobs) {
	RypeFilterJobs jobs;
	jobs.load_index = TakeIndexLoadLocked();
	if (!inline_jobs || !jobs.load_index) {
		jobs.classify = TakeClassifyLocked();
	}
	return jobs;
}

vector<InterruptState> RypeFilterTableFunction::Schedule::TakeWaitersLocked() {
	wakes++;
	vector<InterruptState> woken;
	woken.swap(parked);
	return woken;
}

RypeFilterFollowUp RypeFilterTableFunction::Schedule::PublishLocked(shared_ptr<RypeFilterGroup> group) {
	RypeFilterFollowUp follow_up;
	classifying = false;
	if (group) {
		group->runs_left = group->runs.size();
		for (idx_t run = 0; run < group->runs.size(); run++) {
			units.push_back(RypeFilterUnit {group, run});
		}
	} else {
		// Nothing more will be published: every parked thread can finish.
		input_done = true;
		held_groups--;
	}
	follow_up.wake = TakeWaitersLocked();
	// Carry on with the next group while there is room for it.
	follow_up.jobs.classify = TakeClassifyLocked();
	follow_up.jobs.align = TakeAlignLocked();
	return follow_up;
}

RypeFilterFollowUp
RypeFilterTableFunction::Schedule::IndexLoadedLocked(std::shared_ptr<miint::SharedMinimap2Index> loaded) {
	RypeFilterFollowUp follow_up;
	minimap2_index = std::move(loaded);
	index_loading = false;
	follow_up.wake = TakeWaitersLocked();
	follow_up.jobs.align = TakeAlignLocked();
	return follow_up;
}

RypeFilterFollowUp RypeFilterTableFunction::Schedule::FailLocked(ErrorData failure) {
	RypeFilterFollowUp follow_up;
	if (!error.HasError()) {
		error = std::move(failure);
	}
	follow_up.wake = TakeWaitersLocked();
	return follow_up;
}

bool RypeFilterTableFunction::Schedule::ClaimLocked(RypeFilterUnit &unit) {
	if (!aligned.empty()) {
		unit = std::move(aligned.front());
		aligned.pop_front();
		return true;
	}
	// Published units are the align jobs' when there are any.
	if (max_aligners > 0 || units.empty() || (use_minimap2 && !minimap2_index)) {
		return false;
	}
	unit = std::move(units.front());
	units.pop_front();
	return true;
}

bool RypeFilterTableFunction::Schedule::DrainedLocked() const {
	return input_done && units.empty() && aligned.empty() && units_aligning == 0;
}

bool RypeFilterTableFunction::Schedule::ClaimToAlignLocked(RypeFilterUnit &unit) {
	D_ASSERT(aligners > 0);
	if (error.HasError() || units.empty()) {
		aligners--;
		return false;
	}
	unit = std::move(units.front());
	units.pop_front();
	units_aligning++;
	return true;
}

RypeFilterFollowUp RypeFilterTableFunction::Schedule::AlignedLocked(RypeFilterUnit unit) {
	RypeFilterFollowUp follow_up;
	units_aligning--;
	unit.aligned = true;
	aligned.push_back(std::move(unit));
	follow_up.wake = TakeWaitersLocked();
	return follow_up;
}

RypeFilterFollowUp RypeFilterTableFunction::Schedule::FinishUnitLocked(const RypeFilterUnit &unit, bool inline_jobs) {
	RypeFilterFollowUp follow_up;
	if (--unit.group->runs_left > 0) {
		return follow_up;
	}
	held_groups--;
	if (inline_jobs) {
		// Room for the next group, which a waiting thread can now classify.
		follow_up.wake = TakeWaitersLocked();
	} else {
		// Room for the next group: classify it on the async pool.
		follow_up.jobs.classify = TakeClassifyLocked();
	}
	return follow_up;
}

namespace {

using Schedule = RypeFilterTableFunction::Schedule;

void Wake(Schedule &s, const vector<InterruptState> &parked) {
	s.progress.notify_all();
	// Every parked thread is called back even when one callback throws: one left
	// out would stay parked for good.
	std::exception_ptr first_error;
	for (auto &task : parked) {
		try {
			task.Callback();
		} catch (...) {
			if (!first_error) {
				first_error = std::current_exception();
			}
		}
	}
	if (first_error) {
		std::rethrow_exception(first_error);
	}
}

//! Record the first error and wake everyone, who will throw it.
void Fail(Schedule &s, ErrorData error) {
	RypeFilterFollowUp follow_up;
	{
		std::lock_guard<std::mutex> guard(s.lock);
		follow_up = s.FailLocked(std::move(error));
	}
	Wake(s, follow_up.wake);
}

void StartJobs(TaskExecutor *executor, const shared_ptr<Schedule> &schedule, const RypeFilterJobs &jobs);

//! Classify the next group with one RYpe call and mark the rows it matched.
//! Null when the scan has no rows left.
shared_ptr<RypeFilterGroup> ClassifyGroup(Schedule &s) {
	auto &input = *s.input_stream;
	if (!input.BeginGroup()) {
		return nullptr;
	}
	auto group = make_shared_ptr<RypeFilterGroup>();
	group->first_id = s.next_first_id;

	// The input stream stays owned by the schedule whatever RYpe does with it;
	// see the ownership note in rype_input_stream.hpp. The stream ends after one
	// group, so RYpe sees exactly the rows that group holds.
	ArrowArrayStreamWrapper output;
	if (rype_classify_arrow_ex(s.index, nullptr, &input.stream, s.threshold, s.group_rows,
	                           &output.arrow_array_stream) != 0) {
		throw IOException("RYpe classification failed: %s", RypeErrorText(rype_get_last_error()));
	}
	ArrowSchemaWrapper schema;
	if (output.arrow_array_stream.get_schema(&output.arrow_array_stream, &schema.arrow_schema) != 0) {
		throw IOException("Failed to get RYpe output schema: %s",
		                  RypeErrorText(output.arrow_array_stream.get_last_error(&output.arrow_array_stream)));
	}
	if (schema.arrow_schema.n_children < 1 || std::string(schema.arrow_schema.children[0]->format) != "l") {
		throw IOException("RYpe classify returned an unexpected schema: its first column is not an Int64 query_id");
	}

	// Only hits come back, a read once per bucket it reaches the threshold for.
	group->matched.assign(s.group_rows, false);
	idx_t matched_end = 0;
	for (;;) {
		ArrowArrayWrapper batch;
		if (output.arrow_array_stream.get_next(&output.arrow_array_stream, &batch.arrow_array) != 0) {
			throw IOException("Error getting next batch from RYpe: %s",
			                  RypeErrorText(output.arrow_array_stream.get_last_error(&output.arrow_array_stream)));
		}
		if (!batch.arrow_array.release) {
			break;
		}
		const auto &query_ids = *batch.arrow_array.children[0];
		if (batch.arrow_array.length > 0 && !query_ids.buffers[1]) {
			throw IOException("Arrow array has null data buffer in RYpe classify output");
		}
		const auto ids = static_cast<const int64_t *>(query_ids.buffers[1]);
		for (int64_t r = 0; r < batch.arrow_array.length; r++) {
			const int64_t id = ids[r + batch.arrow_array.offset + query_ids.offset];
			// The ids are the surrogates this group was sent, so anything else means
			// RYpe echoed something it was not given.
			if (id < 0 || static_cast<idx_t>(id) < group->first_id ||
			    static_cast<idx_t>(id) - group->first_id >= s.group_rows) {
				throw IOException("RYpe returned query_id %lld, outside the group it was sent (ids %llu to %llu)",
				                  static_cast<long long>(id), static_cast<unsigned long long>(group->first_id),
				                  static_cast<unsigned long long>(group->first_id + s.group_rows - 1));
			}
			const idx_t row = static_cast<idx_t>(id) - group->first_id;
			group->matched[row] = true;
			matched_end = MaxValue<idx_t>(matched_end, row + 1);
		}
	}

	// RYpe has pulled the group to its end-of-stream, so every row it was sent is
	// held now.
	group->runs = input.TakeGroup();
	idx_t group_size = 0;
	for (const auto &run : group->runs) {
		// A unit is one run and becomes at most one output chunk.
		D_ASSERT(run.to - run.from <= STANDARD_VECTOR_SIZE);
		group_size += run.to - run.from;
	}
	if (group_size == 0) {
		// The same rows would begin the next group, and this would never end.
		throw InternalException("rype_filter: RYpe finished a group without reading any of it");
	}
	D_ASSERT(group->runs.front().first_id == group->first_id);
	if (matched_end > group_size) {
		throw IOException("RYpe returned a query_id past the %llu rows of its group",
		                  static_cast<unsigned long long>(group_size));
	}
	s.next_first_id += group_size;
	return group;
}

//! The classify job: classify groups and publish their runs until MAX_HELD_GROUPS
//! are held or the scan ends. Whoever took the job runs it, so it runs alone.
//! `executor` runs the align jobs a publish starts; null when there are none.
void RunClassify(const shared_ptr<Schedule> &schedule, TaskExecutor *executor) {
	auto &s = *schedule;
	for (;;) {
		shared_ptr<RypeFilterGroup> group;
		try {
			group = ClassifyGroup(s);
		} catch (std::exception &e) {
			Fail(s, ErrorData(e));
			return;
		}
		RypeFilterFollowUp follow_up;
		{
			std::lock_guard<std::mutex> guard(s.lock);
			follow_up = s.PublishLocked(std::move(group));
		}
		Wake(s, follow_up.wake);
		// This job goes on to the next group itself.
		auto started = follow_up.jobs;
		started.classify = false;
		if (started.Any()) {
			StartJobs(executor, schedule, started);
		}
		if (!follow_up.jobs.classify) {
			return;
		}
	}
}

//! The index-load job. minimap2's loader cannot be interrupted, so a query that
//! ends meanwhile waits for it in teardown. A multi-part index is refused once its
//! first part has loaded, which is still before any row goes out: no unit is
//! claimed before the index is in. `executor` runs the align jobs it starts.
void RunIndexLoad(const shared_ptr<Schedule> &schedule, TaskExecutor *executor) {
	auto &s = *schedule;
	std::shared_ptr<miint::SharedMinimap2Index> loaded;
	bool multi_part = false;
	try {
		miint::Minimap2IndexReader reader(s.minimap2_path, FilterMinimap2Config(), s.index_load_threads);
		loaded = reader.ReadNextPart();
		if (!loaded) {
			throw std::runtime_error("it holds no index");
		}
		multi_part = !reader.AtEof();
	} catch (std::exception &e) {
		Fail(s, ErrorData(IOException("Failed to load minimap2 index '%s': %s", s.minimap2_path, e.what())));
		return;
	}
	// Each part would be loaded and mapped against in turn, and the reads held
	// meanwhile: the whole relation, or a second read of it.
	if (multi_part) {
		loaded.reset();
		Fail(s, ErrorData(InvalidInputException(
		            "rype_filter: '%s' is a multi-part minimap2 index, and rype_filter maps against a whole index at "
		            "once. Rebuild it as one part: save_minimap2_index always builds one, and minimap2 -d does when -I "
		            "is larger than the reference.",
		            s.minimap2_path)));
		return;
	}
	RypeFilterFollowUp follow_up;
	{
		std::lock_guard<std::mutex> guard(s.lock);
		follow_up = s.IndexLoadedLocked(std::move(loaded));
	}
	Wake(s, follow_up.wake);
	if (follow_up.jobs.Any()) {
		StartJobs(executor, schedule, follow_up.jobs);
	}
}

//! The rows of `unit` that neither stage removes, into `sel` as positions in its
//! chunk: RYpe did not match them and, given an aligner, minimap2 maps neither
//! mate. Returns how many.
idx_t SelectSurvivors(const Schedule &s, miint::Minimap2Aligner *aligner, const RypeFilterUnit &unit,
                      SelectionVector &sel) {
	const auto &group = *unit.group;
	const auto &run = group.runs[unit.run];
	sel.Initialize(run.to - run.from);
	idx_t count = 0;
	for (idx_t i = run.from; i < run.to; i++) {
		if (!group.matched[run.first_id + (i - run.from) - group.first_id]) {
			sel.set_index(count++, i);
		}
	}
	if (!aligner || count == 0) {
		return count;
	}
	// minimap2 reads the sequences where the chunk holds them, and names each read
	// by its read_id as align_minimap2 does; no name when there is none. Other
	// threads read the chunk meanwhile (the classify job, when two groups share
	// it), which is safe only because the stream's chunks are flat
	// (StreamingQuery): ToUnifiedFormat then changes nothing.
	auto &chunk = *run.chunk;
	UnifiedVectorFormat sequence1;
	UnifiedVectorFormat sequence2;
	D_ASSERT(chunk.data[s.sequence1_column].GetVectorType() == VectorType::FLAT_VECTOR);
	chunk.data[s.sequence1_column].ToUnifiedFormat(sequence1);
	if (s.sequence2_column != DConstants::INVALID_INDEX) {
		D_ASSERT(chunk.data[s.sequence2_column].GetVectorType() == VectorType::FLAT_VECTOR);
		chunk.data[s.sequence2_column].ToUnifiedFormat(sequence2);
	}
	// The names of the rows RYpe kept, in sel's order.
	std::vector<std::string> names;
	std::vector<bool> name_is_null;
	if (s.read_id_column != DConstants::INVALID_INDEX) {
		D_ASSERT(chunk.data[s.read_id_column].GetVectorType() == VectorType::FLAT_VECTOR);
		DataChunk kept_ids;
		kept_ids.InitializeEmpty({s.read_id_type});
		kept_ids.data[0].Slice(chunk.data[s.read_id_column], sel, count);
		kept_ids.SetCardinality(count);
		ExtractIdColumnAsStrings(kept_ids, 0, s.read_id_type, names, name_is_null);
	}
	const auto maps = [&](idx_t j) {
		const auto i = sel.get_index(j);
		const auto s1_idx = sequence1.sel->get_index(i);
		// The stream refused a NULL sequence1 before RYpe saw the row.
		D_ASSERT(sequence1.validity.RowIsValid(s1_idx));
		const auto s1 = UnifiedVectorFormat::GetData<string_t>(sequence1)[s1_idx];
		std::string_view mate;
		if (s.sequence2_column != DConstants::INVALID_INDEX) {
			const auto s2_idx = sequence2.sel->get_index(i);
			if (sequence2.validity.RowIsValid(s2_idx)) {
				const auto s2 = UnifiedVectorFormat::GetData<string_t>(sequence2)[s2_idx];
				mate = std::string_view(s2.GetData(), s2.GetSize());
			}
		}
		const char *name = names.empty() || name_is_null[j] ? nullptr : names[j].c_str();
		return aligner->any_mapped(std::string_view(s1.GetData(), s1.GetSize()), mate, name);
	};
	// Only what RYpe kept is mapped, and what minimap2 keeps stays in place.
	idx_t kept = 0;
	for (idx_t j = 0; j < count; j++) {
		if (!maps(j)) {
			sel.set_index(kept++, sel.get_index(j));
		}
	}
	return kept;
}

//! An align job: map published units until none is left, handing each to
//! DuckDB's threads to emit. Up to max_aligners run at once, each with its own
//! aligner on the one index.
void RunAlign(Schedule &s) {
	try {
		miint::Minimap2Aligner aligner(FilterMinimap2Config());
		RypeFilterUnit unit;
		{
			std::lock_guard<std::mutex> guard(s.lock);
			if (!s.ClaimToAlignLocked(unit)) {
				return;
			}
			aligner.attach_shared_index(s.minimap2_index);
		}
		for (;;) {
			unit.survivor_count = SelectSurvivors(s, &aligner, unit, unit.survivors);
			RypeFilterFollowUp follow_up;
			bool more;
			{
				std::lock_guard<std::mutex> guard(s.lock);
				follow_up = s.AlignedLocked(std::move(unit));
				more = s.ClaimToAlignLocked(unit);
			}
			Wake(s, follow_up.wake);
			if (!more) {
				return;
			}
		}
	} catch (std::exception &e) {
		Fail(s, ErrorData(e));
	}
}

enum class FilterJob { CLASSIFY, LOAD_INDEX, ALIGN, READ_AHEAD };

//! A job scheduled on the async pool. The schedule it runs against is its own
//! reference, so the job never depends on the GlobalState, which drains the
//! executor before it goes anyway.
class FilterJobTask : public BaseExecutorTask {
public:
	FilterJobTask(TaskExecutor &executor, shared_ptr<Schedule> schedule_p, FilterJob job_p)
	    : BaseExecutorTask(executor), schedule(std::move(schedule_p)), job(job_p) {
	}

	void ExecuteTask() override {
		// Each job settles its own work's errors. What escapes it, such as a failure
		// to schedule the jobs it starts, is settled here too: the executor keeps its
		// own record of an error, which the threads waiting on the schedule never read.
		try {
			switch (job) {
			case FilterJob::CLASSIFY:
				RunClassify(schedule, &executor);
				break;
			case FilterJob::LOAD_INDEX:
				RunIndexLoad(schedule, &executor);
				break;
			case FilterJob::ALIGN:
				RunAlign(*schedule);
				break;
			case FilterJob::READ_AHEAD:
				// Keeps a scan error for the classify job that reaches it.
				schedule->input_stream->ReadAhead();
				break;
			}
		} catch (std::exception &e) {
			Fail(*schedule, ErrorData(e));
			throw;
		}
	}

	//! Retired unrun. As the query is torn down, or after another job failed, it is
	//! settled as interrupted, so nothing waits on a job that will never run; an
	//! error already recorded is kept. A job that could not be queued records
	//! nothing, and the scheduling error is thrown to whoever scheduled it: a job
	//! records it on the schedule and wakes everyone, and a DuckDB thread throws it
	//! to DuckDB, which reports it and stops the other threads.
	void Cancel() override {
		if (executor.IsCancelled() || executor.HasError()) {
			Fail(*schedule, ErrorData(InterruptException()));
		}
	}

	string TaskType() const override {
		switch (job) {
		case FilterJob::CLASSIFY:
			return "RypeFilterClassify";
		case FilterJob::LOAD_INDEX:
			return "RypeFilterIndexLoad";
		case FilterJob::ALIGN:
			return "RypeFilterAlign";
		default:
			return "RypeFilterReadAhead";
		}
	}

private:
	shared_ptr<Schedule> schedule;
	FilterJob job;
};

//! Run the jobs a thread took: on the async pool, or here when it has no threads.
void StartJobs(TaskExecutor *executor, const shared_ptr<Schedule> &schedule, const RypeFilterJobs &jobs) {
	if (!executor) {
		// DuckDB's threads map what they emit: no align job is ever due.
		D_ASSERT(jobs.align == 0);
		if (jobs.load_index) {
			RunIndexLoad(schedule, nullptr);
		}
		if (jobs.classify) {
			RunClassify(schedule, nullptr);
		}
		return;
	}
	vector<unique_ptr<BaseExecutorTask>> tasks;
	if (jobs.load_index) {
		tasks.push_back(make_uniq<FilterJobTask>(*executor, schedule, FilterJob::LOAD_INDEX));
	}
	if (jobs.classify) {
		tasks.push_back(make_uniq<FilterJobTask>(*executor, schedule, FilterJob::CLASSIFY));
	}
	for (idx_t i = 0; i < jobs.align; i++) {
		tasks.push_back(make_uniq<FilterJobTask>(*executor, schedule, FilterJob::ALIGN));
	}
	executor->ScheduleTasks(std::move(tasks));
}

//! Slice the rows of `unit` that neither stage removes into `output`, mapping
//! them here unless an align job has. False when it has none, leaving `output`
//! untouched.
bool EmitUnit(const Schedule &s, miint::Minimap2Aligner *aligner, RypeFilterUnit &unit, DataChunk &output) {
	// A unit no align job mapped is mapped here, so with minimap2 on this thread
	// must have an aligner (InitLocal gives it one when there are no align jobs).
	D_ASSERT(unit.aligned || !s.use_minimap2 || aligner);
	if (!unit.aligned) {
		unit.survivor_count = SelectSurvivors(s, aligner, unit, unit.survivors);
	}
	if (unit.survivor_count == 0) {
		return false;
	}
	// The output references the held chunk's buffers, so the unit can be dropped
	// once it is emitted.
	output.Slice(*unit.group->runs[unit.run].chunk, unit.survivors, unit.survivor_count);
	return true;
}

//! Count `unit` emitted, and release its group with its last run.
void FinishUnit(RypeFilterTableFunction::GlobalState &gstate, const RypeFilterUnit &unit) {
	auto &s = *gstate.schedule;
	RypeFilterFollowUp follow_up;
	{
		std::lock_guard<std::mutex> guard(s.lock);
		follow_up = s.FinishUnitLocked(unit, !gstate.executor);
	}
	Wake(s, follow_up.wake);
	if (follow_up.jobs.Any()) {
		StartJobs(gstate.executor.get(), gstate.schedule, follow_up.jobs);
	}
}

} // namespace

// ============================================================================
// GlobalState destructor
// ============================================================================
RypeFilterTableFunction::GlobalState::~GlobalState() {
	// A query can end with groups still in flight, after a LIMIT or an error.
	// Settle the schedule as interrupted, so the align jobs stop at their next
	// unit; stop the RYpe call at its next batch, and the scan behind it at its
	// next fetch. The jobs are waited for below.
	if (schedule) {
		Fail(*schedule, ErrorData(InterruptException()));
	}
	if (schedule && schedule->input_stream) {
		schedule->input_stream->Abort();
	}
	if (schedule && schedule->input_connection) {
		schedule->input_connection->Interrupt();
	}
	if (executor) {
		executor->CancelAndDrain();
	}
	// No job runs now, though a finished one may still hold a reference: free what
	// the schedule holds here, on the query's thread.
	if (schedule) {
		schedule->Release();
	}
}

// ============================================================================
// Bind
// ============================================================================
unique_ptr<FunctionData> RypeFilterTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                       vector<LogicalType> &return_types, vector<Identifier> &names) {
	static const char *const ARGUMENTS[] = {"sequence_table", "index_path", "minimap2_index_path"};
	for (idx_t i = 0; i < input.inputs.size(); i++) {
		if (input.inputs[i].IsNull()) {
			throw InvalidInputException("rype_filter: %s cannot be NULL", ARGUMENTS[i]);
		}
	}
	auto data = make_uniq<Data>();
	data->relation = input.inputs[0].ToString();
	RejectCTERelationName(input, data->relation);
	data->index_path = input.inputs[1].ToString();
	if (input.inputs.size() > 2) {
		data->minimap2_index_path = input.inputs[2].ToString();
		if (data->minimap2_index_path.empty()) {
			throw InvalidInputException(
			    "rype_filter: minimap2_index_path cannot be empty; leave it out to filter with RYpe alone");
		}
	}

	data->threshold = ParseRypeThreshold(input);

	// The index paths are left to the loaders in InitGlobal, as rype_classify
	// leaves its own: they report a missing or damaged index better than a check
	// here, and a minimap2 index is not opened at all for long reads.
	auto columns = GetTableOrViewColumns(context, data->relation, "Table");
	for (idx_t i = 0; i < columns.names.size(); i++) {
		const auto &name = columns.names[i];
		if (StringUtil::CIEquals(name, "read_id")) {
			data->read_id_column = i;
			continue;
		}
		const bool is_sequence1 = StringUtil::CIEquals(name, "sequence1");
		if (!is_sequence1 && !StringUtil::CIEquals(name, "sequence2")) {
			continue;
		}
		(is_sequence1 ? data->sequence1_column : data->sequence2_column) = i;
		RequireRypeSequenceType(columns.types[i], is_sequence1 ? "sequence1" : "sequence2", data->relation);
	}
	if (data->sequence1_column == DConstants::INVALID_INDEX) {
		throw BinderException("Table '%s' missing required column 'sequence1'", data->relation);
	}
	// minimap2 is given the read_id as its name; RYpe alone never reads it.
	if (!data->minimap2_index_path.empty() && data->read_id_column != DConstants::INVALID_INDEX &&
	    !IsAllowedIdType(columns.types[data->read_id_column])) {
		throw BinderException("Column '%s' in '%s' must be %s to name its reads to minimap2, not %s",
		                      columns.names[data->read_id_column], data->relation, AllowedIdTypeList(),
		                      columns.types[data->read_id_column].ToString());
	}

	for (idx_t i = 0; i < columns.names.size(); i++) {
		names.emplace_back(columns.names[i]);
		return_types.push_back(columns.types[i]);
	}
	data->column_names = std::move(columns.names);
	data->column_types = std::move(columns.types);
	return data;
}

// ============================================================================
// InitGlobal
// ============================================================================
unique_ptr<GlobalTableFunctionState> RypeFilterTableFunction::InitGlobal(ClientContext &context,
                                                                         TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	gstate->schedule = make_shared_ptr<Schedule>();
	auto &s = *gstate->schedule;
	s.threshold = bind_data.threshold;
	const idx_t group_rows_setting = GetRypeFilterGroupRows(context);

	// The scan opens first: whether minimap2 runs is read from its head, and a
	// minimap2 index path that holds no index is refused before the RYpe index loads.
	auto &db = DatabaseInstance::GetDatabase(context);
	s.input_connection = make_uniq<Connection>(db);
	InheritTempObjects(context, *s.input_connection);
	auto &conn = *s.input_connection;
	ConfigureRypeArrowExport(conn);

	const idx_t batch_bytes = GetRypeArrowBatchBytes(context);

	RypeInputStreamOptions stream_options;
	stream_options.relation_quoted = SQLIdentifier::ToString(bind_data.relation);
	stream_options.source_name = bind_data.relation;
	stream_options.batch_bytes = batch_bytes;
	s.input_stream = BuildRypeFilterStream(conn, bind_data.column_names, std::move(stream_options));
	// The output schema was fixed at bind; a relation redefined since then would
	// hand rows of other types to a plan that expects these.
	if (s.input_stream->ScanTypes() != bind_data.column_types) {
		throw InvalidInputException("'%s' changed after the query was bound: its column types differ",
		                            bind_data.relation);
	}

	// From the head of the stream's own scan, which SampleSizing does not consume:
	// the relation is read once.
	const auto &minimap2_path = bind_data.minimap2_index_path;
	RypeInputStream::Sizing sizing {};
	sizing = s.input_stream->SampleSizing();
	// An empty relation has nothing to map, and no length to tell short reads by.
	const bool use_minimap2 =
	    !minimap2_path.empty() && sizing.sampled && sizing.avg_read_length < LONG_READ_MEAN_LENGTH;
	idx_t minimap2_index_bytes = 0;
	if (use_minimap2) {
		// Only whether it is an index at all, from its first bytes, as align_minimap2
		// checks: minimap2's reader would take any other file for FASTA and index it.
		// Whether it has more parts than one is known once the first has loaded.
		bool is_index = false;
		try {
			is_index = miint::Minimap2Aligner::is_index_file(minimap2_path);
			minimap2_index_bytes =
			    FileSystem::GetFileSystem(context).OpenFile(minimap2_path, FileFlags::FILE_FLAGS_READ)->GetFileSize();
		} catch (const std::exception &e) {
			throw IOException("Failed to load minimap2 index '%s': %s", minimap2_path, e.what());
		}
		if (!is_index) {
			throw IOException("Failed to load minimap2 index '%s': it is not a minimap2 index", minimap2_path);
		}
	}

	s.index = rype_index_load(bind_data.index_path.c_str());
	if (!s.index) {
		throw IOException("Failed to load RYpe index '%s': %s", bind_data.index_path,
		                  RypeErrorText(rype_get_last_error()));
	}

	if (group_rows_setting > 0) {
		s.group_rows = group_rows_setting;
	} else {
		// The minimap2 index is drawn from the same memory, and loads after this
		// reading is taken. 0 means detection failed, which leaves RYpe to size
		// against its own estimate.
		size_t budget = ResolveRypeMemoryBudget(context, 0);
		if (budget > 0) {
			const size_t floor = MinValue<size_t>(budget, RYPE_MIN_MEMORY_BUDGET);
			const size_t reserved = minimap2_index_bytes;
			budget = MaxValue<size_t>(budget > reserved ? budget - reserved : 0, floor);
			if (budget <= RYPE_MIN_MEMORY_BUDGET) {
				miint::EmitWarning(
				    context,
				    "rype_filter: only %s of memory is free for classification, so groups are small and "
				    "each costs a pass over the index, which slows the run. Free memory, or run when the "
				    "machine is quieter.",
				    StringUtil::BytesToHumanReadableString(budget));
			}
		}
		s.group_rows = AutoGroupRows(context, s.index, sizing, budget);
	}
	// The scan runs a group ahead of RYpe, once the async pool can read it ahead
	// (below), so the next group decodes while this one is classified.
	const idx_t read_ahead_bytes = s.group_rows * sizing.held_bytes_per_row;
	s.input_stream->SetGroupRows(s.group_rows);
	s.input_stream->SetBatchRows(s.group_rows);

	// The index loads on a job of its own, alongside the first groups' classifying.
	s.use_minimap2 = use_minimap2;
	s.minimap2_path = minimap2_path;
	s.sequence1_column = bind_data.sequence1_column;
	s.sequence2_column = bind_data.sequence2_column;
	s.read_id_column = bind_data.read_id_column;
	if (s.read_id_column != DConstants::INVALID_INDEX) {
		s.read_id_type = bind_data.column_types[s.read_id_column];
	}

	auto &scheduler = TaskScheduler::GetScheduler(context);
	gstate->max_threads = NumericCast<idx_t>(scheduler.NumberOfThreads());
	// The index loads while RYpe classifies the first groups, and no unit can be
	// mapped before it is in: it takes as many threads as DuckDB has.
	s.index_load_threads = scheduler.NumberOfThreads();
	if (scheduler.NumberOfAsyncThreads() > 0) {
		gstate->executor = make_uniq<TaskExecutor>(context, TaskSchedulerType::ASYNC);
		// minimap2 maps on the async pool too, as many units at once as DuckDB has
		// threads (see Schedule for why).
		if (s.use_minimap2) {
			s.max_aligners = gstate->max_threads;
		}
		// Decoding the relation overlaps classifying it only if something other than
		// RYpe's pulls fetches the scan (see SetReadAhead).
		weak_ptr<Schedule> weak_schedule = gstate->schedule;
		auto &executor = *gstate->executor;
		s.input_stream->SetReadAhead(read_ahead_bytes, [weak_schedule, &executor]() {
			// Asked for by the classify job, which holds the schedule, and which the
			// executor outlives.
			auto schedule = weak_schedule.lock();
			if (schedule) {
				executor.ScheduleTask(make_uniq<FilterJobTask>(executor, std::move(schedule), FilterJob::READ_AHEAD));
			}
		});
		// Start now, not at the first scan: the minimap2 index takes as long to load
		// as RYpe takes over several groups, and both can begin at once.
		RypeFilterJobs jobs;
		{
			std::lock_guard<std::mutex> guard(s.lock);
			jobs = s.TakeJobsLocked(false);
		}
		StartJobs(gstate->executor.get(), gstate->schedule, jobs);
	}
	return gstate;
}

// ============================================================================
// InitLocal
// ============================================================================
unique_ptr<LocalTableFunctionState> RypeFilterTableFunction::InitLocal(ExecutionContext &context,
                                                                       TableFunctionInitInput &input,
                                                                       GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto lstate = make_uniq<LocalState>();
	if (gstate.schedule->use_minimap2 && gstate.schedule->max_aligners == 0) {
		// Attached to the index once it has loaded, when the thread claims its first unit.
		lstate->aligner = make_uniq<miint::Minimap2Aligner>(FilterMinimap2Config());
	}
	return std::move(lstate);
}

// ============================================================================
// Execute
// ============================================================================
void RypeFilterTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();
	auto &s = *gstate.schedule;
	// Parking needs a task to wake later: there is none under a synchronous scan,
	// or where a positional join drives this one.
	const bool can_park =
	    data_p.results_execution_mode == AsyncResultsExecutionMode::TASK_EXECUTOR && data_p.interrupt_state;
	const bool inline_jobs = !gstate.executor;
	for (;;) {
		context.InterruptCheck();
		RypeFilterUnit unit;
		bool claimed = false;
		RypeFilterJobs jobs;
		bool park = false;
		idx_t wakes_seen = 0;
		{
			std::unique_lock<std::mutex> guard(s.lock);
			if (s.error.HasError()) {
				auto error = s.error;
				guard.unlock();
				error.Throw();
			}
			if (s.ClaimLocked(unit)) {
				claimed = true;
				if (lstate.aligner && !lstate.attached) {
					lstate.aligner->attach_shared_index(s.minimap2_index);
					lstate.attached = true;
				}
			} else if (s.DrainedLocked()) {
				// No unit is left, and no group will follow.
				data_p.async_result = AsyncResultType::FINISHED;
				return;
			} else {
				jobs = s.TakeJobsLocked(inline_jobs);
				// A thread that took jobs starts them, then looks again. Only a thread
				// that took none parks: starting a job can fail and wake the parked
				// threads, and waking its own registration would wait on this call's
				// return for good (Executor::RescheduleTask).
				if (can_park && !jobs.Any()) {
					// Registered last, under the lock every wake-up takes, so none is
					// missed; and nothing this call does after it wakes anyone.
					s.parked.push_back(*data_p.interrupt_state);
					park = true;
				}
				wakes_seen = s.wakes;
			}
		}

		if (claimed) {
			const bool emitted = EmitUnit(s, lstate.aligner.get(), unit, output);
			FinishUnit(gstate, unit);
			if (emitted) {
				data_p.async_result = AsyncResultType::HAVE_MORE_OUTPUT;
				return;
			}
			// A unit with no survivors: claim the next one in this same call.
			continue;
		}
		if (jobs.Any()) {
			StartJobs(gstate.executor.get(), gstate.schedule, jobs);
			continue;
		}
		if (park) {
			data_p.async_result = AsyncResultType::BLOCKED;
			return;
		}
		// Cannot park: wait for the background work, looking again now and then so
		// an interrupted query stops.
		std::unique_lock<std::mutex> guard(s.lock);
		s.progress.wait_for(guard, WAIT_FOR_PROGRESS, [&]() { return s.wakes != wakes_seen; });
	}
}

// ============================================================================
// Register
// ============================================================================
void RypeFilterTableFunction::Register(ExtensionLoader &loader) {
	TableFunctionSet set("rype_filter");
	// rype_filter(relation, rype_index) -> RYpe alone;
	// rype_filter(relation, rype_index, minimap2_index) -> RYpe, then minimap2 for short reads.
	for (const auto &args : {vector<LogicalType> {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                         vector<LogicalType> {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}}) {
		TableFunction tf("rype_filter", args, Execute, Bind, InitGlobal, InitLocal);
		AddNamedParameter(tf, "threshold", LogicalType::DOUBLE);
		// Each thread emits the units it claims, so rows come out in no fixed order,
		// as rype_classify's and align_minimap2's do. Kept in order, every unit would
		// be its own batch, which DuckDB's order-keeping sinks hold until enough rows
		// follow to merge them: a run with few survivors per unit would pay for that.
		tf.order_preservation_type = OrderPreservationType::NO_ORDER;
		set.AddFunction(tf);
	}
	loader.RegisterFunction(set);
}

} // namespace duckdb
