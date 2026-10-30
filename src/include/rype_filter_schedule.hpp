#pragma once
//
// What rype_filter's threads and background jobs share (rype_filter.cpp), and the
// state changes they make to it under its lock. The state changes do no work and
// wake no one themselves: each returns what its caller must do once the lock is
// released, which keeps them testable alone, without threads.

#include "rype_filter.hpp"

#include "duckdb/common/error_data.hpp"
#include "duckdb/parallel/interrupt.hpp"

#include <condition_variable>
#include <deque>
#include <mutex>

namespace duckdb {

//! One group RYpe classified: the runs the stream held for its rows, and which
//! of them RYpe matched. Read-only once published, except runs_left.
struct RypeFilterGroup {
	vector<RypeHeldRows> runs;
	//! The Arrow id of the group's first row; matched is indexed by id - first_id.
	idx_t first_id = 0;
	vector<bool> matched;
	//! Runs not yet emitted. Guarded by the schedule's lock.
	idx_t runs_left = 0;
};

//! One held run of a group, emitted by whichever thread claims it.
struct RypeFilterUnit {
	shared_ptr<RypeFilterGroup> group;
	idx_t run = 0;
	//! Set by the align job that mapped the run: its rows that neither stage
	//! removes, as positions in the run's chunk.
	bool aligned = false;
	SelectionVector survivors;
	idx_t survivor_count = 0;
};

//! Background work a thread has taken on under the lock, and must now run or schedule.
struct RypeFilterJobs {
	bool classify = false;
	bool load_index = false;
	//! How many align jobs to start.
	idx_t align = 0;

	bool Any() const {
		return classify || load_index || align > 0;
	}
};

//! What a state change leaves its caller to do once the lock is released: the
//! jobs it took on, and the parked threads to wake. A parked task's callback waits
//! for that task to return, and it may be about to take the lock.
struct RypeFilterFollowUp {
	RypeFilterJobs jobs;
	vector<InterruptState> wake;
};

//! Jobs run as tasks on DuckDB's async pool, apart from the query's pipeline
//! (DuckDB's idle threads take such tasks too): classifying, one RYpe call
//! per group, one group at a time and in scan order; and, when minimap2 runs,
//! loading its index, and mapping. Reading the scan ahead of RYpe is a job too,
//! asked for by the input stream itself and kept out of this state
//! (RypeInputStream::SetReadAhead). Each group's runs are published as units.
//! When minimap2 runs, align jobs map them, as many jobs as DuckDB has threads,
//! and DuckDB's threads claim and emit the mapped units: a sink that takes rows
//! in order, as miint's FASTQ COPY does, runs this function on one DuckDB thread,
//! which would otherwise do all the mapping. The output is declared unordered,
//! as rype_classify's and align_minimap2's are, so the threads' rows interleave
//! freely. A thread with nothing to claim parks where it can, and parked threads
//! are woken when units become ready, at the end of the input, and for an error.
//! Where a thread cannot park (a synchronous scan, or a positional join driving
//! this one) it waits on `progress` instead.
//! With no async threads the jobs run on the DuckDB thread that finds them
//! needed, DuckDB's threads map the units they claim, and threads are also woken
//! when a group is released, since one of them must then classify the next.
struct RypeFilterTableFunction::Schedule {
	//! At most this many groups are held at once: one whose runs are being mapped
	//! and emitted while the next is classified. The scan reads a group's rows
	//! ahead of them besides, so that decoding overlaps classifying.
	static constexpr idx_t MAX_HELD_GROUPS = 2;

	//! The rows of a group, given the batch RYpe would classify within the budget
	//! (batch_size reads in per_batch_memory bytes) and what holding a row costs:
	//! one group's classification and the rows of MAX_HELD_GROUPS groups and the
	//! one read ahead fit in per_batch_memory. At least one row, and no more than
	//! RYpe's batch.
	static idx_t GroupRows(idx_t per_batch_memory, idx_t batch_size, idx_t held_bytes_per_row);

	// Set at InitGlobal, then read-only.
	RypeIndex *index = nullptr;
	//! The connection the scan runs on. The stream's streaming result points into
	//! its context, so it outlives input_stream.
	unique_ptr<Connection> input_connection;
	unique_ptr<RypeInputStream> input_stream;
	double threshold = RYPE_DEFAULT_THRESHOLD;
	idx_t group_rows = 0;
	bool use_minimap2 = false;
	std::string minimap2_path;
	//! Threads that build the minimap2 index's hash tables as it loads.
	int index_load_threads = 1;
	//! Align jobs map at most this many units at once. 0 when minimap2 does not
	//! run, and when DuckDB's threads map the units they emit (no async threads).
	idx_t max_aligners = 0;
	//! Where minimap2's inputs are among the scanned columns (see Data); the
	//! align jobs read them from here.
	idx_t sequence1_column = DConstants::INVALID_INDEX;
	idx_t sequence2_column = DConstants::INVALID_INDEX;
	idx_t read_id_column = DConstants::INVALID_INDEX;
	LogicalType read_id_type;

	//! The Arrow id the next group starts at. Only the one classify running reads
	//! or writes it.
	idx_t next_first_id = 0;

	std::mutex lock;
	std::condition_variable progress;
	// Guarded by lock.
	vector<InterruptState> parked;
	//! Counts wake-ups, so a thread waiting on `progress` sees one it slept through.
	idx_t wakes = 0;
	//! Published runs not yet claimed, in scan order.
	std::deque<RypeFilterUnit> units;
	//! Units mapped by an align job, not yet claimed by a DuckDB thread.
	std::deque<RypeFilterUnit> aligned;
	//! Align jobs running, and the units they have claimed and not yet handed on.
	idx_t aligners = 0;
	idx_t units_aligning = 0;
	//! Groups being classified, or with runs not yet emitted.
	idx_t held_groups = 0;
	bool classifying = false;
	//! The scan has no rows left: no group will follow those published.
	bool input_done = false;
	bool index_loading = false;
	std::shared_ptr<miint::SharedMinimap2Index> minimap2_index;
	ErrorData error;

	// The state changes. Each is called with `lock` held.

	//! The jobs that are due, marked as running. Run inline, a thread takes one,
	//! so that another thread can take the other.
	RypeFilterJobs TakeJobsLocked(bool inline_jobs);
	//! A classify job's group, or null for the end of the input. The job goes on to
	//! the next group when the follow-up's classify is set.
	RypeFilterFollowUp PublishLocked(shared_ptr<RypeFilterGroup> group);
	RypeFilterFollowUp IndexLoadedLocked(std::shared_ptr<miint::SharedMinimap2Index> loaded);
	//! Record the first error; everyone woken throws it.
	RypeFilterFollowUp FailLocked(ErrorData failure);
	//! A DuckDB thread takes the next unit to emit: one an align job mapped, or,
	//! when there are none, one published once minimap2, if it runs, has its index.
	bool ClaimLocked(RypeFilterUnit &unit);
	//! No unit is left to emit, and none will follow: a DuckDB thread can finish.
	bool DrainedLocked() const;
	//! An align job takes the next published unit to map. False when none is left
	//! or the query is failing: the job then ends, and no longer counts as running.
	bool ClaimToAlignLocked(RypeFilterUnit &unit);
	//! An align job mapped `unit`: it is ready to emit.
	RypeFilterFollowUp AlignedLocked(RypeFilterUnit unit);
	//! `unit` is emitted. With its group's last unit the group is released, and the
	//! next group's classify is due: taken here when async threads will run it, and
	//! otherwise left to a woken thread.
	RypeFilterFollowUp FinishUnitLocked(const RypeFilterUnit &unit, bool inline_jobs);

	//! Free the RYpe index, the held rows, the scan and its connection, and drop
	//! the schedule's reference to the minimap2 index (an aligner a DuckDB thread
	//! holds keeps its own). The GlobalState calls it once the jobs are drained: a job's task
	//! is dropped, with its reference to the schedule, only after the executor has
	//! counted it finished, so the last reference can go on any of the pool's
	//! threads, after the query. ~Schedule calls it too; a second call does nothing.
	void Release();

	~Schedule();

private:
	bool TakeClassifyLocked();
	bool TakeIndexLoadLocked();
	//! Align jobs to start for the published units: no more than there are units,
	//! nor than the limit allows beside those running.
	idx_t TakeAlignLocked();
	//! Every parked thread, and a wake-up for those waiting on `progress`.
	vector<InterruptState> TakeWaitersLocked();
};

} // namespace duckdb
