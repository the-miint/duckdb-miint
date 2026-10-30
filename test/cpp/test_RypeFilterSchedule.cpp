// rype_filter's schedule bookkeeping, driven without RYpe, minimap2 or a query.
//
// Which background job is due, and which parked threads to wake, are decided under the schedule's
// lock by a few state changes (rype_filter_schedule.hpp). Most of them show in a query's result
// when wrong: a lost wake-up hangs the query, and a unit dropped or emitted twice changes its rows
// (test/sql/rype_filter_parallel.test). These do not, or only now and then, and are pinned here:
//  - At most two groups are held: one being mapped while the next is classified. Classifying is
//    faster than mapping, so without the cap the classify job would run ahead through the whole
//    input, holding every read it passed in memory, and the rows would still come out right.
//  - The classify job goes straight on to the second group. Otherwise the next classify would
//    start only once a thread ran out of units, and every thread would then wait on it: mapping
//    and classifying would take turns instead of overlapping, again with the right rows.
//  - Releasing a group starts the next one's classify. Besides keeping the two overlapped, this is
//    what wakes the parked threads when the thread that released the group stops pulling rows, as
//    it does once a LIMIT is met: with no classify started they wait for one forever. Whether a
//    query meets that depends on timing, so its hang shows only in some runs.
//  - With no async threads, releasing a group wakes the parked threads instead, since one of them
//    must classify the next group itself. Missing that hangs a query the same way, as rarely.
//  - With async threads, minimap2 maps on align jobs there, as many as DuckDB has threads, and
//    DuckDB's threads only emit what they mapped. A sink that takes rows in order (miint's FASTQ
//    COPY) runs this function on one DuckDB thread, so mapping there would take one core: the rows
//    would be right, many times slower.
//  - An align job goes on until no published unit is left, and the next publish starts jobs again.
//    A unit published with no job to map it waits forever, which hangs the query only when the
//    last job stops just before that publish.
//  - A DuckDB thread finishes only once no unit is still being mapped. Every thread finishing while
//    a job maps the last units drops their rows, and only when the jobs run late.
//  - A group is sized so that three groups' rows fit beside RYpe's own cost: the one being emitted,
//    the one RYpe classifies, and the next, read ahead so that decoding it overlaps classifying.
//    Sized for two, the run holds a group's rows more than its budget allows, and nothing fails.

#include <catch2/catch_test_macros.hpp>

#include "rype_filter_schedule.hpp"

namespace {

using namespace duckdb;
using Schedule = RypeFilterTableFunction::Schedule;

//! A classified group of `runs` runs, as the classify job publishes it. The runs hold no rows:
//! the bookkeeping only counts them.
shared_ptr<RypeFilterGroup> Group(idx_t runs) {
	auto group = make_shared_ptr<RypeFilterGroup>();
	group->runs.resize(runs);
	return group;
}

vector<RypeFilterUnit> ClaimAll(Schedule &s) {
	vector<RypeFilterUnit> claimed;
	RypeFilterUnit unit;
	while (s.ClaimLocked(unit)) {
		claimed.push_back(std::move(unit));
	}
	return claimed;
}

//! A loaded minimap2 index, as far as the bookkeeping can tell: it is only ever handed on.
std::shared_ptr<miint::SharedMinimap2Index> Index() {
	return std::make_shared<miint::SharedMinimap2Index>(miint::Minimap2IndexPtr(), mm_mapopt_t {},
	                                                    std::vector<std::string>());
}

//! A schedule that maps on up to `aligners` async jobs, with its index loaded and its classify job
//! running.
void StartMapping(Schedule &s, idx_t aligners) {
	s.use_minimap2 = true;
	s.max_aligners = aligners;
	REQUIRE(s.TakeJobsLocked(false).classify);
	REQUIRE(s.IndexLoadedLocked(Index()).jobs.align == 0);
}

} // namespace

TEST_CASE("rype_filter classifies one group ahead of the one being mapped, and no further", "[rype_filter]") {
	Schedule s;
	REQUIRE(s.TakeJobsLocked(false).classify);
	// The job publishes the first group and goes straight on to the second.
	REQUIRE(s.PublishLocked(Group(2)).jobs.classify);
	// With two groups held it stops, and no thread can start a third.
	REQUIRE_FALSE(s.PublishLocked(Group(1)).jobs.classify);
	REQUIRE_FALSE(s.TakeJobsLocked(false).classify);
	// Claiming units releases nothing: they are held until emitted.
	auto units = ClaimAll(s);
	REQUIRE(units.size() == 3);
	REQUIRE_FALSE(s.TakeJobsLocked(false).classify);
}

TEST_CASE("rype_filter starts the next group's classify as a group is released", "[rype_filter]") {
	Schedule s;
	REQUIRE(s.TakeJobsLocked(false).classify);
	REQUIRE(s.PublishLocked(Group(2)).jobs.classify);
	REQUIRE_FALSE(s.PublishLocked(Group(1)).jobs.classify);
	auto units = ClaimAll(s);
	REQUIRE(units.size() == 3);
	// The first group's first unit leaves it held.
	REQUIRE_FALSE(s.FinishUnitLocked(units[0], false).jobs.Any());
	// Its last releases it: the third group's classify is due at once, while the second group's unit
	// is still being mapped, and it is handed to the thread that released the group.
	REQUIRE(s.FinishUnitLocked(units[1], false).jobs.classify);
	REQUIRE_FALSE(s.TakeJobsLocked(false).classify);
}

TEST_CASE("rype_filter without async threads wakes the parked threads as a group is released", "[rype_filter]") {
	Schedule s;
	REQUIRE(s.TakeJobsLocked(true).classify);
	REQUIRE(s.PublishLocked(Group(1)).jobs.classify);
	REQUIRE_FALSE(s.PublishLocked(Group(1)).jobs.classify);
	auto units = ClaimAll(s);
	REQUIRE(units.size() == 2);
	s.parked.push_back(InterruptState());
	auto follow_up = s.FinishUnitLocked(units[0], true);
	// The job is not the releasing thread's: a woken thread takes it, and runs it itself.
	REQUIRE_FALSE(follow_up.jobs.Any());
	REQUIRE(follow_up.wake.size() == 1);
	REQUIRE(s.parked.empty());
	REQUIRE(s.TakeJobsLocked(true).classify);
}

TEST_CASE("rype_filter maps on as many align jobs as DuckDB has threads, once the index is in", "[rype_filter]") {
	Schedule s;
	s.use_minimap2 = true;
	s.max_aligners = 2;
	auto jobs = s.TakeJobsLocked(false);
	REQUIRE(jobs.load_index);
	REQUIRE(jobs.classify);
	// Published before the index is in, the units wait for it.
	REQUIRE(s.PublishLocked(Group(3)).jobs.align == 0);
	// Then two jobs map the three, not one per unit.
	REQUIRE(s.IndexLoadedLocked(Index()).jobs.align == 2);
	// The next group's units go to those two, which are still running.
	REQUIRE(s.PublishLocked(Group(2)).jobs.align == 0);
	// DuckDB's threads emit only what has been mapped.
	RypeFilterUnit unit;
	REQUIRE_FALSE(s.ClaimLocked(unit));
}

TEST_CASE("rype_filter's mapped units go to a DuckDB thread, which is woken for them", "[rype_filter]") {
	Schedule s;
	StartMapping(s, 1);
	REQUIRE(s.PublishLocked(Group(2)).jobs.align == 1);
	RypeFilterUnit unit;
	REQUIRE(s.ClaimToAlignLocked(unit));
	REQUIRE(unit.run == 0);
	s.parked.push_back(InterruptState());
	auto follow_up = s.AlignedLocked(std::move(unit));
	REQUIRE(follow_up.wake.size() == 1);
	RypeFilterUnit emitted;
	REQUIRE(s.ClaimLocked(emitted));
	REQUIRE(emitted.run == 0);
	// The job maps the group's other unit next; nothing else is ready meanwhile.
	REQUIRE_FALSE(s.ClaimLocked(emitted));
	REQUIRE(s.ClaimToAlignLocked(unit));
	REQUIRE(unit.run == 1);
}

TEST_CASE("rype_filter's align jobs stop when they run out, and the next group starts them again", "[rype_filter]") {
	Schedule s;
	StartMapping(s, 2);
	// No more jobs than units.
	REQUIRE(s.PublishLocked(Group(1)).jobs.align == 1);
	RypeFilterUnit unit;
	REQUIRE(s.ClaimToAlignLocked(unit));
	s.AlignedLocked(std::move(unit));
	// Nothing left: the job stops, and no longer counts as running.
	REQUIRE_FALSE(s.ClaimToAlignLocked(unit));
	REQUIRE(s.PublishLocked(Group(3)).jobs.align == 2);
}

TEST_CASE("rype_filter's DuckDB threads finish only once no unit is being mapped", "[rype_filter]") {
	Schedule s;
	StartMapping(s, 1);
	REQUIRE(s.PublishLocked(Group(1)).jobs.align == 1);
	s.PublishLocked(nullptr);
	RypeFilterUnit unit;
	REQUIRE(s.ClaimToAlignLocked(unit));
	// The input is done and nothing is queued, but the job's unit is still to come.
	REQUIRE_FALSE(s.DrainedLocked());
	s.AlignedLocked(std::move(unit));
	REQUIRE_FALSE(s.DrainedLocked());
	RypeFilterUnit emitted;
	REQUIRE(s.ClaimLocked(emitted));
	// The thread that claimed it emits it; the others can finish.
	REQUIRE(s.DrainedLocked());
}

TEST_CASE("rype_filter's align jobs stop at a failure", "[rype_filter]") {
	Schedule s;
	StartMapping(s, 2);
	REQUIRE(s.PublishLocked(Group(3)).jobs.align == 2);
	RypeFilterUnit unit;
	REQUIRE(s.ClaimToAlignLocked(unit));
	s.FailLocked(ErrorData(IOException("mapping failed")));
	// The query is failing: nobody would emit what they map.
	REQUIRE_FALSE(s.ClaimToAlignLocked(unit));
	REQUIRE(s.PublishLocked(Group(1)).jobs.align == 0);
}

TEST_CASE("rype_filter without async threads maps on DuckDB's threads", "[rype_filter]") {
	Schedule s;
	s.use_minimap2 = true;
	// The index load is the one job a thread takes inline; the next thread takes the classify.
	REQUIRE(s.TakeJobsLocked(true).load_index);
	REQUIRE(s.TakeJobsLocked(true).classify);
	REQUIRE(s.PublishLocked(Group(1)).jobs.align == 0);
	RypeFilterUnit unit;
	REQUIRE_FALSE(s.ClaimLocked(unit));
	REQUIRE(s.IndexLoadedLocked(Index()).jobs.align == 0);
	REQUIRE(s.ClaimLocked(unit));
}

TEST_CASE("rype_filter sizes a group so that three groups' rows fit beside RYpe's cost", "[rype_filter]") {
	// RYpe would classify 1000 reads in 1 MB, 1000 bytes a read, and a held row costs 100.
	// 769 reads cost RYpe 769,000 bytes, and three groups of them hold 230,700: 999,700.
	REQUIRE(Schedule::GroupRows(1000000, 1000, 100) == 769);
	// Rows that cost nothing to hold leave RYpe its whole batch.
	REQUIRE(Schedule::GroupRows(1000000, 1000, 0) == 1000);
	// Rows too large to hold three of still make groups of one.
	REQUIRE(Schedule::GroupRows(1000000, 1000, 1000000) == 1);
}
