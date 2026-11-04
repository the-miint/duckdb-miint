#include <catch2/catch_test_macros.hpp>

#include "interruptible_sleep.hpp"

#include <chrono>

// The HTTP clients back off for seconds between retries, and BLAST waits minutes for a
// search. They wait through InterruptibleSleep with a check that throws on Ctrl-C or a
// passed max_execution_time, so the wait must poll that check while it lasts, not only
// when it ends.

namespace {

struct Interrupted {};

using std::chrono::milliseconds;
using std::chrono::steady_clock;

} // namespace

TEST_CASE("InterruptibleSleep polls its check every slice and lasts the whole wait", "[interruptible_sleep]") {
	int checks = 0;
	const auto start = steady_clock::now();
	miint::InterruptibleSleep(milliseconds(250), [&] { checks++; });
	// A backoff whose check never throws must still wait its full length before the retry.
	CHECK(steady_clock::now() - start >= milliseconds(250));
	// 100 + 100 + 50 ms.
	CHECK(checks == 3);
}

TEST_CASE("InterruptibleSleep checks after its last slice", "[interruptible_sleep]") {
	// The caller sends its next request as soon as the wait returns. A cancel during the
	// last slice must stop that request, which for a Webin submission is a POST that
	// creates objects, so the last check comes after all the sleeping.
	const auto start = steady_clock::now();
	auto last_check = start;
	miint::InterruptibleSleep(milliseconds(250), [&] { last_check = steady_clock::now(); });
	CHECK(last_check - start >= milliseconds(250));
}

TEST_CASE("InterruptibleSleep ends when its check throws", "[interruptible_sleep]") {
	int checks = 0;
	const auto start = steady_clock::now();
	REQUIRE_THROWS_AS(miint::InterruptibleSleep(milliseconds(60000),
	                                            [&] {
		                                            if (++checks == 3) {
			                                            throw Interrupted();
		                                            }
	                                            }),
	                  Interrupted);
	CHECK(checks == 3);
	// Three slices, not the minute; the bound is loose so a loaded machine cannot fail it.
	CHECK(steady_clock::now() - start < milliseconds(10000));
}
