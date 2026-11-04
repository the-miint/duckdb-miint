#pragma once

#include <algorithm>
#include <chrono>
#include <functional>
#include <thread>

namespace miint {

// The longest a wait runs on after its check would have thrown.
constexpr std::chrono::milliseconds INTERRUPT_POLL_INTERVAL {100};

// Sleeps for `duration` in slices of at most INTERRUPT_POLL_INTERVAL, calling `check`
// after each one. A check that throws ends the wait with its exception. The HTTP
// clients pass ClientContext::BlockingInterruptCheck, so Ctrl-C and max_execution_time
// stop a retry backoff or a BLAST wait within one slice instead of when it ends. The
// last check follows the last slice, so a cancel during the wait also stops the request
// the caller would send next.
inline void InterruptibleSleep(std::chrono::milliseconds duration, const std::function<void()> &check) {
	while (duration > std::chrono::milliseconds::zero()) {
		const auto slice = std::min(duration, INTERRUPT_POLL_INTERVAL);
		std::this_thread::sleep_for(slice);
		duration -= slice;
		check();
	}
}

} // namespace miint
