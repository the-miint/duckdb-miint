#pragma once
/*
 * Shared debug logging utilities for sharded alignment table functions.
 * Used by align_minimap2_sharded and align_bowtie2_sharded.
 */

#include <chrono>
#include <cstdio>
#include <cstring>
#include <thread>

namespace duckdb {

// Read current process RSS from /proc/self/status (Linux only, returns MB)
inline long GetRSSMB() {
	FILE *f = fopen("/proc/self/status", "r");
	if (!f) {
		return -1;
	}
	long rss_kb = -1;
	char line[256];
	while (fgets(line, sizeof(line), f)) {
		if (strncmp(line, "VmRSS:", 6) == 0) {
			sscanf(line + 6, "%ld", &rss_kb);
			break;
		}
	}
	fclose(f);
	return rss_kb > 0 ? rss_kb / 1024 : -1;
}

} // namespace duckdb

// NOLINTNEXTLINE: macro for debug logging with elapsed time and thread ID
#define SHARD_DBG(gstate, ...)                                                                                         \
	do {                                                                                                               \
		if ((gstate).debug) {                                                                                          \
			auto _elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() -   \
			                                                                      (gstate).start_time);                \
			auto _tid = std::hash<std::thread::id> {}(std::this_thread::get_id()) % 10000;                             \
			fprintf(stderr, "[%7ldms T%04zu] ", static_cast<long>(_elapsed.count()), _tid);                            \
			fprintf(stderr, __VA_ARGS__);                                                                              \
			fprintf(stderr, "\n");                                                                                     \
		}                                                                                                              \
	} while (0)

// Like SHARD_DBG but appends current RSS
#define SHARD_DBG_MEM(gstate, ...)                                                                                     \
	do {                                                                                                               \
		if ((gstate).debug) {                                                                                          \
			auto _elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() -   \
			                                                                      (gstate).start_time);                \
			auto _tid = std::hash<std::thread::id> {}(std::this_thread::get_id()) % 10000;                             \
			fprintf(stderr, "[%7ldms T%04zu] ", static_cast<long>(_elapsed.count()), _tid);                            \
			fprintf(stderr, __VA_ARGS__);                                                                              \
			fprintf(stderr, " [RSS=%ldMB]\n", duckdb::GetRSSMB());                                                     \
		}                                                                                                              \
	} while (0)
