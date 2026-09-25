#pragma once

// Single-threaded stand-ins for the vsearch entry points that always start threads.
//
// dust_all(), search_batch(), chimera_detect_batch() and cluster_assign_batch()
// each pthread_create opt_threads workers, even when opt_threads is 1. Where no
// thread can be started (a wasm build without -pthread) that is fatal: vsearch's
// fatal() calls exit(1). Each stand-in takes the same arguments as the function it
// replaces and does the same work on the calling thread, in input order, calling
// the per-item function each worker calls. Wrappers select them through their
// params' `serial` flag (default: kVsearchSerialByDefault in vsearch_utils.hpp).
// They are built everywhere so the native C++ tests (test_VsearchSerial.cpp) can
// require that they return exactly what the threaded functions return.

#include "vsearch_api.h"

#include <cstddef>
#include <cstdint>

namespace miint {

// dust_all(): its workers dust one database sequence at a time.
inline void DustAllSerial() {
	for (uint64_t seqno = 0; seqno < db_getsequencecount(); seqno++) {
		dust(db_getsequence(seqno), static_cast<int>(db_getsequencelen(seqno)));
	}
}

// search_batch(): its workers run each query through their own search state;
// here one session serves every query in turn.
inline void SearchBatchSerial(const char **query_seqs, const char **query_heads, const int *query_lens,
                              const int *query_sizes, int query_count, search_result_s *results,
                              int max_results_per_query, int *result_counts) {
	search_session_s *ss = search_session_alloc();
	search_session_init(ss);
	for (int qi = 0; qi < query_count; qi++) {
		search_session_single(ss, query_seqs[qi], query_heads[qi], query_lens[qi], query_sizes[qi],
		                      results + static_cast<size_t>(qi) * max_results_per_query, max_results_per_query,
		                      &result_counts[qi]);
	}
	search_session_cleanup(ss);
	search_session_free(ss);
}

// chimera_detect_batch(): the same session bracket, including saving and
// restoring the option globals chimera_session_init() overwrites, around one
// chimera_detect_single() per query.
inline void ChimeraDetectBatchSerial(const char **query_seqs, const char **query_heads, const int *query_lens,
                                     const int *query_sizes, int query_count, chimera_result_s *results) {
	if (query_count <= 0) {
		return;
	}
	auto const saved_maxaccepts = opt_maxaccepts;
	auto const saved_maxrejects = opt_maxrejects;
	auto const saved_id = opt_id;
	auto const saved_weak_id = opt_weak_id;
	auto const saved_self = opt_self;
	auto const saved_selfid = opt_selfid;
	auto const saved_maxsizeratio = opt_maxsizeratio;

	chimera_session_init();
	chimera_info_s *ci = chimera_info_alloc();
	chimera_detect_thread_init(ci);
	for (int qi = 0; qi < query_count; qi++) {
		chimera_detect_single(ci, query_seqs[qi], query_heads[qi], query_lens[qi], query_sizes[qi], &results[qi]);
	}
	chimera_detect_thread_cleanup(ci);
	chimera_info_free(ci);
	chimera_session_cleanup();

	opt_maxaccepts = saved_maxaccepts;
	opt_maxrejects = saved_maxrejects;
	opt_id = saved_id;
	opt_weak_id = saved_weak_id;
	opt_self = saved_self;
	opt_selfid = saved_selfid;
	opt_maxsizeratio = saved_maxsizeratio;
}

// cluster_assign_batch(): assigns sequences start_seqno .. start_seqno+count-1
// one at a time, in order, as cluster_assign_single() requires.
inline void ClusterAssignBatchSerial(cluster_session_s *cs, int start_seqno, int count, cluster_result_s *results) {
	for (int i = 0; i < count; i++) {
		cluster_assign_single(cs, start_seqno + i, &results[i]);
	}
}

} // namespace miint
