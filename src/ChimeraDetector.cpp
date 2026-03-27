#include "ChimeraDetector.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace miint {

static bool is_gap(char c) {
	return c == '-';
}

// Returns true for IUPAC ambiguity codes (N, R, Y, S, W, K, M, B, D, H, V)
// and any non-ACGT non-gap character. Gap ('-') is NOT ambiguous — it has
// its own handling in classify_diffs.
static bool is_ambiguous(char c) {
	char upper = std::toupper(static_cast<unsigned char>(c));
	switch (upper) {
	case 'A':
	case 'C':
	case 'G':
	case 'T':
		return false;
	default:
		return !is_gap(c);
	}
}

static char to_upper(char c) {
	return std::toupper(static_cast<unsigned char>(c));
}

StarAlignment build_star_alignment(const std::string &query_aligned_a, const std::string &subject_aligned_a,
                                   const std::string &query_aligned_b, const std::string &subject_aligned_b) {
	StarAlignment result;

	// Walk both pairwise alignments using ungapped query positions as the anchor.
	// Both alignments are global (alignEnd2End) on the same query, so they consume
	// exactly the same query bases. The only structural difference is where gaps
	// appear (insertions in one subject but not the other).
	//
	// When one alignment has a gap in the query (= insertion in that subject),
	// the other arm gets a gap column inserted to maintain positional correspondence.

	size_t ia = 0; // position in alignment A
	size_t ib = 0; // position in alignment B

	while (ia < query_aligned_a.size() && ib < query_aligned_b.size()) {
		bool a_query_gap = is_gap(query_aligned_a[ia]);
		bool b_query_gap = is_gap(query_aligned_b[ib]);

		if (a_query_gap && b_query_gap) {
			// Both have insertion at this position
			result.query_row += '-';
			result.parent_a_row += subject_aligned_a[ia];
			result.parent_b_row += subject_aligned_b[ib];
			ia++;
			ib++;
		} else if (a_query_gap) {
			// Insertion in A's subject only — gap in B
			result.query_row += '-';
			result.parent_a_row += subject_aligned_a[ia];
			result.parent_b_row += '-';
			ia++;
		} else if (b_query_gap) {
			// Insertion in B's subject only — gap in A
			result.query_row += '-';
			result.parent_a_row += '-';
			result.parent_b_row += subject_aligned_b[ib];
			ib++;
		} else {
			// Both have a query base — emit the column
			result.query_row += query_aligned_a[ia];
			result.parent_a_row += subject_aligned_a[ia];
			result.parent_b_row += subject_aligned_b[ib];
			ia++;
			ib++;
		}
	}

	// Drain remaining insertions from whichever alignment hasn't finished.
	// After the main loop, one alignment may have trailing insertions (query gaps).
	// These are subject-only insertions that must be represented.
	while (ia < query_aligned_a.size()) {
		// Remaining columns in A must be query gaps (subject insertions)
		if (!is_gap(query_aligned_a[ia])) {
			throw std::runtime_error("build_star_alignment: alignment A has unconsumed query bases after B finished");
		}
		result.query_row += '-';
		result.parent_a_row += subject_aligned_a[ia];
		result.parent_b_row += '-';
		ia++;
	}
	while (ib < query_aligned_b.size()) {
		if (!is_gap(query_aligned_b[ib])) {
			throw std::runtime_error("build_star_alignment: alignment B has unconsumed query bases after A finished");
		}
		result.query_row += '-';
		result.parent_a_row += '-';
		result.parent_b_row += subject_aligned_b[ib];
		ib++;
	}

	return result;
}

void classify_diffs(StarAlignment &star) {
	size_t len = star.query_row.size();
	star.diffs.resize(len, DiffType::IGNORE);

	for (size_t i = 0; i < len; i++) {
		char q = to_upper(star.query_row[i]);
		char a = to_upper(star.parent_a_row[i]);
		char b = to_upper(star.parent_b_row[i]);

		// Skip if any position has a gap
		if (is_gap(q) || is_gap(a) || is_gap(b)) {
			star.diffs[i] = DiffType::IGNORE;
			continue;
		}

		// Skip if adjacent to a gap (i-1 or i+1 has a gap in any row).
		// Gap chars ('-') are case-insensitive so raw string access is safe here.
		bool adjacent_gap = false;
		if (i > 0) {
			adjacent_gap =
			    is_gap(star.query_row[i - 1]) || is_gap(star.parent_a_row[i - 1]) || is_gap(star.parent_b_row[i - 1]);
		}
		if (!adjacent_gap && i + 1 < len) {
			adjacent_gap =
			    is_gap(star.query_row[i + 1]) || is_gap(star.parent_a_row[i + 1]) || is_gap(star.parent_b_row[i + 1]);
		}
		if (adjacent_gap) {
			star.diffs[i] = DiffType::IGNORE;
			continue;
		}

		// Skip if any base is ambiguous (N, R, Y, etc.)
		if (is_ambiguous(q) || is_ambiguous(a) || is_ambiguous(b)) {
			star.diffs[i] = DiffType::IGNORE;
			continue;
		}

		// Classify based on which parent the query matches
		if (a == b) {
			if (q == a) {
				star.diffs[i] = DiffType::IGNORE; // All identical
			} else {
				star.diffs[i] = DiffType::NO_VOTE; // Both parents agree, query differs
			}
		} else {
			// Parents differ
			if (q == a) {
				star.diffs[i] = DiffType::MATCH_A;
			} else if (q == b) {
				star.diffs[i] = DiffType::MATCH_B;
			} else {
				star.diffs[i] = DiffType::ABSTAIN;
			}
		}
	}
}

// ============================================================================
// Phase 4: Smoothed parent selection
// ============================================================================

std::vector<int> compute_match_profile(const std::string &query_aligned, const std::string &subject_aligned) {
	size_t len = query_aligned.size();
	std::vector<int> profile(len, 0);
	for (size_t i = 0; i < len; i++) {
		if (!is_gap(query_aligned[i]) && !is_gap(subject_aligned[i]) &&
		    to_upper(query_aligned[i]) == to_upper(subject_aligned[i])) {
			profile[i] = 1;
		}
	}
	return profile;
}

std::vector<int> compute_smoothed(const std::vector<int> &match_profile, int window_size) {
	size_t len = match_profile.size();
	std::vector<int> smoothed(len, 0);
	if (len == 0) {
		return smoothed;
	}

	// Running sum
	int sum = 0;
	for (size_t i = 0; i < len; i++) {
		sum += match_profile[i];
		if (i >= static_cast<size_t>(window_size)) {
			sum -= match_profile[i - window_size];
		}
		smoothed[i] = sum;
	}
	return smoothed;
}

std::optional<ParentPair> select_parents(const std::string &query, const std::vector<uint32_t> &candidate_indices,
                                         const std::vector<std::string> &ref_sequences, WFA2Aligner &aligner) {
	if (candidate_indices.size() < 2) {
		return std::nullopt;
	}

	// Align query to each candidate and compute smoothed identity profiles.
	struct CandidateAlignment {
		uint32_t idx;
		WFA2FullResult alignment;
		std::vector<int> smoothed;
	};

	std::vector<CandidateAlignment> candidates;
	candidates.reserve(candidate_indices.size());

	for (uint32_t cidx : candidate_indices) {
		auto result = aligner.align_full(query, ref_sequences[cidx]);
		if (!result.has_value()) {
			continue;
		}
		auto match = compute_match_profile(result->query_aligned, result->subject_aligned);
		auto smooth = compute_smoothed(match, SMOOTHING_WINDOW);
		candidates.push_back({cidx, std::move(*result), std::move(smooth)});
	}

	if (candidates.size() < 2) {
		return std::nullopt;
	}

	// Use the minimum smoothed profile length across all candidates.
	// Alignments may differ in length due to indels — don't read past the shortest.
	size_t align_len = candidates[0].smoothed.size();
	for (size_t c = 1; c < candidates.size(); c++) {
		align_len = std::min(align_len, candidates[c].smoothed.size());
	}

	// Start win-counting at the first position where the smoothing window is full.
	// Positions 0..SMOOTHING_WINDOW-2 have partial sums and are not comparable.
	size_t win_start = (align_len >= static_cast<size_t>(SMOOTHING_WINDOW)) ? SMOOTHING_WINDOW - 1 : 0;

	// Round 1: find Parent A — candidate with most "wins" (positions where it has max smooth)
	std::vector<int> wins(candidates.size(), 0);
	for (size_t pos = win_start; pos < align_len; pos++) {
		int max_smooth = 0;
		for (size_t c = 0; c < candidates.size(); c++) {
			max_smooth = std::max(max_smooth, candidates[c].smoothed[pos]);
		}
		if (max_smooth == 0) {
			continue;
		}
		for (size_t c = 0; c < candidates.size(); c++) {
			if (candidates[c].smoothed[pos] == max_smooth) {
				wins[c]++;
			}
		}
	}

	// Parent A = candidate with most wins. Ties broken by lowest index.
	size_t parent_a_local = 0;
	for (size_t c = 1; c < candidates.size(); c++) {
		if (wins[c] > wins[parent_a_local]) {
			parent_a_local = c;
		}
	}

	// Round 2: wipe Parent A's winning positions, recompute for Parent B.
	// "Wipe" means: for all positions where A had the max smooth, set A's smooth to 0.
	// Then recompute wins excluding A.
	std::vector<int> wins2(candidates.size(), 0);
	for (size_t pos = win_start; pos < align_len; pos++) {
		int max_smooth = 0;
		for (size_t c = 0; c < candidates.size(); c++) {
			if (c == parent_a_local) {
				continue;
			}
			max_smooth = std::max(max_smooth, candidates[c].smoothed[pos]);
		}
		if (max_smooth == 0) {
			continue;
		}
		for (size_t c = 0; c < candidates.size(); c++) {
			if (c == parent_a_local) {
				continue;
			}
			if (candidates[c].smoothed[pos] == max_smooth) {
				wins2[c]++;
			}
		}
	}

	// Parent B = candidate with most wins in round 2 (excluding A).
	size_t parent_b_local = (parent_a_local == 0) ? 1 : 0;
	for (size_t c = 0; c < candidates.size(); c++) {
		if (c == parent_a_local) {
			continue;
		}
		if (wins2[c] > wins2[parent_b_local]) {
			parent_b_local = c;
		}
	}

	return ParentPair {candidates[parent_a_local].idx, candidates[parent_b_local].idx,
	                   std::move(candidates[parent_a_local].alignment),
	                   std::move(candidates[parent_b_local].alignment)};
}

// ============================================================================
// Phase 5: Breakpoint sweep + classification
// ============================================================================

BreakpointResult sweep_breakpoints(const std::vector<DiffType> &diffs, const UchimeParams &params) {
	// Count total diffs by type
	int sum_a = 0, sum_b = 0, sum_n_q = 0; // N and ? combined as "abstain/neutral"
	for (auto d : diffs) {
		switch (d) {
		case DiffType::MATCH_A:
			sum_a++;
			break;
		case DiffType::MATCH_B:
			sum_b++;
			break;
		case DiffType::NO_VOTE:
		case DiffType::ABSTAIN:
			sum_n_q++;
			break;
		case DiffType::IGNORE:
			break;
		}
	}

	if (sum_a + sum_b + sum_n_q == 0) {
		return {}; // No informative diffs
	}

	// Initialize: everything in the right segment.
	// In the normal (non-reversed) configuration:
	//   Left segment: A-diffs are yes votes, B-diffs are no votes
	//   Right segment: B-diffs are yes votes, A-diffs are no votes
	int left_y = 0, left_n = 0, left_a = 0;
	int right_y = sum_b, right_n = sum_a, right_a = sum_n_q;

	BreakpointResult best;

	for (size_t i = 0; i < diffs.size(); i++) {
		DiffType d = diffs[i];
		if (d == DiffType::IGNORE) {
			continue;
		}

		// Transfer this diff from right to left
		switch (d) {
		case DiffType::MATCH_A:
			left_y++;
			right_n--;
			break;
		case DiffType::MATCH_B:
			left_n++;
			right_y--;
			break;
		case DiffType::NO_VOTE:
		case DiffType::ABSTAIN:
			left_a++;
			right_a--;
			break;
		case DiffType::IGNORE:
			break;
		}

		// Normal configuration: left side has A-majority, right side has B-majority
		if (left_y > left_n && right_y > right_n) {
			double lh = static_cast<double>(left_y) / (params.xn * (left_n + params.dn) + left_a);
			double rh = static_cast<double>(right_y) / (params.xn * (right_n + params.dn) + right_a);
			double h = lh * rh;
			if (h > best.best_h) {
				best.best_h = h;
				best.best_pos = static_cast<int>(i);
				best.reversed = false;
				best.left_yes = left_y;
				best.left_no = left_n;
				best.left_abstain = left_a;
				best.right_yes = right_y;
				best.right_no = right_n;
				best.right_abstain = right_a;
			}
		}

		// Reversed configuration: left has B-majority, right has A-majority (parents swapped)
		if (left_n > left_y && right_n > right_y) {
			double lh = static_cast<double>(left_n) / (params.xn * (left_y + params.dn) + left_a);
			double rh = static_cast<double>(right_n) / (params.xn * (right_y + params.dn) + right_a);
			double h = lh * rh;
			if (h > best.best_h) {
				best.best_h = h;
				best.best_pos = static_cast<int>(i);
				best.reversed = true;
				best.left_yes = left_n;
				best.left_no = left_y;
				best.left_abstain = left_a;
				best.right_yes = right_n;
				best.right_no = right_y;
				best.right_abstain = right_a;
			}
		}
	}

	return best;
}

// ============================================================================
// Identity computation
// ============================================================================

double ChimeraDetector::compute_identity(const std::string &aligned_a, const std::string &aligned_b) {
	int matches = 0, total = 0;
	for (size_t i = 0; i < aligned_a.size(); i++) {
		if (is_gap(aligned_a[i]) || is_gap(aligned_b[i])) {
			continue;
		}
		total++;
		if (to_upper(aligned_a[i]) == to_upper(aligned_b[i])) {
			matches++;
		}
	}
	if (total == 0) {
		return 0.0;
	}
	return 100.0 * matches / total;
}

double ChimeraDetector::compute_model_identity(const StarAlignment &star, int breakpoint, bool reversed) {
	int matches = 0, total = 0;
	for (size_t i = 0; i < star.diffs.size(); i++) {
		char q = to_upper(star.query_row[i]);
		char a = to_upper(star.parent_a_row[i]);
		char b = to_upper(star.parent_b_row[i]);

		if (is_gap(q) || is_gap(a) || is_gap(b)) {
			continue;
		}
		if (is_ambiguous(q) || is_ambiguous(a) || is_ambiguous(b)) {
			continue;
		}

		total++;
		// Model: parent A on left side, parent B on right side (or reversed)
		char model_base;
		if (!reversed) {
			model_base = (static_cast<int>(i) <= breakpoint) ? a : b;
		} else {
			model_base = (static_cast<int>(i) <= breakpoint) ? b : a;
		}
		if (q == model_base) {
			matches++;
		}
	}
	if (total == 0) {
		return 0.0;
	}
	return 100.0 * matches / total;
}

// ============================================================================
// ChimeraDetector class
// ============================================================================

ChimeraDetector::ChimeraDetector(const UchimeParams &params) : params_(params) {
}

void ChimeraDetector::set_reference(const std::vector<std::string> &labels, const std::vector<std::string> &sequences) {
	ref_labels_ = labels;
	ref_sequences_ = sequences;
	kmer_index_ = KmerIndex();
	for (size_t i = 0; i < sequences.size(); i++) {
		kmer_index_.add_sequence(sequences[i]);
	}
}

void ChimeraDetector::add_to_reference(const std::string &label, const std::string &sequence, int64_t abundance) {
	ref_labels_.push_back(label);
	ref_sequences_.push_back(sequence);
	ref_abundances_.push_back(abundance);
	kmer_index_.add_sequence(sequence);
}

UchimeResult ChimeraDetector::detect(const std::string &query_label, const std::string &query_sequence,
                                     WFA2Aligner &aligner) const {
	UchimeResult result;
	result.query_label = query_label;

	// Step 1: Find candidate parents via k-mer index
	auto candidates = kmer_index_.find_candidates(query_sequence);
	if (candidates.size() < 2) {
		// No parents found
		result.parent_a_label = "*";
		result.parent_b_label = "*";
		result.closest_parent_label = "*";
		return result;
	}

	// Step 2: Select best two parents via smoothed identity
	auto parents = select_parents(query_sequence, candidates, ref_sequences_, aligner);
	if (!parents.has_value()) {
		result.parent_a_label = "*";
		result.parent_b_label = "*";
		result.closest_parent_label = "*";
		return result;
	}

	result.parent_a_label = ref_labels_[parents->parent_a_idx];
	result.parent_b_label = ref_labels_[parents->parent_b_idx];

	// Step 3: Build 3-way star alignment and classify diffs
	auto star = build_star_alignment(parents->align_a.query_aligned, parents->align_a.subject_aligned,
	                                 parents->align_b.query_aligned, parents->align_b.subject_aligned);
	classify_diffs(star);

	// Step 4: Sweep breakpoints
	auto bp = sweep_breakpoints(star.diffs, params_);

	result.score = bp.best_h;
	result.left_yes = bp.left_yes;
	result.left_no = bp.left_no;
	result.left_abstain = bp.left_abstain;
	result.right_yes = bp.right_yes;
	result.right_no = bp.right_no;
	result.right_abstain = bp.right_abstain;

	// Step 5: Compute identities
	result.id_query_a = compute_identity(parents->align_a.query_aligned, parents->align_a.subject_aligned);
	result.id_query_b = compute_identity(parents->align_b.query_aligned, parents->align_b.subject_aligned);

	// A-B identity: align parents to each other
	auto ab_align = aligner.align_full(ref_sequences_[parents->parent_a_idx], ref_sequences_[parents->parent_b_idx]);
	if (ab_align.has_value()) {
		result.id_a_b = compute_identity(ab_align->query_aligned, ab_align->subject_aligned);
	}

	// Model identity, closest parent, and divergence
	result.id_query_top = std::max(result.id_query_a, result.id_query_b);
	result.closest_parent_label =
	    (result.id_query_a >= result.id_query_b) ? result.parent_a_label : result.parent_b_label;

	if (bp.best_pos >= 0) {
		result.id_query_model = compute_model_identity(star, bp.best_pos, bp.reversed);
		result.divergence = result.id_query_model - result.id_query_top;
	}
	// When no valid breakpoint, id_query_model and divergence remain 0.0 (default)

	// Step 6: Classify
	int sum_left = bp.left_yes + bp.left_no + bp.left_abstain;
	int sum_right = bp.right_yes + bp.right_no + bp.right_abstain;

	if (bp.best_h >= params_.minh) {
		if (result.divergence >= params_.mindiv && sum_left >= params_.mindiffs && sum_right >= params_.mindiffs) {
			result.flag = "Y";
		} else {
			result.flag = "?";
		}
	} else {
		// Non-chimeric: match vsearch convention of reporting * for parents/identities
		result.flag = "N";
		result.score = 0.0;
		result.parent_a_label = "*";
		result.parent_b_label = "*";
		result.closest_parent_label = "*";
		result.id_query_model = 0.0;
		result.id_query_a = 0.0;
		result.id_query_b = 0.0;
		result.id_a_b = 0.0;
		result.id_query_top = 0.0;
		result.divergence = 0.0;
		result.left_yes = 0;
		result.left_no = 0;
		result.left_abstain = 0;
		result.right_yes = 0;
		result.right_no = 0;
		result.right_abstain = 0;
	}

	return result;
}

UchimeResult ChimeraDetector::detect_denovo(const std::string &query_label, const std::string &query_sequence,
                                            int64_t query_abundance, WFA2Aligner &aligner) const {
	UchimeResult result;
	result.query_label = query_label;

	// Step 1: Find candidate parents via k-mer index
	auto all_candidates = kmer_index_.find_candidates(query_sequence);

	// Step 1b: Filter by abundance skew — candidate must have abundance >= abskew * query_abundance
	double min_abundance = params_.abskew * static_cast<double>(query_abundance);
	std::vector<uint32_t> candidates;
	for (uint32_t idx : all_candidates) {
		if (idx < ref_abundances_.size() && static_cast<double>(ref_abundances_[idx]) >= min_abundance) {
			candidates.push_back(idx);
		}
	}

	if (candidates.size() < 2) {
		result.parent_a_label = "*";
		result.parent_b_label = "*";
		result.closest_parent_label = "*";
		return result;
	}

	// Steps 2-6: same as detect()
	auto parents = select_parents(query_sequence, candidates, ref_sequences_, aligner);
	if (!parents.has_value()) {
		result.parent_a_label = "*";
		result.parent_b_label = "*";
		result.closest_parent_label = "*";
		return result;
	}

	result.parent_a_label = ref_labels_[parents->parent_a_idx];
	result.parent_b_label = ref_labels_[parents->parent_b_idx];

	auto star = build_star_alignment(parents->align_a.query_aligned, parents->align_a.subject_aligned,
	                                 parents->align_b.query_aligned, parents->align_b.subject_aligned);
	classify_diffs(star);

	auto bp = sweep_breakpoints(star.diffs, params_);

	result.score = bp.best_h;
	result.left_yes = bp.left_yes;
	result.left_no = bp.left_no;
	result.left_abstain = bp.left_abstain;
	result.right_yes = bp.right_yes;
	result.right_no = bp.right_no;
	result.right_abstain = bp.right_abstain;

	result.id_query_a = compute_identity(parents->align_a.query_aligned, parents->align_a.subject_aligned);
	result.id_query_b = compute_identity(parents->align_b.query_aligned, parents->align_b.subject_aligned);

	auto ab_align = aligner.align_full(ref_sequences_[parents->parent_a_idx], ref_sequences_[parents->parent_b_idx]);
	if (ab_align.has_value()) {
		result.id_a_b = compute_identity(ab_align->query_aligned, ab_align->subject_aligned);
	}

	result.id_query_top = std::max(result.id_query_a, result.id_query_b);
	result.closest_parent_label =
	    (result.id_query_a >= result.id_query_b) ? result.parent_a_label : result.parent_b_label;

	if (bp.best_pos >= 0) {
		result.id_query_model = compute_model_identity(star, bp.best_pos, bp.reversed);
		result.divergence = result.id_query_model - result.id_query_top;
	}

	int sum_left = bp.left_yes + bp.left_no + bp.left_abstain;
	int sum_right = bp.right_yes + bp.right_no + bp.right_abstain;

	if (bp.best_h >= params_.minh) {
		if (result.divergence >= params_.mindiv && sum_left >= params_.mindiffs && sum_right >= params_.mindiffs) {
			result.flag = "Y";
		} else {
			result.flag = "?";
		}
	} else {
		result.flag = "N";
		result.score = 0.0;
		result.parent_a_label = "*";
		result.parent_b_label = "*";
		result.closest_parent_label = "*";
		result.id_query_model = 0.0;
		result.id_query_a = 0.0;
		result.id_query_b = 0.0;
		result.id_a_b = 0.0;
		result.id_query_top = 0.0;
		result.divergence = 0.0;
		result.left_yes = 0;
		result.left_no = 0;
		result.left_abstain = 0;
		result.right_yes = 0;
		result.right_no = 0;
		result.right_abstain = 0;
	}

	return result;
}

} // namespace miint
