#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Parameters for greedy clustering.
struct ClusterParams {
	double id = 0.0; // minimum identity threshold (0.0-1.0)
	int strand = 1;  // 1 = plus only, 2 = both strands
};

// Single clustering assignment result.
struct ClusterResult {
	std::string read_id;
	bool is_centroid = false;
	int cluster_id = 0;
	std::string centroid_id; // read_id of the centroid
	double identity = 0.0;   // percent identity (0-100); 100.0 if centroid
	std::string cigar;       // CIGAR alignment string (empty if centroid)
	bool cigar_truncated = false;
};

// Wrapper around vsearch's greedy clustering library API.
//
// Clusters sequences in the order provided (caller controls sort order).
// For cluster_fast behavior, sort by length descending before calling.
// For cluster_size behavior, sort by abundance descending.
//
// Single-threaded: clustering is inherently sequential because each new
// centroid must be indexed before the next sequence is processed.
class VsearchClusterWrapper {
public:
	explicit VsearchClusterWrapper(const ClusterParams &params = ClusterParams {});
	~VsearchClusterWrapper();

	VsearchClusterWrapper(const VsearchClusterWrapper &) = delete;
	VsearchClusterWrapper &operator=(const VsearchClusterWrapper &) = delete;
	VsearchClusterWrapper(VsearchClusterWrapper &&) noexcept;
	VsearchClusterWrapper &operator=(VsearchClusterWrapper &&) noexcept;

	// Load sequences into vsearch's global DB, apply DUST masking, and
	// prepare the k-mer index (empty — centroids indexed incrementally).
	// Acquires vsearch session mutex.
	void set_sequences(const std::vector<std::string> &labels, const std::vector<std::string> &sequences);

	// Cluster all loaded sequences sequentially.
	// Must be called AFTER set_sequences().
	// Returns one ClusterResult per input sequence.
	std::vector<ClusterResult> cluster_all();

private:
	ClusterParams params_;

	struct VsearchState;
	std::unique_ptr<VsearchState> state_;
};

} // namespace miint
