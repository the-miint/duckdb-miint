#pragma once
//
// The long-form feature-table cell produced by ReadFeatureTable.
//
// Split out of unifrac_support_biom.hpp so that consumers of the generic table
// readers (community_distances, cluster_kmeans, cluster_upgma) do not transitively
// include unifrac-binaries' api.hpp for the sake of a three-field struct. That
// include was the only thing tying those functions to the UniFrac feature group.
//
// The `miint::unifrac` namespace is retained so the many existing users need no
// change; the type itself is not UniFrac-specific.

#include <string>

namespace miint::unifrac {

struct CooRow {
	std::string sample_id;
	std::string feature_id;
	double count;
};

} // namespace miint::unifrac
