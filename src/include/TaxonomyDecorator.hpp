#pragma once

#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"

#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace miint {

// The scoring rule that decides which node of a tree a taxonomy name is assigned to:
// tax2tree's F1 measure. For a node, c is the number of informative tips below it carrying
// the name, m the number of informative tips below it, and total the number of informative
// tips in the tree carrying the name. Precision is c/m and recall c/total, so
// F1 = 2pr/(p+r) = 2c/(total+m), and nodes are compared by exact integer cross-multiplication.
// tax2tree compares F1 as doubles, so two nodes whose F1 is the same rational can come out
// unequal there: c=2, m=2 against c=3, m=5 at total=4 gives 0.6666666666666666 against
// 0.6666666666666665. Here they tie.
//
// Contract. Decoration finds each name's lowest eligible nodes by searching at most two root
// paths, and that search relies on both of these:
//   1. A node is eligible for a name exactly when its count c reaches a threshold that depends
//      only on (total, min_count) and is at least 1 and at least ceil(total/2). Counts never
//      decrease toward the root, so a name's eligible nodes are closed upward. A threshold of
//      at least 1 keeps every node without the name ineligible, and one of at least
//      ceil(total/2) means two disjoint eligible nodes hold at least total between them, so
//      there are at most two lowest ones, and two only when each holds exactly half of an
//      even total.
//   2. The score depends on (c, m, total) only.
// A rule that breaks (1), for example one under which an ancestor is ineligible while a
// descendant is eligible, breaks that search. Decoration would then need per-node count maps
// merged up the tree, small-to-large: O(n log n) map operations to build, against O(n) work to
// build this index and two binary searches per count. Nothing here supports that.
//
// The F-beta family generalizes F1 as (1+b^2)c/(b^2 total+m). For beta = 2 or 1/2 its cross
// products reach 5L^2 for L tips, against F1's 2L^2, which overflows uint64 at kMaxTips: it
// would need a cap of 2^30 tips, or wider products.
struct F1Score {
	// The most tips a tree may have. Every count is at most this, so the cross products in
	// Compare, at most kMaxTips * 2 * kMaxTips = 2^63, fit in uint64.
	static constexpr uint32_t kMaxTips = uint32_t(1) << 31;

	// The count a node needs to be eligible for a name: max(min_count, ceil(total/2)), so
	// recall is at least one half (2c >= total) and c is at least min_count. Throws
	// std::invalid_argument when min_count is 0, which would let a node with no informative
	// tips be eligible, with precision 0/0.
	static uint32_t Threshold(uint32_t total, uint32_t min_count);

	// Orders two nodes' scores for one name by F1: greater when (c1, m1) scores higher than
	// (c2, m2). Requires 1 <= c <= m <= kMaxTips and c <= total <= kMaxTips for each.
	static std::strong_ordering Compare(uint32_t c1, uint32_t m1, uint32_t c2, uint32_t m2, uint32_t total);

	// F1 as a double, for reporting only; nodes are ordered by Compare.
	static double ScoreForReport(uint32_t c, uint32_t m, uint32_t total);
};

// Counts of taxonomy names over a tree's subtrees: for a node, a rank and a name, how many
// informative tips below the node carry the name, in two binary searches.
//
// Each subtree is a contiguous range of preorder positions (TreeIndex), so the tips carrying
// a name, listed by preorder position, hold the ones below a node as one contiguous run.
// Memory is O(n) per-node arrays, one entry per (informative tip, named rank), and one offset
// per name in the store.
//
// A tip is informative when it carries a taxonomy row that kept a lineage after cleaning
// (TaxonomyStore::IsInformative). A tip with no row, or with a row cleaning emptied, counts
// in TipCount only.
//
// Accessors take nodes, ranks (kLineageRanks indexes) and NameIds in range for the tree and
// store the index was built from; kMissingName is not a NameId.
class CladeCountIndex {
public:
	// `row_of_node` has one entry per node of `tree`: the store row a tip carries, or nullopt.
	// Throws std::invalid_argument when its size is not the tree's, when an internal node
	// carries a row, when a row is out of range for `store`, or when the tree has more than
	// F1Score::kMaxTips tips.
	CladeCountIndex(TreeIndex tree, const TaxonomyStore &store, const std::vector<std::optional<size_t>> &row_of_node);

	const TreeIndex &Tree() const;

	// Informative tips below `node` (itself included) carrying `name` at `rank`.
	uint32_t CladeCount(uint32_t node, size_t rank, NameId name) const;
	// Informative tips below `node`: the denominator of precision. tax2tree calls it NumTips.
	uint32_t InformativeTipCount(uint32_t node) const;
	// All tips below `node`, informative or not: what tax2tree breaks ties between equally
	// scored nodes with.
	uint32_t TipCount(uint32_t node) const;
	// Informative tips in the tree carrying `name` at `rank`: the denominator of recall.
	uint32_t Total(size_t rank, NameId name) const;

	// The number of names at `rank` in the store the index was built from; NameIds at that
	// rank are 0..NameCount-1.
	size_t NameCount(size_t rank) const;
	// The preorder positions of the informative tips carrying `name` at `rank`, ascending.
	std::span<const uint32_t> Occurrences(size_t rank, NameId name) const;
	// The node at preorder position `p`.
	uint32_t NodeAtPreorder(uint32_t p) const;

private:
	TreeIndex tree_;
	// The node at each preorder position.
	std::vector<uint32_t> node_at_;
	// Over preorder positions: entry p counts the tips (informative tips) at positions less than p.
	std::vector<uint32_t> tip_prefix_;
	std::vector<uint32_t> informative_prefix_;
	// Per rank, the preorder positions of the informative tips carrying each name, ascending:
	// name k's are occurrences_[rank][offsets_[rank][k] .. offsets_[rank][k + 1]).
	std::array<std::vector<uint32_t>, kLineageRanks.size()> offsets_;
	std::array<std::vector<uint32_t>, kLineageRanks.size()> occurrences_;
};

// A name's lowest eligible nodes at `rank`: the nodes whose count reaches `threshold` while no
// child's does, in ascending preorder. There are none when even the root falls short of
// `threshold`, and two only when each holds exactly half of an even total. Throws
// std::invalid_argument for a threshold below 1 or below ceil(Total/2), which
// F1Score::Threshold never gives.
std::vector<uint32_t> LowestEligibleNodes(const CladeCountIndex &index, size_t rank, NameId name, uint32_t threshold);

// The min_count rank `rank` uses: `min_count`, except 1 at strain, the last rank. A strain on
// one tip can be held alone only by nodes spanning that tip, and under a min_count of 2 it would
// be eligible nowhere.
uint32_t MinCountAt(size_t rank, uint32_t min_count);

// Per node, bit r is set when rank r is safe there: exactly one name is eligible at rank r, each
// name's threshold being F1Score::Threshold(Total, MinCountAt(r, min_count)). Tips get 0, since
// no name is assigned to a tip. Throws std::invalid_argument when min_count is 0.
std::vector<uint8_t> RankSafeMasks(const CladeCountIndex &index, uint32_t min_count);

// The ranks a node takes names at, from its safe mask: the run of safe ranks that starts at its
// shallowest safe rank and stops at the first unsafe rank after it (tax2tree's pick_names).
uint8_t PickedRanks(uint8_t safe_mask);

// A name assigned to a node, and the counts its score came from: c of the node's m informative
// tips carry the name, of `total` in the tree.
struct NameScore {
	size_t rank;
	NameId name;
	uint32_t node;
	uint32_t c;
	uint32_t m;
	uint32_t total;

	bool operator==(const NameScore &) const = default;
};

// Why a name some informative tip carries was assigned to no node.
enum class NotPlacedReason {
	// Fewer tips carry it than its rank's min_count.
	kNeverEligible,
	// Eligible, but at every internal node where it is, so is another name at its rank.
	kNeverRankSafe,
	// Its rank is safe somewhere, but at each such node an unsafe shallower rank ended the
	// node's run of picked ranks first.
	kBlockedByContiguity,
};

struct NotPlacedName {
	size_t rank;
	NameId name;
	uint32_t total;
	NotPlacedReason reason;

	bool operator==(const NotPlacedName &) const = default;
};

// The names AssignNames gave each internal node, with the counts each score came from, and the
// names it gave no node, with why.
class NameAssignment {
public:
	// The names assigned to `node`, by ascending rank; empty for most nodes and every tip.
	std::span<const NameScore> NamesAt(uint32_t node) const;
	// The deepest of NamesAt(node), or nullptr when it is empty (tax2tree's
	// set_preliminary_name_and_rank).
	const NameScore *Preliminary(uint32_t node) const;
	// Every assigned name, by (rank, name).
	const std::vector<NameScore> &Scores() const;
	// Every name with an informative tip that was assigned nowhere, by (rank, name).
	const std::vector<NotPlacedName> &NotPlaced() const;

private:
	friend NameAssignment AssignNames(const CladeCountIndex &index, uint32_t min_count);
	// Keeps `scores` and `not_placed`, both by (rank, name), and groups the scores by node.
	NameAssignment(std::vector<NameScore> scores, std::vector<NotPlacedName> not_placed, uint32_t node_count);

	std::vector<NameScore> scores_;
	std::vector<NotPlacedName> not_placed_;
	// The scores again, grouped by node, each node's by ascending rank: node v's are
	// by_node_[node_offsets_[v] .. node_offsets_[v + 1]).
	std::vector<NameScore> by_node_;
	std::vector<uint32_t> node_offsets_;
};

// Assigns each name to the node that scores best for it (tax2tree's pick_names and
// name_node_score_fold). A node's candidates are its picked ranks (PickedRanks of its
// RankSafeMasks entry), each for the one name eligible there. A name's best candidate has the
// highest F1 (F1Score::Compare), then the fewest tips (TipCount), then is the latest in
// postorder, children visited in ascending index as TreeIndex visits them. Throws
// std::invalid_argument when min_count is 0.
NameAssignment AssignNames(const CladeCountIndex &index, uint32_t min_count);

// A name on a node at a rank (a kLineageRanks index).
struct RankName {
	size_t rank;
	NameId name;

	bool operator==(const RankName &) const = default;
};

// A node's label: its names by ascending rank, at most one per rank.
using NodeLabel = std::vector<RankName>;

// A node whose preliminary name is at a shallower rank than its nearest named ancestor's
// (tax2tree's levels < 0). tax2tree empties such a node's label, losing a name AssignNames
// placed, unless the name is a domain, which it never compares; here the node keeps its
// preliminary name, and this records that the ranks were inverted, at the domain too.
struct InvertedBackfill {
	uint32_t node;
	size_t rank;
	NameId name;
	uint32_t ancestor;
	size_t ancestor_rank;

	bool operator==(const InvertedBackfill &) const = default;
};

struct Backfill {
	// One label per node of the tree.
	std::vector<NodeLabel> labels;
	// By ascending node.
	std::vector<InvertedBackfill> inverted;
};

// Each node's label after tax2tree's set_preliminary_name_and_rank and backfill_names_gap. An
// internal node keeps its preliminary name (NameAssignment::Preliminary), preceded by the
// names the hierarchy (TaxonomyStore::Parent) gives it at the ranks between that name's rank
// and the preliminary rank of its nearest named ancestor, the nearest ancestor with a
// preliminary name. The node's other assigned names are dropped; where the fill reaches their
// ranks, the hierarchy's names stand in their place.
//   - Tips, and internal nodes with no preliminary name, get an empty label.
//   - With no named ancestor, the root included, the fill reaches rank 0. tax2tree skips such
//     a node and keeps its name alone, so the shallower names AssignNames placed on it are lost:
//     a root named down to the phylum leaves every lineage without a domain.
//   - With the ancestor's rank equal to the node's or one above it, there is no gap to fill.
//   - With the ancestor's rank deeper than the node's, the node keeps its preliminary name
//     alone and `inverted` records it. tax2tree empties the label.
// `assignment` must come from AssignNames over an index built on `tree` and `store`.
Backfill BackfillNames(const NameAssignment &assignment, const TaxonomyStore &store, const TreeIndex &tree);

// A name without its polyphyletic suffix, a trailing underscore and one or more capital letters:
// "Bacillus_A" and "Bacillus_AB" give "Bacillus". A name without one is returned whole, as is a
// name that is nothing but one ("_A"). Digits and lower-case letters are no suffix, so "Foo_1",
// as tax2tree's make_names_unique numbers a name it finds on several nodes, "Foo_A1" and
// "Candidatus_Foo" are returned whole. This is the suffix tax2tree's promotion looks for, but
// anchored at the end of the name; tax2tree's pattern is anchored only at the start, so it finds
// a suffix in "Foo_A1" and "Candidatus_Foo".
std::string_view PolyphyleticBase(std::string_view name);

// tax2tree's commonname_promotion over the labels BackfillNames gave a tree's nodes. In preorder,
// each internal node finds its nearest descendants with a non-empty label, searching below the
// empty ones, and takes their first names. While those agree, it appends the name to its own
// label and removes it from each of theirs, stopping after the round that leaves one of them with
// one name or none, so a descendant's last name can be taken. With one such descendant, agreement
// is trivial, so a name moves up to the topmost node whose only labeled descendant holds it.
//   - Two names agree when they are equal, or when one is the other's PolyphyleticBase, and the
//     suffixed one is promoted. tax2tree matches by substring, so Baci and Bacillus_A agree there.
//   - A name no deeper than every rank the node's path from the root holds, its own label
//     included, stops the node's promotion before it is appended: a promoted name says the node is
//     that taxon, which must fit below what the node's lineage already says it is, and a node with
//     an empty label says what its ancestors say. Such a name can only come from an inverted-rank
//     node, so an inverted name stays on the node it was placed on, and an inverted-rank node takes
//     no name at a rank its path holds above it. tax2tree appends the name, which moves an inverted
//     name up, or gives a label two names at one rank or a shallower name after deeper ones.
// Labels stay by ascending rank. `labels` must come from BackfillNames over `tree` and `store`.
std::vector<NodeLabel> PromoteCommonNames(std::vector<NodeLabel> labels, const TaxonomyStore &store,
                                          const TreeIndex &tree);

// A label as text: its names in rank order, each after its rank's letter and two underscores
// ("d__Bacteria", "s__Escherichia coli"), joined by "; ". An empty label gives "".
std::string FormatLabel(const NodeLabel &label, const TaxonomyStore &store);

// Decoration: AssignNames, then BackfillNames, then PromoteCommonNames.
struct DecorateResult {
	NameAssignment assignment;
	// One label per node, after promotion.
	std::vector<NodeLabel> labels;
	// The inversions BackfillNames recorded.
	std::vector<InvertedBackfill> inverted;
};

// Throws std::invalid_argument when min_count is 0.
DecorateResult Decorate(const CladeCountIndex &index, const TaxonomyStore &store, uint32_t min_count);

// A name of a label read from text, at a rank (a kLineageRanks index).
struct LabelPart {
	size_t rank;
	std::string name;

	bool operator==(const LabelPart &) const = default;
};

// A label read from text, as FormatLabel writes it or another tool did. Its components are
// separated by ';' and trimmed of whitespace (IsLineageSpace), and each is a rank letter
// (dpcofgst, or k for the domain), two underscores and a name. The parts come by ascending rank.
// An empty label has no parts. Throws std::invalid_argument, naming the label and the component,
// for an empty component, one not of that form, an unknown letter, an empty name, a name that
// starts with whitespace ("f__ A") or with a rank prefix of its own ("c__c__A"), and a rank given
// twice. Unlike ParseLineageString, which reads a whole lineage, the components may come in any
// order and skip ranks, as in a label another tool wrote (tax2tree can put a shallower name after
// a deeper one); either k or d names the domain; and an empty name is an error, not a missing
// rank.
std::vector<LabelPart> ParseNodeLabel(std::string_view label);

// A rank that a node's label holds and a label above it on its path already holds; the tips
// below get the ancestor's name at that rank.
struct RankRepeat {
	uint32_t node;
	size_t rank;
	// The shallowest node above `node` whose label holds `rank`.
	uint32_t ancestor;

	bool operator==(const RankRepeat &) const = default;
};

// Each tip's lineage from the labels on its path (tax2tree's pull_consensus_strings), given as
// the nodes whose labels supply each rank, so names are not copied once per tip.
struct EmittedLineages {
	static constexpr uint32_t kNone = UINT32_MAX;
	// The tips, by ascending node.
	std::vector<uint32_t> tips;
	// Per tip, per rank: the shallowest node on the path from the root to the tip whose label
	// holds the rank, since tax2tree, writing from the tip up, lets the shallower name win; kNone
	// where no label on the path holds it.
	std::vector<std::array<uint32_t, kLineageRanks.size()>> holders;
	// By ascending node, then rank.
	std::vector<RankRepeat> repeats;
};

// `labels` has one entry per node of `tree`, each with at most one part per rank, as
// ParseNodeLabel gives them. A tip's label is ignored, its name being its id. Throws
// std::invalid_argument when `labels` and `tree` differ in size, or when an internal node's label
// holds a rank twice or a rank past the last.
EmittedLineages EmitLineages(const TreeIndex &tree, const std::vector<std::vector<LabelPart>> &labels);

} // namespace miint
