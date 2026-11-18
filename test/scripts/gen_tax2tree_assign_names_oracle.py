#!/usr/bin/env python
"""Generate the tax2tree parity fixture for AssignNames (src/TaxonomyDecorator.cpp).

AssignNames ports how tax2tree gives a tree's internal nodes names: its pick_names and
name_node_score_fold (tax2tree, BSD-3-Clause; Copyright (c) 2011, The tax2tree project;
https://github.com/biocore/tax2tree). This script runs tax2tree's own stages, unchanged and in
the order its `t2t decorate` command runs them, on small random trees and taxonomies, and writes
the names each internal node keeps as the expected result:

  determine_rank_order, load_consensus_map, load_tree, collect_names_at_ranks_counts,
  decorate_ntips, decorate_name_relative_freqs, set_ranksafe, pick_names,
  name_node_score_fold

tax2tree is never a build- or test-time dependency: the outputs are committed, and
test/cpp/test_TaxonomyDecoratorParity.cpp reads them. That test's header lists where AssignNames
differs from tax2tree and how these draws keep clear of each difference. One of them, equal F1
from different counts, which tax2tree orders by floating-point rounding, is handled here: a draw
where some name's candidates include such a pair is rejected and drawn again.

Each case is drawn from an RNG seeded from (SEED, the case's name, the attempt), so adding a case
never changes another. Nodes are numbered at random, except that each node's children ascend
left to right, the order AssignNames reads siblings in; so neither node index nor its reverse
is postorder.

The script stops with an error, before writing anything, when the tax2tree checkout has tracked
changes, when tax2tree's floating-point scores order two of a name's candidates differently
from their exact F1, and when the fixture lacks any of: a name decided by fewer tips; one
decided by postorder; and among those, one whose node has a smaller index than another tied
node's, and one whose node has a larger index. Without the last two, a tie-break that compared
indices, in one direction or the other, would match.

Outputs, under data/tax2tree/assign_names/, each table gzip-compressed with no timestamp or
file name in its header, so the same tables give the same bytes:
  cases.tsv.gz          case, min_count
  trees.tsv.gz          case, node_index, parent_index (empty for the root), name (tips only);
                        each case's nodes in index order
  taxonomy.tsv.gz       case, feature_id, then the eight rank columns (empty for NULL)
  expected.tsv.gz       case, node_index, rank, name: every name an internal node keeps after
                        the fold, by node, then rank
  provenance.txt        the tax2tree commit, the scikit-bio, numpy and Python versions, the
                        seed, the redraws, how many names each tie-break decided and by which
                        index, the command

Usage, from the repository root:
  conda run -n miint-beta-oracle python -I test/scripts/gen_tax2tree_assign_names_oracle.py \\
      --tax2tree PATH_TO_TAX2TREE_CHECKOUT
Running it again must leave `git diff data/tax2tree/assign_names` empty.
"""
import argparse
import gzip
import hashlib
import os
import platform
import random
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(HERE, "..", "..", "data", "tax2tree", "assign_names"))
SEED = 20261117
N_CASES = 20
RANK_LETTERS = "dpcofgs"
COLUMNS = ["domain", "phylum", "class", "order", "family", "genus", "species", "strain"]


def case_rng(name, attempt):
    """An RNG for one draw of one case, independent of every other case."""
    digest = hashlib.sha256(f"{SEED}:{name}:{attempt}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "little"))


def draw_shape(rng):
    """Children lists of a random rooted tree (shape node 0 the root) with 8 to 60 tips, and
    the tips that must carry no lineage.

    Each node's parent is an earlier node, often the one just before it, so unary chains
    appear, and otherwise any earlier node, so multifurcations do. Then some internal nodes get
    one more tip, which carries no lineage: such a node has its other children's informative
    tips and one tip more, so where it ties with one of them on F1, fewer tips decides. Each
    child list is then shuffled, so left-to-right order is not creation order.
    """
    while True:
        n = rng.randint(12, 90)
        children = [[] for _ in range(n)]
        for k in range(1, n):
            parent = k - 1 if rng.random() < 0.3 else rng.randrange(k)
            children[parent].append(k)
        rowless = set()
        for v in range(n):
            if children[v] and rng.random() < 0.15:
                rowless.add(len(children))
                children[v].append(len(children))
                children.append([])
        tips = sum(1 for ch in children if not ch)
        if 8 <= tips <= 60:
            break
    for ch in children:
        rng.shuffle(ch)
    return children, rowless


def draw_numbering(rng, children):
    """Each shape node's index: a random permutation, re-sorted within each set of siblings so
    that they ascend left to right."""
    index = list(range(len(children)))
    rng.shuffle(index)
    for ch in children:
        for child, idx in zip(ch, sorted(index[c] for c in ch)):
            index[child] = idx
    return index


def postorder(children):
    """Shape nodes in postorder, children visited left to right."""
    order = []
    stack = [(0, False)]
    while stack:
        v, expanded = stack.pop()
        if expanded:
            order.append(v)
            continue
        stack.append((v, True))
        for child in reversed(children[v]):
            stack.append((child, False))
    return order


def draw_hierarchy(rng):
    """parent_name[r][k]: the index at rank r - 1 of name k at rank r's parent."""
    parent_name = []
    for r in range(len(RANK_LETTERS)):
        count = 1 + rng.randint(0, r + 1)
        parent_name.append([0 if r == 0 else rng.randrange(len(parent_name[r - 1])) for _ in range(count)])
    return parent_name


def draw_lineage(rng, parent_name):
    """Seven names, domain first, with None below a random depth of at least two."""
    depth = rng.randint(2, len(RANK_LETTERS))
    names = [None] * len(RANK_LETTERS)
    k = rng.randrange(len(parent_name[depth - 1]))
    for r in range(depth - 1, -1, -1):
        names[r] = f"N{r}_{k}"
        k = parent_name[r][k]
    return names


def draw_case(rng):
    """A tree, as (index, parent index, is_tip) rows by index, the lineage each tip index
    carries, and the indices in postorder.

    Lineages follow clades: a few internal nodes, the root always among them, each get a
    lineage, and a tip takes its nearest such ancestor's with probability 0.8, or a random one.
    About one tip in eight carries none, as do the tips draw_shape added for that purpose.
    """
    children, rowless = draw_shape(rng)
    n = len(children)
    index = draw_numbering(rng, children)
    parent = [None] * n
    for v, ch in enumerate(children):
        for child in ch:
            parent[child] = v
    hierarchy = draw_hierarchy(rng)
    seeded = {0: draw_lineage(rng, hierarchy)}
    for v in range(1, n):
        if children[v] and rng.random() < 0.25:
            seeded[v] = draw_lineage(rng, hierarchy)
    lineage_of = {}
    for v in range(n):
        if children[v] or v in rowless or rng.random() < 0.12:
            continue
        a = v
        while a not in seeded:
            a = parent[a]
        lineage_of[index[v]] = seeded[a] if rng.random() < 0.8 else draw_lineage(rng, hierarchy)
    tree_rows = sorted((index[v], None if parent[v] is None else index[parent[v]], not children[v]) for v in range(n))
    return tree_rows, lineage_of, [index[v] for v in postorder(children)]


def build_skbio_tree(tree_rows, expected_postorder):
    """An skbio TreeNode with each node's index in attribute `idx`, tips named T<index>.

    Internal nodes stay unnamed: load_tree would read a numeric internal name, such as an index,
    as a support value and drop it.
    """
    from skbio import TreeNode

    nodes = {}
    for idx, _, is_tip in tree_rows:
        node = TreeNode(name=f"T{idx}" if is_tip else None)
        node.idx = idx
        nodes[idx] = node
    # Siblings ascend left to right, so appending each node to its parent in index order keeps
    # them in that order.
    root = None
    for idx, parent, _ in tree_rows:
        if parent is None:
            root = nodes[idx]
        else:
            nodes[parent].append(nodes[idx])
    assert [n.idx for n in root.postorder(include_self=True)] == expected_postorder
    return root


def consensus_lines(lineage_of):
    lines = []
    for idx in sorted(lineage_of):
        fields = [f"{RANK_LETTERS[r]}__{name or ''}" for r, name in enumerate(lineage_of[idx])]
        lines.append(f"T{idx}\t" + "; ".join(fields))
    return lines


def candidates(tree):
    """After pick_names: rank -> name -> the nodes that picked it."""
    by_rank = {}
    for node in tree.non_tips(include_self=True):
        for rank, name in enumerate(node.RankNames):
            if name is not None:
                by_rank.setdefault(rank, {}).setdefault(name, []).append(node)
    return by_rank


def score_names(nl, tree, counts):
    """After pick_names: each name's candidates with their exact F1 and the float score the fold
    sorts them by, computed from the same ValidRelFreq and ConsensusRelFreq values the fold
    reads.

    Returns (rational_tie, decided, postorder_ties): rational_tie is True when two of a name's
    candidates have equal F1 from different counts; decided counts the names whose best
    candidates tie on F1, by the rule that separates them ("fewer tips" or "postorder");
    postorder_ties lists (rank, name, the tied nodes' indices) for each name postorder decides.
    Raises when the float scores order two candidates differently from their exact F1.
    """
    decided = {"fewer tips": 0, "postorder": 0}
    postorder_ties = []
    for rank, names in candidates(tree).items():
        for name, nodes in names.items():
            total = counts[rank][name]
            scored = []
            for node in nodes:
                c = sum(1 for tip in node.tips() if tip.Consensus[rank] == name)
                exact = Fraction(2 * c, total + node.NumTips)
                floating = nl.fmeasure(node.ValidRelFreq[rank][name], node.ConsensusRelFreq[rank][name])
                scored.append((exact, floating, c, node.NumTips, len(list(node.tips())), node.idx))
            for i in range(len(scored)):
                for j in range(i + 1, len(scored)):
                    ei, fi, ci, mi, _, _ = scored[i]
                    ej, fj, cj, mj, _, _ = scored[j]
                    if ei == ej and (ci, mi) != (cj, mj):
                        return True, None, None
                    if (ei < ej) != (fi < fj) or (ei == ej) != (fi == fj):
                        raise RuntimeError(f"rank {rank}, {name}: float scores order candidates unlike exact F1")
            best = max(s[0] for s in scored)
            tied_tips = [s[4] for s in scored if s[0] == best]
            if len(tied_tips) > 1:
                fewest = min(tied_tips)
                if tied_tips.count(fewest) == 1:
                    decided["fewer tips"] += 1
                else:
                    decided["postorder"] += 1
                    postorder_ties.append((rank, name, [s[5] for s in scored if s[0] == best and s[4] == fewest]))
    return False, decided, postorder_ties


def run_tax2tree(nl, tree_rows, lineage_of, expected_postorder, min_count):
    """The names each internal node keeps, or None for a draw to reject, and the tie-break
    counts, with the names postorder decided counted by whether a smaller or a larger index won
    (a winner between two tied nodes counts as both)."""
    lines = consensus_lines(lineage_of)
    nl.determine_rank_order(lines[0].split("\t")[1])
    tipname_map = nl.load_consensus_map(lines, False)
    tree = nl.load_tree(build_skbio_tree(tree_rows, expected_postorder), tipname_map)
    counts = nl.collect_names_at_ranks_counts(tree)
    nl.decorate_ntips(tree)
    nl.decorate_name_relative_freqs(tree, counts, min_count)
    nl.set_ranksafe(tree)
    nl.pick_names(tree)
    rational_tie, decided, postorder_ties = score_names(nl, tree, counts)
    if rational_tie:
        return None, None
    nl.name_node_score_fold(tree)
    kept = []
    holder = {}
    for node in tree.non_tips(include_self=True):
        for rank, name in enumerate(node.RankNames):
            if name is not None:
                kept.append((node.idx, rank, name.split("__", 1)[1]))
                holder[(rank, name)] = node.idx
    decided["postorder, won by a smaller index"] = 0
    decided["postorder, won by a larger index"] = 0
    for rank, name, tied in postorder_ties:
        winner = holder[(rank, name)]
        assert winner in tied, (rank, name, winner, tied)
        decided["postorder, won by a smaller index"] += winner < max(tied)
        decided["postorder, won by a larger index"] += winner > min(tied)
    return sorted(kept), decided


def write_tsv_gz(path, header, rows):
    text = "\t".join(header) + "\n"
    for row in rows:
        text += "\t".join("" if x is None else str(x) for x in row) + "\n"
    with open(path, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
        f.write(text.encode())


def git(args, cwd):
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tax2tree", required=True, help="a checkout of tax2tree")
    args = parser.parse_args()
    tax2tree = os.path.abspath(os.path.expanduser(args.tax2tree))
    if git(["status", "--porcelain", "--untracked-files=no"], tax2tree):
        sys.exit(f"{tax2tree} has tracked changes; generate from a clean checkout")
    sys.path.insert(0, tax2tree)
    import numpy
    import skbio
    import t2t.nlevel as nl

    redraws = 0
    decided_total = {}
    drawn = []
    stats = []
    for i in range(N_CASES):
        name = f"case{i:02d}"
        attempt = 0
        while True:
            rng = case_rng(name, attempt)
            tree_rows, lineage_of, expected_postorder = draw_case(rng)
            min_count = rng.choice([1, 2])
            kept, decided = run_tax2tree(nl, tree_rows, lineage_of, expected_postorder, min_count)
            if kept is not None:
                break
            redraws += 1
            attempt += 1
        for rule, count in decided.items():
            decided_total[rule] = decided_total.get(rule, 0) + count
        drawn.append((name, min_count, tree_rows, lineage_of, kept))
        root = next(idx for idx, parent, _ in tree_rows if parent is None)
        tips = sum(1 for row in tree_rows if row[2])
        below_root = sum(1 for idx, _, _ in kept if idx != root)
        stats.append(
            f"{name}: {len(tree_rows)} nodes, {tips} tips, {len(lineage_of)} with a lineage, min_count {min_count}, "
            f"{len(kept)} names kept ({below_root} below the root); tie-breaks decided {decided['fewer tips']} "
            f"by fewer tips, {decided['postorder']} by postorder; drawn {attempt + 1} time(s)"
        )
    for rule, count in decided_total.items():
        if count == 0:
            sys.exit(f"no name in the fixture is decided by {rule}; change SEED or N_CASES")

    os.makedirs(OUTDIR, exist_ok=True)
    write_tsv_gz(os.path.join(OUTDIR, "cases.tsv.gz"), ["case", "min_count"], [(name, mc) for name, mc, *_ in drawn])
    write_tsv_gz(
        os.path.join(OUTDIR, "trees.tsv.gz"),
        ["case", "node_index", "parent_index", "name"],
        [
            (name, idx, parent, f"T{idx}" if is_tip else None)
            for name, _, tree_rows, _, _ in drawn
            for idx, parent, is_tip in tree_rows
        ],
    )
    write_tsv_gz(
        os.path.join(OUTDIR, "taxonomy.tsv.gz"),
        ["case", "feature_id", *COLUMNS],
        [
            (name, f"T{idx}", *lineage_of[idx], None)
            for name, _, _, lineage_of, _ in drawn
            for idx in sorted(lineage_of)
        ],
    )
    write_tsv_gz(
        os.path.join(OUTDIR, "expected.tsv.gz"),
        ["case", "node_index", "rank", "name"],
        [(name, idx, COLUMNS[rank], text) for name, _, _, _, kept in drawn for idx, rank, text in kept],
    )

    with open(os.path.join(OUTDIR, "provenance.txt"), "w", newline="\n") as f:
        f.write("Generated by test/scripts/gen_tax2tree_assign_names_oracle.py.\n")
        f.write(f"tax2tree: commit {git(['rev-parse', 'HEAD'], tax2tree)}, no tracked changes\n")
        f.write(f"scikit-bio {skbio.__version__}, numpy {numpy.__version__}, Python {platform.python_version()}\n")
        f.write(f"seed {SEED}, {N_CASES} cases\n")
        f.write(f"redrawn for equal F1 from different counts: {redraws}\n")
        f.write(
            f"names decided by a tie-break: {decided_total['fewer tips']} by fewer tips, "
            f"{decided_total['postorder']} by postorder ({decided_total['postorder, won by a smaller index']} "
            f"won by a smaller index than a tied node's, {decided_total['postorder, won by a larger index']} "
            "by a larger)\n"
        )
        f.write(
            "command: conda run -n miint-beta-oracle python -I test/scripts/gen_tax2tree_assign_names_oracle.py "
            "--tax2tree PATH_TO_TAX2TREE_CHECKOUT\n"
        )
        f.write("\n")
        for line in stats:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
