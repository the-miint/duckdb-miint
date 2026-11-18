#!/usr/bin/env python
"""Generate the tax2tree placement parity fixture for AssignNames (src/TaxonomyDecorator.cpp).

AssignNames ports tax2tree's primary placement (tax2tree, BSD-3-Clause; Copyright (c) 2011,
The tax2tree project; https://github.com/biocore/tax2tree). This script runs tax2tree's own
placement stages, unchanged and in the order its `t2t decorate` command runs them, on small
random trees and taxonomies, and writes the names each internal node keeps as the expected
result:

  determine_rank_order, load_consensus_map, load_tree, collect_names_at_ranks_counts,
  decorate_ntips, decorate_name_relative_freqs, set_ranksafe, pick_names,
  name_node_score_fold

tax2tree is never a build- or test-time dependency: the outputs are committed, and
test/cpp/test_TaxonomyDecoratorParity.cpp reads them.

Each case is drawn from an RNG seeded from (SEED, the case's name), so adding a case never
changes another. The draws stay clear of the places where AssignNames deliberately differs
from tax2tree (the divergence ledger in the parity test lists them):
  - seven ranks, so no strain;
  - no name holds a bad-name token, and every lineage has a domain and a phylum, so tax2tree's
    input cleaning and TaxonomyStore's agree;
  - tip names are unique and hold no apostrophe;
  - a draw in which some name's candidates include two nodes with equal F1 from different
    counts, which tax2tree orders by floating-point rounding, is rejected and drawn again. So is
    one where tax2tree's floating-point scores order two candidates differently from their exact
    F1, which these tree sizes should never produce. provenance.txt counts both.

Outputs, under data/tax2tree/placement/:
  <case>.tree.tsv       node_index, parent_index (empty for the root), name (tips only).
                        Nodes are numbered in postorder, so siblings ascend left to right.
  <case>.taxonomy.tsv   feature_id, then the eight rank columns (empty for NULL)
  <case>.expected.tsv   node_index, rank, name: every name an internal node keeps after the
                        fold, by node, then rank
  cases.tsv             case, min_count
  provenance.txt        the tax2tree commit, the scikit-bio, numpy and Python versions, the
                        seed, the redraw counts and the command

Usage, from the repository root:
  conda run -n miint-beta-oracle python -I test/scripts/gen_tax2tree_placement_oracle.py \\
      --tax2tree PATH_TO_TAX2TREE_CHECKOUT
Running it again must leave `git diff data/tax2tree/placement` empty.
"""
import argparse
import hashlib
import os
import platform
import random
import subprocess
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(HERE, "..", "..", "data", "tax2tree", "placement"))
SEED = 20261117
N_CASES = 20
RANK_LETTERS = "dpcofgs"
COLUMNS = ["domain", "phylum", "class", "order", "family", "genus", "species", "strain"]


def case_rng(name, attempt):
    """An RNG for one draw of one case, independent of every other case."""
    digest = hashlib.sha256(f"{SEED}:{name}:{attempt}".encode()).digest()
    return random.Random(int.from_bytes(digest[:8], "little"))


def draw_shape(rng):
    """Children lists of a random rooted tree (node 0 the root) with 8 to 60 tips.

    Each node's parent is an earlier node, often the one just before it, so unary chains
    appear, and otherwise any earlier node, so multifurcations do. Each child list is then
    shuffled, so left-to-right order is not creation order.
    """
    while True:
        n = rng.randint(12, 90)
        children = [[] for _ in range(n)]
        for k in range(1, n):
            parent = k - 1 if rng.random() < 0.3 else rng.randrange(k)
            children[parent].append(k)
        tips = sum(1 for ch in children if not ch)
        if 8 <= tips <= 60:
            break
    for ch in children:
        rng.shuffle(ch)
    return children


def postorder_numbers(children):
    """Each shape node's postorder position, children visited left to right."""
    number = [None] * len(children)
    stack = [(0, False)]
    count = 0
    while stack:
        v, expanded = stack.pop()
        if expanded:
            number[v] = count
            count += 1
            continue
        stack.append((v, True))
        for child in reversed(children[v]):
            stack.append((child, False))
    return number


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
    """A tree and the lineage each tip carries, numbered in postorder.

    Lineages follow clades: a few internal nodes, the root always among them, each get a
    lineage, and a tip takes its nearest such ancestor's with probability 0.8, or a random one.
    About one tip in eight carries none.
    """
    children = draw_shape(rng)
    n = len(children)
    number = postorder_numbers(children)
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
        if children[v] or rng.random() < 0.12:
            continue
        a = v
        while a not in seeded:
            a = parent[a]
        lineage_of[number[v]] = seeded[a] if rng.random() < 0.8 else draw_lineage(rng, hierarchy)
    tree_rows = []
    for v in range(n):
        tree_rows.append((number[v], None if parent[v] is None else number[parent[v]], not children[v]))
    tree_rows.sort()
    return tree_rows, lineage_of


def build_skbio_tree(tree_rows):
    """An skbio TreeNode with each node's index in attribute `idx`, tips named T<index>.

    Internal nodes stay unnamed: load_tree reads an internal name as a support value.
    """
    from skbio import TreeNode

    nodes = {}
    for idx, _, is_tip in tree_rows:
        node = TreeNode(name=f"T{idx}" if is_tip else None)
        node.idx = idx
        nodes[idx] = node
    # Rows are in postorder, so appending each node to its parent in index order keeps siblings
    # left to right.
    root = None
    for idx, parent, _ in tree_rows:
        if parent is None:
            root = nodes[idx]
        else:
            nodes[parent].append(nodes[idx])
    assert [n.idx for n in root.postorder(include_self=True)] == [row[0] for row in tree_rows]
    return root


def consensus_lines(lineage_of):
    lines = []
    for idx in sorted(lineage_of):
        fields = [f"{RANK_LETTERS[r]}__{name or ''}" for r, name in enumerate(lineage_of[idx])]
        lines.append(f"T{idx}\t" + "; ".join(fields))
    return lines


def check_scores(nl, tree, counts):
    """After pick_names: whether every name's candidates order the same by tax2tree's float
    scores as by exact F1, with equal F1 only from equal counts.

    Returns (problem, genuine_ties): problem is None, "rational tie" or "float order";
    genuine_ties counts names whose best candidates tie on equal counts.
    """
    genuine_ties = 0
    for rank in range(len(RANK_LETTERS)):
        by_name = {}
        for node in tree.non_tips(include_self=True):
            name = node.RankNames[rank]
            if name is not None:
                by_name.setdefault(name, []).append(node)
        for name, nodes in by_name.items():
            total = counts[rank][name]
            scored = []
            for node in nodes:
                c = sum(1 for tip in node.tips() if tip.Consensus[rank] == name)
                m = node.NumTips
                exact = Fraction(2 * c, total + m)
                floating = nl.fmeasure(float(c) / m, float(c) / total)
                scored.append((exact, floating, c, m))
            for i in range(len(scored)):
                for j in range(i + 1, len(scored)):
                    ei, fi, ci, mi = scored[i]
                    ej, fj, cj, mj = scored[j]
                    if ei == ej and (ci, mi) != (cj, mj):
                        return "rational tie", 0
                    if (ei < ej) != (fi < fj) or (ei == ej) != (fi == fj):
                        return "float order", 0
            best = max(s[0] for s in scored)
            if sum(1 for s in scored if s[0] == best) > 1:
                genuine_ties += 1
    return None, genuine_ties


def run_tax2tree(nl, tree_rows, lineage_of, min_count):
    lines = consensus_lines(lineage_of)
    nl.determine_rank_order(lines[0].split("\t")[1])
    tipname_map = nl.load_consensus_map(lines, False)
    tree = nl.load_tree(build_skbio_tree(tree_rows), tipname_map)
    counts = nl.collect_names_at_ranks_counts(tree)
    nl.decorate_ntips(tree)
    nl.decorate_name_relative_freqs(tree, counts, min_count)
    nl.set_ranksafe(tree)
    nl.pick_names(tree)
    problem, genuine_ties = check_scores(nl, tree, counts)
    if problem is not None:
        return None, problem, 0
    nl.name_node_score_fold(tree)
    kept = []
    for node in tree.non_tips(include_self=True):
        for rank, name in enumerate(node.RankNames):
            if name is not None:
                kept.append((node.idx, rank, name.split("__", 1)[1]))
    return sorted(kept), None, genuine_ties


def write_tsv(path, header, rows):
    with open(path, "w", newline="\n") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join("" if x is None else str(x) for x in row) + "\n")


def git(args, cwd):
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tax2tree", required=True, help="a checkout of tax2tree")
    args = parser.parse_args()
    tax2tree = os.path.abspath(os.path.expanduser(args.tax2tree))
    sys.path.insert(0, tax2tree)
    import numpy
    import skbio
    import t2t.nlevel as nl

    os.makedirs(OUTDIR, exist_ok=True)
    redraws = {"rational tie": 0, "float order": 0}
    cases = []
    stats = []
    for i in range(N_CASES):
        name = f"case{i:02d}"
        attempt = 0
        while True:
            rng = case_rng(name, attempt)
            tree_rows, lineage_of = draw_case(rng)
            min_count = rng.choice([1, 2])
            kept, problem, genuine_ties = run_tax2tree(nl, tree_rows, lineage_of, min_count)
            if problem is None:
                break
            redraws[problem] += 1
            attempt += 1
        root = next(idx for idx, parent, _ in tree_rows if parent is None)
        write_tsv(
            os.path.join(OUTDIR, f"{name}.tree.tsv"),
            ["node_index", "parent_index", "name"],
            [(idx, parent, f"T{idx}" if is_tip else None) for idx, parent, is_tip in tree_rows],
        )
        write_tsv(
            os.path.join(OUTDIR, f"{name}.taxonomy.tsv"),
            ["feature_id", *COLUMNS],
            [(f"T{idx}", *lineage_of[idx], None) for idx in sorted(lineage_of)],
        )
        write_tsv(
            os.path.join(OUTDIR, f"{name}.expected.tsv"),
            ["node_index", "rank", "name"],
            [(idx, COLUMNS[rank], text) for idx, rank, text in kept],
        )
        cases.append((name, min_count))
        tips = sum(1 for row in tree_rows if row[2])
        below_root = sum(1 for idx, _, _ in kept if idx != root)
        stats.append(
            f"{name}: {len(tree_rows)} nodes, {tips} tips, {len(lineage_of)} with a lineage, min_count {min_count}, "
            f"{len(kept)} names kept ({below_root} below the root), {genuine_ties} names decided by a tie-break, "
            f"drawn {attempt + 1} time(s)"
        )
    write_tsv(os.path.join(OUTDIR, "cases.tsv"), ["case", "min_count"], cases)

    status = git(["status", "--porcelain", "--untracked-files=no"], tax2tree)
    with open(os.path.join(OUTDIR, "provenance.txt"), "w", newline="\n") as f:
        f.write("Generated by test/scripts/gen_tax2tree_placement_oracle.py.\n")
        f.write(f"tax2tree: commit {git(['rev-parse', 'HEAD'], tax2tree)}")
        f.write(", no tracked changes\n" if not status else ", WITH tracked changes\n")
        f.write(f"scikit-bio {skbio.__version__}, numpy {numpy.__version__}, Python {platform.python_version()}\n")
        f.write(f"seed {SEED}, {N_CASES} cases\n")
        f.write(
            f"redrawn: {redraws['rational tie']} for equal F1 from different counts, "
            f"{redraws['float order']} for float scores out of exact order\n"
        )
        f.write(
            "command: conda run -n miint-beta-oracle python -I test/scripts/gen_tax2tree_placement_oracle.py "
            "--tax2tree PATH_TO_TAX2TREE_CHECKOUT\n"
        )
        f.write("\n")
        for line in stats:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
