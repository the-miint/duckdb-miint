# tax2tree decoration parity fixture

Fixtures for `test/cpp/test_TaxonomyDecoratorParity.cpp`, which checks that `Decorate`
(`src/TaxonomyDecorator.cpp`) labels every internal node as tax2tree's stages after the fold
do (`set_preliminary_name_and_rank`, `backfill_names_gap`, `commonname_promotion`), and that
`EmitLineages` gives every tip tax2tree's lineage (`pull_consensus_strings`). tax2tree's
correction stage is not run, as `Decorate` does not correct. That test's header numbers the
differences between miint and tax2tree; a case that hits one is left out, and `provenance.txt`
lists which, by those numbers.

## Files

Each table is gzip-compressed tab-separated text with a header line, and every row but those of
`cases.tsv.gz` starts with its case.

| File | Contents |
|---|---|
| `cases.tsv.gz` | `case, min_count`: 20 random cases (those of `data/tax2tree/assign_names/` that hit no difference, then `dec00` on) and `golden` |
| `trees.tsv.gz` | `case, node_index, parent_index, name`: each case's tree, its nodes in index order; the root's parent and internal nodes' names are empty. Siblings ascend left to right |
| `taxonomy.tsv.gz` | `case, feature_id`, then the eight rank columns, empty for NULL |
| `expected_labels.tsv.gz` | `case, node_index, label`: each labeled internal node's label, as tax2tree writes it |
| `expected_lineages.tsv.gz` | `case, feature_id`, then the eight rank columns: every tip's lineage, bare names, empty where tax2tree writes `x__` |
| `provenance.txt` | the tax2tree commit, the scikit-bio, numpy and Python versions, the cases left out and the differences each hits, the rejected draws, what the cases exercise, and the command |

`golden` is Greengenes-shaped: two domains, GTDB-style names and genome ids, a genus and its
`_A` variant that promotion merges, a genus whose larger clade leaves the smaller one to
promotion, and one genome with no taxonomy.

## Regenerating

From the repository root, with a clean tax2tree checkout:

```bash
conda run -n miint-beta-oracle python -I test/scripts/gen_tax2tree_assign_names_oracle.py \
    --tax2tree PATH_TO_TAX2TREE_CHECKOUT
```

The same script writes `data/tax2tree/assign_names/`. The output is deterministic for the
tax2tree commit and the scikit-bio, numpy and Python versions `provenance.txt` records: running
it again with those leaves `git diff data/tax2tree` empty. tax2tree is never a build- or test-time dependency.
