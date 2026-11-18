# tax2tree name-assignment parity fixture

Fixtures for `test/cpp/test_TaxonomyDecoratorParity.cpp`, which checks that `AssignNames`
(`src/TaxonomyDecorator.cpp`) gives every internal node the same names as tax2tree's own
`pick_names` and `name_node_score_fold`. That test's header lists where `AssignNames` differs
from tax2tree, and how these cases keep clear of each difference.

## Files

| File | Contents |
|---|---|
| `cases.tsv` | `case, min_count`: the 20 cases |
| `<case>.tree.tsv` | `node_index, parent_index, name`: the tree; the root's parent and internal nodes' names are empty. Siblings ascend left to right |
| `<case>.taxonomy.tsv` | `feature_id`, then the eight rank columns, empty for NULL: the lineage each named tip carries |
| `<case>.expected.tsv` | `node_index, rank, name`: every name an internal node keeps, by node, then rank |
| `provenance.txt` | the tax2tree commit, the scikit-bio, numpy and Python versions, the seed, how many draws were rejected, how many names each tie-break decided, and the command |

## Regenerating

From the repository root, with a clean tax2tree checkout:

```bash
conda run -n miint-beta-oracle python -I test/scripts/gen_tax2tree_assign_names_oracle.py \
    --tax2tree PATH_TO_TAX2TREE_CHECKOUT
```

The output is deterministic: running it again leaves `git diff data/tax2tree/assign_names`
empty. tax2tree is never a build- or test-time dependency.
