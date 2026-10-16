# GFA fixtures

## `spec_example.gfa`

A hand-written GFA 1.1 file containing every GFA 1.1 record type, each of which `read_gfa`
splits into positional fields and tags, plus the cases that make the split matter. The
GFA 1.2 jump record `J` is not in it; the C++ parser tests cover it.

| Line | Record | What it covers |
|---|---|---|
| 1 | `H` | a header that is all tags; `VN:Z:1.1` |
| 2 | `#` | a comment line, which is skipped but still counts toward line numbers |
| 3 | `S` with a sequence | two `i` tags |
| 4 | `S` with `*` | a missing sequence; an `A` tag whose value is `+` |
| 5 | `S` with a sequence | one tag of every other type: a negative `i`, an exponent `f`, a `Z` containing a space, `J`, `H`, `B:i` with a negative element, and `B:f` |
| 6 | `L` with a CIGAR overlap | a `Z` tag |
| 7 | `L` with `*` overlap | no tags at all, so the tag map is empty rather than NULL |
| 8 | `C` | six positional fields |
| 9 | `P` | three positional fields |
| 10 | `W` | six positional fields, including the 0-based `seq_start`/`seq_end` that stay as written |
| 11 | `x` | a record type `read_gfa` does not know, so every field stays positional; its last field, `with:no:tags`, looks like a tag and must not be parsed as one |

Positional field counts follow the GFA 1.1 specification
(https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).

Real assembler output is covered by `data/myloasm/ecoli_sim/final_contig_graph.gfa`.
Malformed inputs are written by the tests themselves, so each defect sits next to the
assertion that names it.
