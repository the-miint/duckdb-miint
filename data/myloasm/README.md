# myloasm fixtures

## `ecoli_sim/`

Real myloasm output for simulated HiFi-like reads. Each read's name records where it came
from, so the myloasm reader tests can check their results against how the reads were built
rather than against myloasm's own description of its output.

Only the files miint reads are kept, at the paths myloasm writes them:

| File | What it is |
|---|---|
| `assembly_primary.fa` | Polished primary contigs |
| `alternate_assemblies/assembly_alternate.fa` | Contigs myloasm set aside as alternates of a primary; sequences included |
| `alternate_assemblies/duplicated_contigs.fa` | Contigs myloasm set aside as duplicates; sequences included |
| `final_contig_graph.gfa` | Contig graph: `S`, read-placement `a` lines, and `L` links, on the unpolished contigs |
| `3-mapping/map_to_unitigs.paf.gz` | The final mapping of every read to the unpolished contigs |
| `0-cleaning_and_unitigs/read_coverages.txt.gz` | Round-1 chimera breakpoints |
| `0-cleaning_and_unitigs/remap_temp/read_coverages.txt.gz` | Round-2 chimera breakpoints (empty here; see below) |
| `myloasm.log` | The run's log, renamed from myloasm's `myloasm_<timestamp>.log`, with the `SYSTEM HOST NAME` value replaced by `redacted` |

### The reads

The 1006 reads (11.6 Mb) are not committed. The generator rebuilds them from E. coli K-12
W3110 (`data/sylph/tiny_refs/e.coli-K12.fasta.gz`). Each name is
`s<index>_<kind>_<fields>_<f|r>`, with coordinates 0-based half-open on that genome. `r`
means the read is the reverse complement of the described sequence; odd indices are
reversed. The comment block at the top of `test/scripts/generate_myloasm_fixture.sh` defines
every kind:

| Kind | Reads | Built from | Purpose |
|---|---|---|---|
| `lina` | 500 | strain A of `[1000000, 1200000)`, about 30× | primary contig |
| `linb` | 333 | strain B: the same region with a substitution every 150 bp, about 20× | alternate contig |
| `circ` | 67 | a 40 kb circle at `[2000000, 2040000)`, reads wrapping past its end | circular contig |
| `circb` | 50 | a variant of that circle with a substitution every 1000 bp | duplicated contig |
| `low` | 13 | `[3000000, 3060000)`, about 2.5× | low-depth contig |
| `iso` | 3 | three non-overlapping reads from `[3300000, 3400000)` | reads myloasm assembles into nothing; absent from every output |
| `chim` | 20 | strain A then a circle read, joined | chimera splits at the junction |
| `ovh` | 20 | 9 kb of strain A then 3 kb of uncovered genome | a split that leaves a partial alignment, in both orientations |

### What myloasm produced

| Contig | File | Origin | Notes |
|---|---|---|---|
| `u491ctg` | primary | strain A | 196340 bp polished, `LN` 196466 |
| `u826ctg` | alternate | strain B | `secondary_to` `u491ctg`, ANI 99.18 |
| `u912ctg` | primary | low-coverage region | depth 2 at every tier |
| `u982ctg` | primary | the circle variant `circb` | `circular-yes`, self-link `L` |
| `u840ctg` | duplicated | the circle `circ` | `secondary_to` `u982ctg`, ANI 99.91; `circular-yes`, self-link `L` |

The PAF has 1085 rows over all five contigs. The round-1 breakpoint log has 40 lines, one
breakpoint each:
- 20 chimeras with `COND:1`
- 10 forward overhang reads with `COND:-3`
- 10 reverse overhang reads with `COND:-2`

PAF read names carry no suffix, `+split0` or `+split1`.

**Two myloasm behaviours do not occur in this fixture:**
- **A contig in the GFA but in no FASTA.** myloasm drops a contig this way when, after the
  GFA is written, more than 10% of it has no read coverage at depth ≤ 5, or when polishing
  skips it. Clean simulated reads did not trigger either case.
- **Round-2 breakpoints, and the two-level `+splitA+splitB` names they produce.**

### Format facts these files show

These are checked against myloasm 0.7.0's source.
- PAF tags are untyped (`shared_mini:N`, `diff_snp:N`, `shared_snp:N`).
- For `-` alignments, PAF query coordinates are on the reverse complement of the read piece.
- Columns 10 and 11 of the PAF are anchor spans, not residue matches and block length.
- `a`-line slices tile each contig exactly: offsets accumulate slice lengths up to `LN`.
- `a`-line read ranges are on the reverse complement for `-` reads.
- Contig coordinates in the GFA and PAF are on the unpolished contig. The FASTA sequences
  are polished, and circular contigs are trimmed: `u840ctg` is 51274 bp unpolished and
  40000 bp polished.

## Regenerating

Generated with **myloasm 0.7.0** (bioconda build `h7f49ad2_1`, linux-64), as
`myloasm reads.fa -o run1 --hifi -t 1`. On other platforms, which have their own 0.7.0 builds,
drop the build string from the command below; whether another build reproduces these files
exactly has not been checked.

```bash
conda create -n myloasm-070 -c conda-forge -c bioconda myloasm=0.7.0=h7f49ad2_1   # linux-64
bash build.sh
bash test/scripts/generate_myloasm_fixture.sh /path/to/fresh/work_dir
git diff --stat data/myloasm/   # expect a diff only in myloasm.log's timestamps
```

The generator refuses any other myloasm version. It runs myloasm twice and fails unless every
file above except the log is identical between the runs. It uses one thread because myloasm
writes the PAF and the breakpoint logs in thread order. Both full runs are left in the work
directory.
