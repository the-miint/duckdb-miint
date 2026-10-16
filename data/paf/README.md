# PAF fixtures

Inputs for the `read_paf` tests. The reads are cut from a reference at known coordinates
and their names say where each one must align, so a test can derive the expected PAF row
from how the read was built instead of copying it from minimap2's output.

## Files

| File | Contents |
|---|---|
| `k12_ref.fa` | One record, `k12_100000_120000`: E. coli K-12 W3110 (`data/sylph/tiny_refs/e.coli-K12.fasta.gz`) bases `[100000, 120000)`, 0-based half-open |
| `k12_reads.fa` | Seven reads built from that reference; see below |
| `k12_minimap2.paf` | minimap2's PAF for the reads against the reference |
| `k12_minimap2.paf.gz` | The same bytes, gzip-compressed with `gzip -n` |

## Reads

Coordinates are 0-based half-open on `k12_100000_120000`, the same convention PAF uses.

| Read | Construction | Alignment it must produce |
|---|---|---|
| `fwd_1000_3000` | reference `[1000, 3000)` | `+`, query `[0, 2000)`, target `[1000, 3000)` |
| `fwd_5000_9000` | reference `[5000, 9000)` | `+`, query `[0, 4000)`, target `[5000, 9000)` |
| `rev_10000_13000` | reverse complement of `[10000, 13000)` | `-`, query `[0, 3000)`, target `[10000, 13000)` |
| `rev_14000_18000` | reverse complement of `[14000, 18000)` | `-`, query `[0, 4000)`, target `[14000, 18000)` |
| `ovhfwd_16000_20000` | `[16000, 20000)` then 1500 bp of genome `[3500000, 3501500)` | `+`, query `[0, 4000)` of 5500, target `[16000, 20000)` |
| `ovhrev_17000_20000` | reverse complement of (`[17000, 20000)` then genome `[3510000, 3511500)`) | `-`, query `[1500, 4500)` of 4500, target `[17000, 20000)` |
| `unmapped` | genome `[3600000, 3603000)`, outside the reference | no hit |

The two overhang reads are what make a partial alignment visible. On the reverse one, PAF's
forward-strand query coordinates put the aligned block at the read's 3' end. A reader that
mirrored the coordinates would report `[0, 3000)` instead. Both overhang slices end at the
reference's end, so chance matches cannot extend an alignment into the foreign tail.

minimap2 writes the unmapped read (`--paf-no-hit`) with its query length (3000), strand `*`,
target `*`, and zeros for every other coordinate, length and count. Its only tag is `rl:i:0`.

## Regenerating

Generated with **minimap2 2.31-r1302** as `minimap2 -x map-hifi -c --paf-no-hit -t 1`. This is
the command-line minimap2, not the 2.30 that miint embeds for `align_minimap2`: the tests compare
the two, so the fixture comes from a separate build. The generator refuses any other version.

```bash
bash build.sh
bash test/scripts/generate_paf_fixture.sh
git diff --stat data/paf/        # expect no diff with minimap2 2.31
```
