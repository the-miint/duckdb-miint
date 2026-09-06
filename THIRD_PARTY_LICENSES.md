# Third-Party Licenses

This file contains licenses for third-party software referenced or used
by this project.

**Policy.** Full license text is inlined below for permissively licensed
dependencies (MIT / BSD / Apache 2.0) whose source is embedded in this
repository. For copyleft dependencies (LGPL-*, MPL) we instead identify
the SPDX expression and point at the license file that ships alongside
the upstream source — the source is the authoritative copy and
reproducing it here would create drift. Downstream distributors of a
binary that links any copyleft dependency listed below must ship the
corresponding full license text themselves (LGPL §6 / §10 obligations
are not satisfied by this file alone).

Dependencies whose source is **not** vendored in this repository — those
fetched at build time by vcpkg, and the OpenMP runtime supplied by the
system toolchain — are likewise identified by SPDX expression with a
pointer to where the authoritative license text ships (the vcpkg build
tree / the installed toolchain), rather than inlined here, for the same
anti-drift reason.

---

## HTSlib

SAM/BAM/CRAM file parsing. Used by `read_alignments`, `read_sequences_sam`,
`COPY FORMAT SAM`, `COPY FORMAT BAM`, and alignment analysis functions.

- Repository: https://github.com/samtools/htslib
- Version: 1.22.1
- License: MIT/Expat (main library), Modified BSD 3-Clause (cram/ subdirectory)

### Citation

Bonfield JK, Marshall J, Danecek P, Li H, Ohan V, Whitwham A, Keane T, Davies RM.
"HTSlib: C library for reading/writing high-throughput sequencing data."
*GigaScience*, 2021; 10(2):giab007.
doi: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007)

### MIT/Expat License

[Files outwith the cram/ subdirectory]

Copyright (C) 2012-2025 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

### Modified BSD 3-Clause License

[Files within the cram/ subdirectory]

Copyright (C) 2012-2025 Genome Research Ltd.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
   nor the names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## minimap2

Sequence-to-reference alignment. Used by `align_minimap2`,
`align_minimap2_sharded`, and `save_minimap2_index`. The bundled KSW2 alignment
library (`ksw2_*` sources, MIT-licensed) is the backend for `align_pairwise_ksw2_*`,
`align_pairwise_ksw2_dual_affine_*`, and `align_pairwise_ksw2_splice_*`.

- Repository: https://github.com/lh3/minimap2
- License: MIT

### Citation

Li H. "Minimap2: pairwise alignment for nucleotide sequences."
*Bioinformatics*, 2018; 34(18):3094-3100.
doi: [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

### MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
              2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## WFA2-lib

Pairwise sequence alignment using the Wavefront Alignment Algorithm (WFA).
Used by `align_pairwise_wfa2_score`, `align_pairwise_wfa2_cigar`, and `align_pairwise_wfa2_full`.

- Repository: https://github.com/smarco/WFA2-lib
- Version: v2.3.6
- License: MIT

### Citation

Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa.
"Fast gap-affine pairwise alignment using the wavefront algorithm."
*Bioinformatics*, 2021; 37(4):456-463.
doi: [10.1093/bioinformatics/btaa777](https://doi.org/10.1093/bioinformatics/btaa777)

### MIT License

Copyright (c) 2017 Santiago Marco-Sola

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## VSEARCH

Embedded as a static library (`ext/vsearch`, branch `v2.30.5-miint` from
`the-miint/vsearch` fork). Used by `detect_chimera_uchime`,
`detect_chimera_uchime_denovo`, `search_sequences`, `cluster_sequences`,
`merge_pairs`, and `mask_dust`. We use the BSD-2-Clause license (vsearch is
dual-licensed BSD-2-Clause / GPL-3.0).

- Upstream: https://github.com/torognes/vsearch
- Fork: https://github.com/the-miint/vsearch (branch `v2.30.5-miint`)
- License: BSD-2-Clause (chosen from dual-license with GPL-3.0)

### Citations

Rognes T, Flouri T, Nichols B, Quince C, Mahe F.
"VSEARCH: a versatile open source tool for metagenomics."
*PeerJ*, 2016; 4:e2584.
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

Edgar RC, Haas BJ, Clemente JC, Quince C, Knight R.
"UCHIME improves sensitivity and speed of chimera detection."
*Bioinformatics*, 2011; 27(16):2194-2200.
doi: [10.1093/bioinformatics/btr381](https://doi.org/10.1093/bioinformatics/btr381)

### BSD 2-Clause License

Copyright (C) 2014-2026, Torbjorn Rognes, Frederic Mahe and Tomas Flouri
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## kseq++

Modern C++ FASTA/FASTQ parser (header-only). Used by `read_fastx`,
`read_sequences_sff`, and sequence reading infrastructure.

- Repository: https://github.com/cartoonist/kseqpp
- License: MIT

### MIT License

Copyright (c) 2018, Ali Ghaffaari
Max-Planck-Institut fuer Informatik

This source code is released under the terms of the MIT License.

---

## LBFGS++

Header-only L-BFGS quasi-Newton optimizer (Eigen-backed). Drives the
maximum-a-posteriori fit in `mmvec`. Vendored as tracked file copies at
`ext/LBFGSpp/` (the same treatment as kseq++ and concurrentqueue: header-only,
no build step). Only the unbounded solver is vendored; the bounded `LBFGSB`
solver is omitted as unused. Every vendored file is unmodified — see
`ext/LBFGSpp/PROVENANCE.md` for per-file checksums.

- Repository: https://github.com/yixuan/LBFGSpp
- Version: v0.4.0, tag commit `c524a407fb85b74807f53de5a3ca2ddbcc164e54`
- License: MIT (also vendored verbatim at `ext/LBFGSpp/LICENSE.md`)

Portions derive from Jorge Nocedal's original Fortran L-BFGS and from
libLBFGS (Naoaki Okazaki); all copyright holders are listed below.

### MIT License

Copyright (c) 1990 Jorge Nocedal

Copyright (c) 2007-2010 Naoaki Okazaki

Copyright (c) 2016-2023 Yixuan Qiu

Copyright (c) 2018-2023 Dirk Toewe

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

## Biopython

The SFF binary parser (`src/SFFReader.cpp`) was developed with reference to
Biopython's `Bio.SeqIO.SffIO` module for understanding the SFF file format.

- Repository: https://github.com/biopython/biopython
- License: Biopython License Agreement (also available under BSD 3-Clause)

### Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.

### BSD 3-Clause License

Copyright (c) 1999-2025, The Biopython Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## SciPy

The Procrustes core (`src/procrustes_core.cpp`, `src/include/procrustes_core.hpp`)
is ported/adapted from SciPy: `scipy.spatial.procrustes` and
`scipy.linalg.orthogonal_procrustes`. As a source-derived work it inlines the
full license text below.

- Repository: https://github.com/scipy/scipy
- License: BSD 3-Clause

### BSD 3-Clause License

Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## Pyteomics

Monoisotopic element masses in `src/formula_function.cpp` are derived from the
`_nist_mass` data table in pyteomics (`pyteomics/auxiliary/constants.py`).

- Repository: https://github.com/levitsky/pyteomics
- Version: commit 3f1fd4afb51a5033222851666bef585c9253cd68
- License: Apache License 2.0

### Citations

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013)
"Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software
Prototyping in Proteomics", Journal of The American Society for Mass Spectrometry,
24(2), 301-304. doi: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018)
"Pyteomics 4.0: five years of development of a Python proteomics framework",
Journal of Proteome Research.
doi: 10.1021/acs.jproteome.8b00717

### Apache License 2.0

Copyright (c) Lev Levitsky and Joshua Klein

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---

## MAFFT

Multiple sequence alignment. Used by `align_mafft` (PartTree algorithm).
Embedded as a statically linked C library built from a modified fork
with a library-callable API (`splittbfast_library()`).

- Repository: https://github.com/GSLBiotech/mafft (upstream), https://mafft.cbrc.jp/alignment/software/ (official)
- License: BSD 3-Clause

### Citations

Katoh K, Toh H.
"PartTree: an algorithm to build an approximate tree from a large number of
unaligned sequences."
*Bioinformatics*, 2007; 23(3):372-374.
doi: [10.1093/bioinformatics/btl592](https://doi.org/10.1093/bioinformatics/btl592)

Katoh K, Standley DM.
"MAFFT Multiple Sequence Alignment Software Version 7: Improvements in
Performance and Usability."
*Molecular Biology and Evolution*, 2013; 30(4):772-780.
doi: [10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

### BSD 3-Clause License

MAFFT: multiple sequence alignment program
Copyright (c) 2009 Kazutaka Katoh

Redistribution and use in source and binary forms,
with or without modification, are permitted provided
that the following conditions are met:

Redistributions of source code must retain the
above copyright notice, this list of conditions
and the following disclaimer.  Redistributions in
binary form must reproduce the above copyright
notice, this list of conditions and the following
disclaimer in the documentation and/or other
materials provided with the distribution.

The name of the author may not be used to endorse
or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.

---

## krepp

k-mer LSH index and maximum pseudo-likelihood phylogenetic placement. Used
by `place_krepp`. Embedded as a statically linked C++ library built from
11 of the upstream `src/*.cpp` files; `src/krepp.cpp` (the CLI entry point)
is excluded. miint reads an index built by krepp's own command-line tool.
Pinned at v0.9.1.

- Repository: https://github.com/bo1929/krepp
- License: MIT

### MIT License

Copyright (c) 2026 Ali Osman Berk Şapcı

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### Code vendored inside krepp

krepp's `src/` carries several third-party files that are compiled into
miint along with the rest of the library. Upstream confirmed their
provenance (Ali Osman Berk Şapcı, 2026-09-05); the files themselves carry
author attribution but, apart from MurmurHash3, no inline license text, so
it is recorded here.

| File(s) | Origin | License |
|---|---|---|
| `src/hyperloglog.hpp` | [cpp-HyperLogLog](https://github.com/hideo55/cpp-HyperLogLog), Hideaki Ohno | MIT |
| `src/kdq.h`, `src/kalloc.h`, `src/ketopt.h`, `src/kseq.h`, `src/kvec.h` | [klib](https://github.com/attractivechaos/klib), Attractive Chaos / Heng Li | MIT |
| `src/MurmurHash3.{cpp,hpp}` | MurmurHash3, Austin Appleby | Public domain — the file states "placed in the public domain. The author hereby disclaims copyright to this source code." (upstream described it as MIT; the vendored file's own header is the record here) |

---

## abPOA

Adaptive banded partial order alignment for multiple sequence alignment
and consensus generation. Used by `align_abpoa` and `consensus_abpoa`.
Embedded as a statically linked C library via `add_subdirectory` from a
fork with symbol namespacing, CMake SIMD dispatch, and subdirectory
friendliness.

- Repository: https://github.com/yangao07/abPOA (upstream)
- License: MIT

### Citation

Gao Y, Liu Y, Ma Y, Liu B, Wang Y, Xing Y.
"abPOA: an SIMD-based C library for fast partial order alignment using
adaptive band."
*Bioinformatics*, 2021; 37(15):2209-2211.
doi: [10.1093/bioinformatics/btaa963](https://doi.org/10.1093/bioinformatics/btaa963)

### MIT License

Copyright (c) 2020 Yan Gao

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## SortMeRNA

rRNA read filtering and alignment. Used by `align_sortmerna` and
`align_sortmerna_rrna`. Embedded as a statically linked C++ library
from a fork with a streaming C API (`smr_run_seqs_with_index`).
Full source is available under `ext/sortmerna/` in this repository,
which satisfies the LGPL-3.0 source-availability requirement for
distributed binaries that link against it.

- Repository (upstream): https://github.com/biocore/sortmerna
- Repository (fork used here): https://github.com/the-miint/sortmerna — branch `v4.4.0-miint`, pinned at the submodule SHA recorded in `.gitmodules` / `git submodule status ext/sortmerna`
- Version: 4.4.0
- License: LGPL-3.0-or-later (SPDX)

### Citation

Kopylova E, Noé L, Touzet H.
"SortMeRNA: Fast and accurate filtering of ribosomal RNAs in
metatranscriptomic data."
*Bioinformatics*, 2012; 28(24):3211-3217.
doi: [10.1093/bioinformatics/bts611](https://doi.org/10.1093/bioinformatics/bts611)

### GNU Lesser General Public License v3.0

SortMeRNA is dual-licensed under LGPL-3.0 and GPL-3.0; this project
uses it under LGPL-3.0-or-later. The full LGPL-3.0 text is in
`ext/sortmerna/LICENSE.LESSER.txt`. The required GPL-3.0 text
(referenced by LGPL-3.0 §0) is in `ext/sortmerna/LICENSE.txt`.

---

## RocksDB

Embedded key-value store used by SortMeRNA for index persistence.
Sourced via vcpkg and linked statically into `libsortmerna_bundle.a`.

- Repository: https://github.com/facebook/rocksdb
- License: `GPL-2.0-only OR Apache-2.0` (SPDX expression, dual); used here under Apache 2.0. See RocksDB's `LICENSE.Apache` for full text.

### Apache License 2.0 (Summary)

Copyright (c) Meta Platforms, Inc. and affiliates.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
implied. See the License for the specific language governing
permissions and limitations under the License.

---

## cmph

Minimal perfect hash functions. Transitive dependency of SortMeRNA,
bundled into `libsortmerna_bundle.a` via the post-build `ar`-append
step in `cmake/bundle_sortmerna.cmake`.

- Repository: http://cmph.sourceforge.net/
- Vendored location in this repo: `ext/sortmerna/3rdparty/cmph/` (submodule; same pin as SortMeRNA)
- License: `LGPL-2.0-only OR MPL-1.1` (SPDX expression, dual)

The cmph sources vendored inside SortMeRNA carry no per-file license
headers. The upstream SourceForge project declares the dual license
above. Downstream distributors of a binary linking cmph should include
the full LGPL-2.0 and/or MPL-1.1 text; inlining them here would
duplicate content that the upstream project does not itself ship
alongside the sources.

---

## alp (ascending ladder algorithm for p-values)

NCBI statistical library used by SortMeRNA for alignment p-value
computation. Transitive dependency bundled into
`libsortmerna_bundle.a`.

- Origin: NCBI BLAST+ toolkit
- License: Public domain (U.S. government work per 17 U.S.C. §105)

### NCBI Public Domain Notice

This software/database is a "United States Government Work" under
the terms of the United States Copyright Act. It was written as part
of the author's official duties as a United States Government
employee and thus cannot be copyrighted. This software/database is
freely available to the public for use. The National Library of
Medicine and the U.S. Government have not placed any restriction on
its use or reproduction.

Although all reasonable efforts have been taken to ensure the
accuracy and reliability of the software and data, the NLM and the
U.S. Government do not and cannot warrant the performance or
results that may be obtained by using this software or data. The
NLM and the U.S. Government disclaim all warranties, express or
implied, including warranties of performance, merchantability or
fitness for any particular purpose.

---

## fastp

Read-level QC algorithms (adapter trimming, paired-end overlap analysis,
polyG/polyX tail trimming, sliding-window quality trimming, per-read
quality/length/N-base filtering) used by the `qc_*` scalar function family
(`trim_adapters`, `trim_adapters_pe`, `trim_polyg`, `trim_polyx`,
`trim_quality_3p`, `trim_quality_5p`, `trim_quality_sliding`, `filter_read`).
The ported routines include fastp's `OverlapAnalysis` (the paired-end
overlap-based adapter trimming behind `trim_adapters_pe`). fastp is **not** a
build dependency; specific algorithms have been ported into
`src/qc_algorithms.cpp` with per-file SPDX attribution. Adapter auto-detection
is **not** ported — adapters are user-supplied.

- Repository: https://github.com/OpenGene/fastp
- License: MIT

### Citation

Chen S, Zhou Y, Chen Y, Gu J.
"fastp: an ultra-fast all-in-one FASTQ preprocessor."
*Bioinformatics*, 2018; 34(17):i884–i890.
doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

### MIT License

Copyright (c) 2016 OpenGene - Open Source Genetics Toolbox

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## concurrentqueue

Lock-free concurrent queue used by SortMeRNA's producer/consumer
pipeline. Vendored as a header-only include at
`ext/concurrentqueue/concurrentqueue.h` (pinned to a specific commit
for build reproducibility; not carried as a submodule to keep CI
checkout simple).

- Repository: https://github.com/cameron314/concurrentqueue
- License: Simplified BSD / Boost Software License 1.0 (dual);
  used under the Simplified BSD terms.

### Simplified BSD License

Copyright (c) 2013-2016, Cameron Desrochers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the
  distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

---

## unifrac-binaries (libssu)

UniFrac phylogenetic distances and Faith's PD engine. Embedded as a
statically linked archive (`libssu_inmem.a` native / `libssu_wasm.a`
Emscripten). Used by `unifrac_distances`, `unifrac_pcoa`,
`unifrac_permanova`, and `unifrac_faith_pd`. Reported as the `unifrac`
row in `miint_versions()`.

- Repository: https://github.com/biocore/unifrac-binaries (branch `main`, pinned at the submodule SHA in `.gitmodules` / `git submodule status ext/unifrac-binaries`)
- License: BSD 3-Clause

### BSD 3-Clause License

BSD 3-Clause License

Copyright (c) 2022, biocore
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## scikit-bio-binaries (libskbb)

Principal Coordinates Analysis (randomized FSVD) and PERMANOVA pseudo-F.
Embedded as a statically linked archive (`libskbb_inmem.a` native /
`libskbb_wasm.a` Emscripten); backs the PCoA/PERMANOVA paths of
`unifrac_pcoa` and `unifrac_permanova` (and the metric-agnostic `pcoa`
/ `permanova` functions). Reported as the `scikit-bio-binaries` row in
`miint_versions()`.

- Repository: https://github.com/scikit-bio/scikit-bio-binaries (branch `main`, pinned at the submodule SHA in `.gitmodules` / `git submodule status ext/scikit-bio-binaries`)
- License: BSD 3-Clause

### BSD 3-Clause License

Copyright (c) 2013--, scikit-bio development team.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the names scikit-bio, skbio, or biocore nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## rype

Rust implementations of k-mer classification, minimizer extraction, and
log-ratio operations, consumed via the Arrow C Data Interface. Its
symbols ship inside the merged `libmiint_rust_glue.a` static archive
(built from the `ext/miint-rust-glue/` umbrella crate) and are linked
into the extension. Used by `rype_classify`, `rype_extract_minimizer_set`,
`rype_extract_strand_minimizers`, `rype_index_create`, and
`rype_log_ratio`. Reported as the `rype` row in `miint_versions()`.

- Repository: https://github.com/the-miint/rype (pinned at the submodule SHA in `.gitmodules` / `git submodule status ext/rype`)
- License: BSD 3-Clause (Modified BSD)

### Modified BSD License

Copyright (c) 2024-, The RYpe Development Team <damcdonald@ucsd.edu>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the RYpe development team nor the names of its
      contributors may be used to endorse or promote products derived from this
      software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE RYpe DEVELOPMENT TEAM BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

---

## sylph

FracMinHash sketch-based relative-abundance profiling of microbial
communities, consumed via the Arrow C Data Interface. Its symbols ship
inside the same merged `libmiint_rust_glue.a` archive as rype. Used by
the `sylph_profile` and `sylph_index_create` table functions. Reported
as the `sylph` row in `miint_versions()`.

- Repository: https://github.com/the-miint/sylph (branch `v0.9.0-miint`, pinned at the submodule SHA in `.gitmodules` / `git submodule status ext/sylph`)
- License: `MIT OR Apache-2.0` (SPDX expression, dual). The MIT text ships at `ext/sylph/LICENSE` and is reproduced below.

### MIT License

MIT License

Copyright (c) 2023 Jim Shaw

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## microtar

Minimal TAR archive reader. Compiled as C directly into the extension
(`third_party/microtar/microtar.c`), no feature flag. Used by
`src/taxdump_archive.cpp` to read NCBI taxdump `.tar` archives from an
in-memory buffer (read side only).

- Repository: https://github.com/rxi/microtar
- License: MIT

### MIT License

Copyright (c) 2017 rxi

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## nanoarrow / nanoarrow_ipc / flatcc

Arrow IPC stream byte-serialization used by the GPL-boundary transport
(`src/gpl_boundary/arrow_ipc.cpp`). DuckDB exposes only the Arrow C Data
Interface, not IPC byte serialization, so the official `nanoarrow_ipc`
(plus the `nanoarrow` C runtime and the `flatcc` FlatBuffers runtime it
depends on) is vendored under `third_party/nanoarrow_ipc/` and compiled
as the `miint_nanoarrow_ipc` object library linked into the extension.
Only built when `MIINT_ENABLE_GPL_BOUNDARY` is on (auto-off on
Emscripten and Windows). Symbols are namespaced (`NANOARROW_NAMESPACE=miint`)
to avoid colliding with DuckDB's bundled nanoarrow.

- Repository: https://github.com/apache/arrow-nanoarrow (pinned tag `apache-arrow-nanoarrow-0.8.0`)
- License: Apache License 2.0 (both nanoarrow/nanoarrow_ipc and the bundled flatcc runtime)

### Apache License 2.0 (Summary)

Full text ships in-tree at `third_party/nanoarrow_ipc/LICENSE` (nanoarrow)
and `third_party/nanoarrow_ipc/LICENSE.flatcc` (flatcc runtime).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---

## OpenMP runtime (libgomp / libomp)

The UniFrac archives (`libssu_inmem.a` / `libskbb_inmem.a`) are compiled
with `-fopenmp`, so the final extension must resolve their OpenMP runtime
symbols. This runtime is **not** vendored in this repository — it is
supplied by the platform toolchain and linked in via the `miint_openmp`
INTERFACE target (see `CMakeLists.txt`). Only present when
`MIINT_ENABLE_UNIFRAC` is on and not building for Emscripten.

- **Linux:** GCC's `libgomp`, `find_package(OpenMP)` → `OpenMP::OpenMP_CXX`.
  License: `GPL-3.0-or-later WITH GCC-exception-3.1` (SPDX) — GPLv3 with the
  GCC Runtime Library Exception. The exception explicitly permits linking
  the runtime into a program **regardless of that program's own license**,
  so it imposes no copyleft obligation on this BSD-licensed extension. Full
  text ships with the installed GCC toolchain (e.g.
  `/usr/share/gcc/COPYING.RUNTIME` / `COPYING3`).
- **macOS:** LLVM's `libomp`, linked statically (`libomp.a`) from Homebrew's
  `libomp` keg. License: `Apache-2.0 WITH LLVM-exception` (SPDX). Full text
  ships with the installed keg (`$(brew --prefix libomp)/LICENSE.TXT`).

Neither runtime causes copyleft to propagate to the extension; both carry
explicit linking exceptions that authorize inclusion in non-GPL binaries.

---

## System libraries (via vcpkg)

The following libraries are fetched at configure time by vcpkg (see
`vcpkg.json`) and compiled into or linked against the extension. Their source is
not vendored in this repository; vcpkg places each library's authoritative
license text under its build tree (e.g.
`vcpkg/buildtrees/<port>/src/.../COPYING`) and installed share directory
(`vcpkg_installed/<triplet>/share/<port>/copyright`). All are permissive except
Eigen, which is weak-copyleft (MPL-2.0) and is therefore restricted at compile
time to its pure-MPL2 subset — see its row below.

| Library | Purpose | SPDX / License | Notes |
|---|---|---|---|
| zlib | DEFLATE (de)compression for gzip I/O | `Zlib` | |
| zstd | Zstandard (de)compression | `BSD-3-Clause` (dual with `GPL-2.0-only`; used under BSD-3-Clause) | |
| libdeflate | Fast DEFLATE for bgzf gzip COPY | `MIT` | Linux only (`vcpkg.json` `platform: linux`) |
| expat | XML parsing (mzML / ENA payloads) | `MIT` | |
| curl | HTTP(S) transfers (INSDC/ENA paths) | `curl` (MIT/X-derivative) | Linked only when `MIINT_ENABLE_CURL` is on (auto-disabled on macOS to avoid an MD5/SHA1 symbol clash with vsearch) |
| OpenSSL | TLS backend | `Apache-2.0` (OpenSSL 3.x) | Pulled in transitively by curl; present only where curl is linked |
| HDF5 | HDF5 container I/O (fast5, etc.) | `BSD-3-Clause`-style (The HDF Group license) | Not built for Emscripten (`vcpkg.json` `platform: !emscripten`) |
| Eigen | Linear algebra: `SelfAdjointEigenSolver` for the SYM/ARD Mk ancestral-state reconstruction, and the matrix backend for the vendored LBFGS++ used by `mmvec` | `MPL-2.0` | Header-only, so nothing is linked — the headers are compiled in. `EIGEN_MPL2_ONLY` is defined globally (`CMakeLists.txt`), which confines the build to Eigen's pure-MPL2 subset and keeps no stronger-copyleft file reachable. Authoritative text ships as `COPYING.MPL2` in the vcpkg source tree. A pinned, SHA-verified header-only Eigen 3.4.0 is fetched directly when the vcpkg CONFIG package is absent (the no-vcpkg "Tidy Check" CI lane only). |

Per the policy above, MPL-2.0 is identified by SPDX expression with a pointer to
the upstream licence file rather than inlined. MPL-2.0 obligations attach
per-file to Eigen's own sources, which are unmodified; it imposes no condition on
miint's own Modified-BSD code.

RocksDB — also sourced via vcpkg — is documented separately above (it is
dual `GPL-2.0-only OR Apache-2.0`, used under Apache-2.0).

---

## Validation oracles and algorithm references (not linked)

The projects in this section are **not** dependencies of the built extension —
nothing here is compiled, vendored, or linked. They are listed because miint
code was developed with reference to them, or because they generated committed
golden data used by the parity tests. This mirrors the Biopython entry above,
where a reference-only consultation is likewise recorded.

All seven Python projects are distributed under the 3-clause BSD license; the
shared license text appears once at the end of this section with each project's
copyright notice listed alongside it.

### cogent3 (and its predecessor PyCogent)

The eight community β-diversity metrics in `src/community_distances.cpp` are
independent C++ implementations of the metrics used by Kuczynski et al. 2010,
whose reference implementation is `cogent3.maths.distance_transform` (formerly
`cogent.maths.distance_transform` in PyCogent). No code was copied. Two
zero-variance / zero-row-sum conventions that the published formulas leave
undefined — `dist_pearson` on a constant profile, `dist_chisq` on a zero-sum row
— follow that module's behavior so the reproduction is faithful; both were
verified numerically against cogent3 (agreement ≤ 4.4e-16 across all eight
metrics on `data/simsurvey/beta_distance_fixture.csv`).

- Repository: https://github.com/cogent3/cogent3
- Version referenced: cogent3 2026.7.6a0
- License: `BSD-3-Clause` — Copyright 2019-2021 Gavin Huttley
- Authoritative license text: `LICENSE` in the cogent3 source distribution

#### Citations

Knight, R.; Maxwell, P.; Birmingham, A.; Carnes, M.; Caporaso, J.G.; Easton, B.C.;
Eaton, M.; Hamady, M.; Lindsay, H.; Liu, Z.; Lozupone, C.; McDonald, D.;
Robeson, M.; Sammut, R.; Smit, S.; Wakefield, M.J.; Widmann, J.; Wikman, S.;
Wilson, S.; Ying, H.; and Huttley, G.A. (2007) "PyCogent: a toolkit for making
sense from sequence", Genome Biology, 8(8), R171. doi: 10.1186/gb-2007-8-8-r171

cogent3 is cited via Zenodo: doi: 10.5281/zenodo.15067121

### NumPy

The seeded generator `numpy.random.default_rng(20260817)` produced the input
samples committed in `data/ks2samp/ks_2samp_fixture.csv` (1902 rows), which
`test/sql/ks_2samp_parity.test` reads. NumPy runs only inside
`test/scripts/generate_ks2samp_oracle.py`, when those goldens are regenerated; no
NumPy code is vendored, translated, or executed by miint itself.

- Repository: https://github.com/numpy/numpy
- Version referenced: 2.5.1
- License: `BSD-3-Clause` — Copyright (c) 2005-2025, NumPy Developers.
  (The published NumPy wheel additionally bundles vendored components under
  0BSD, MIT, Zlib and CC0-1.0. None of NumPy is distributed with miint, so only
  its own BSD-3-Clause terms are reproduced here.)

#### Citation

Harris, C.R.; Millman, K.J.; van der Walt, S.J.; Gommers, R.; Virtanen, P.;
Cournapeau, D.; et al. (2020) "Array programming with NumPy", Nature, 585(7825),
357-362. doi: 10.1038/s41586-020-2649-2

### SciPy

`scipy.cluster.hierarchy.linkage(method='average')` and `cophenet` generated the
committed cophenetic goldens in `data/simsurvey/cluster_upgma_oracle.csv`, which
`test/sql/cluster_upgma_parity.test` checks `cluster_upgma` against. SciPy
primitives also generated `data/simsurvey/beta_distance_oracle.csv`.

SciPy is additionally the *behavioral specification* for `Linregress` in
`src/absquant.cpp`. pysyndna fits its models by calling `scipy.stats.linregress`,
so reproducing pysyndna means reproducing scipy — including several behaviors a
textbook least-squares implementation does not have (biased covariances, a
`TINY = 1e-20` guard inside the t statistic, a special case at n = 2, and
`intercept_stderr = stderr·sqrt(ssxm + xmean²)`). Those were determined by
reading scipy's source and confirmed numerically; the divergences are recorded in
the doc comment on `Linregress` and in `docs/absolute_quantification.md`.

`scipy.stats.t.sf` generated `data/syndna/studentt_sf_oracle.csv`, the grid that
`StudentTSurvival` is checked against. The regularized incomplete beta underneath
it is written from the standard continued fraction (DLMF 8.17.22) evaluated by
Lentz's method — **not** transcribed from scipy's vendored Cephes `incbet.c`,
which carries its own separate provenance.
`scipy.stats.ks_2samp` generated the committed parity goldens
`data/ks2samp/ks_2samp_oracle.csv` (32 cases x 2 methods = 64 rows) and
`data/ks2samp/ks_2samp_large_oracle.csv` (5 cases at n = 5000-10000, 10 rows),
which `test/sql/ks_2samp_parity.test` checks miint's `ks_2samp` against. Unusually
for this section the generator itself is committed, at
`test/scripts/generate_ks2samp_oracle.py`: the small grid is drawn from a seeded
NumPy generator, and a seeded draw sequence cannot be reconstructed from a prose
description of the grid, so without the script a reviewer facing a SciPy bump
could not tell a real regression from a different set of draws.

**miint's `ks_2samp` is not derived from SciPy's.** The exact two-sample p-value
is worked out from the lattice-path formulation in Hodges (1958) -- following
Drion and Gnedenko-Korolyuk -- and is computed differently from SciPy's, as the
probability mass escaping an absorbing band rather than as `1 - P(stay inside)`.
SciPy is the oracle the result is checked against, not its source. `method :=
'asymp'` is deliberately not implemented, so none of SciPy's `_ksstats.py` region
selection is reproduced either. The same statement appears next to the algorithm,
in `src/include/KsTwoSample.hpp`.

- Repository: https://github.com/scipy/scipy
- Versions referenced: 1.18.0 (`data/simsurvey/`), 1.17.1 (`data/syndna/`)
- License: `BSD-3-Clause` — Copyright (c) 2001-2002 Enthought, Inc.; 2003-, SciPy Developers

#### Citation

Virtanen, P.; Gommers, R.; Oliphant, T.E.; et al. (2020) "SciPy 1.0: fundamental
algorithms for scientific computing in Python", Nature Methods, 17(3), 261-272.
doi: 10.1038/s41592-019-0686-2

### pysyndna

The three `absquant_*` functions are an independent C++ reimplementation of
pysyndna, which realizes the synDNA spike-in method of Zaramela et al. 2022 for
turning compositional metagenomic read counts into absolute quantities. Each
reimplements one pysyndna entry point:

| miint | pysyndna |
|---|---|
| `absquant_fit_models` | `fit_linear_regression_models` |
| `absquant_cell_counts` | `calc_ogu_cell_counts_biom` |
| `absquant_orf_copies` | `calc_copies_of_ogu_orf_ssrna_per_g_sample_from_dfs` |

They live in `src/absquant.cpp` (the pure core, which also carries `Linregress`
and `StudentTSurvival`), the three DuckDB wrappers
`src/absquant_function.cpp`, `src/absquant_cell_counts_function.cpp` and
`src/absquant_orf_copies_function.cpp`, the shared relation readers in
`src/absquant_readers.cpp`, and the headers `src/include/absquant.hpp` and
`src/include/absquant_readers.hpp`. No code was copied; the port was written from
pysyndna's algorithm as documented and read in its source, and the behaviors it
deliberately does not reproduce are listed under "Differences from pysyndna" in
`docs/absolute_quantification.md`.

pysyndna also generated every input fixture and parity golden under
`data/syndna/` except `studentt_sf_oracle.csv` — including the ORF coordinate,
count and parameter inputs and both ORF goldens (`orf_oracle.csv`,
`orfb_oracle.csv`). It was run **once, offline**, in a dedicated conda
environment pinned to the commit below; only its numbers are committed.
duckdb-miint never invokes pysyndna at build, run, or test time, and does not
depend on it. `data/syndna/README.md` records the provenance and the
regeneration recipe.

- Repository: https://github.com/biocore/pysyndna
- Version referenced: `a64687d4fb37ef7939b1cef8406c0b9758ebb8d7` (version 2026.02.2)
- License: `BSD-3-Clause` — Copyright (c) 2023, Amanda Birmingham
- Note: the repository declares its license in `setup.py` only —
  `license='BSD-3-Clause'` plus the header "Distributed under the terms of the
  Modified BSD License". That header points at "the file LICENSE, distributed
  with this software", but no such file is present in the repository at the
  referenced commit, and no other file carries a license header. The
  declaration in the package metadata is therefore the authoritative statement,
  and the standard 3-clause BSD text at the end of this section is reproduced on
  that basis.

#### Citation

Zaramela, L.S.; Tjuanta, M.; Moyne, O.; Neal, M.; and Zengler, K. (2022)
"synDNA—a Synthetic DNA Spike-in Method for Absolute Quantification of Shotgun
Metagenomic Sequencing", mSystems, 7(6), e00447-22.
doi: 10.1128/msystems.00447-22

### scikit-learn

`sklearn.cluster.KMeans` generated the committed goldens in
`data/simsurvey/cluster_kmeans_oracle.csv`, which
`test/sql/cluster_kmeans_parity.test` checks `cluster_kmeans` against. Its
documented behavior on duplicate points (fewer than `k` non-empty clusters) also
defines the contract recorded in `src/include/cluster_kmeans.hpp`.

- Repository: https://github.com/scikit-learn/scikit-learn
- Version referenced: 1.9.0
- License: `BSD-3-Clause` — Copyright (c) 2007-2024 The scikit-learn developers

#### Citation

Pedregosa, F.; Varoquaux, G.; Gramfort, A.; Michel, V.; Thirion, B.; Grisel, O.;
Blondel, M.; Prettenhofer, P.; Weiss, R.; Dubourg, V.; Vanderplas, J.; Passos, A.;
Cournapeau, D.; Brucher, M.; Perrot, M.; and Duchesnay, E. (2011)
"Scikit-learn: Machine Learning in Python", Journal of Machine Learning Research,
12, 2825-2830.

### scikit-bio (Python package)

Used alongside SciPy to cross-check the community β-diversity metrics when
generating `data/simsurvey/beta_distance_oracle.csv`, and for classical PCoA
during validation. Distinct from the **scikit-bio-binaries (libskbb)** C++
library documented above, which *is* linked into the extension.

`mmvec` (`src/mmvec.cpp`) is an independent C++ reimplementation of
`skbio.stats.ordination.mmvec`, written to reproduce its results. Every expected
value in `test/cpp/mmvec_oracle.hpp` was produced by executing scikit-bio 0.7.3
(the module's sha256 is recorded in that header), and the synthetic input
fixtures in `data/mmvec/` were generated by its private test simulator
`random_multimodal`. The `soils` fixtures were unpivoted from the TSVs committed
in scikit-bio's own test data. No scikit-bio code is vendored, translated
line-by-line, or executed by miint at any point — see `data/mmvec/README.md` for
the full derivation.

- Repository: https://github.com/scikit-bio/scikit-bio
- Version referenced: 0.7.3
- License: `BSD-3-Clause` — Copyright (c) 2013--, scikit-bio development team

#### Citation

Aton, M.; McDonald, D.; Cañardo Alastuey, J.; et al. (2025) "Scikit-bio: a
fundamental Python library for biological omic data analysis", Nature Methods.
doi: 10.1038/s41592-025-02981-z

### mmvec (biocore/mmvec)

The original implementation of the microbe–metabolite co-occurrence method that
scikit-bio's `mmvec` — and therefore miint's — derives from. Consulted as an
algorithm reference; no code is vendored or executed.

The cystic-fibrosis fixtures `data/mmvec/cf_*.parquet` are derived from the BIOM
tables and metadata this project distributes under `examples/cf/`; upstream file
checksums and the derivation are recorded in `data/mmvec/README.md`, along with
the citations for the two studies (`soils`, `cf`) whose data the fixtures carry.

- Repository: https://github.com/biocore/mmvec
- License: `BSD-3-Clause` — Copyright (c) 2018, Jamie Morton

#### Citation

Morton, J. T.; Aksenov, A. A.; Nothias, L. F.; Foulds, J. R.; Quinn, R. A.;
Badri, M. H.; Swenson, T. L.; Van Goethem, M. W.; Northen, T. R.; Vazquez-Baeza,
Y.; Wang, M.; Bokulich, N. A.; Watters, A.; Song, S. J.; Bonneau, R.;
Dorrestein, P. C.; and Knight, R. (2019) "Learning representations of
microbe-metabolite interactions", Nature Methods, 16(12), 1306-1314.
doi: 10.1038/s41592-019-0616-3

### species_abund_sim / ord_survey (Kuczynski et al. 2010 reference code)

`simulate_gradient_otus` and `simulate_cluster_otus`
(`src/simulate_resemblance.cpp`) reproduce the generative models in the
`ord_survey` reference code accompanying Kuczynski et al. 2010. The committed
empirical abundance vectors `data/simsurvey/soil_abundances.csv` (576 species)
and `data/simsurvey/keyboard_abundances.csv` (684 species) were extracted from
that code's `noah_soil_data.py` / `keyboard_data.py`, and the statistical
tolerance bands in `data/simsurvey/gradient_oracle.csv` /
`cluster_oracle.csv` were produced by executing it.

- Author: Justin Kuczynski (2010); the code carries no license statement
- **Used with the author's permission.**

#### Citation

Kuczynski, J.; Liu, Z.; Lozupone, C.; McDonald, D.; Fierer, N.; and Knight, R.
(2010) "Microbial community resemblance methods differ in their ability to detect
biologically relevant patterns", Nature Methods, 7(10), 813-819.
doi: 10.1038/nmeth.1499

### BSD 3-Clause License

Applies to cogent3, NumPy, SciPy, pysyndna, scikit-learn, scikit-bio, and mmvec,
each with its own copyright notice as listed above.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
