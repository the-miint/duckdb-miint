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
`align_minimap2_sharded`, and `save_minimap2_index`.

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
Used by `align_pairwise_score`, `align_pairwise_cigar`, `align_pairwise_full`,
`uchime_ref`, and `uchime_denovo`.

- Repository: https://github.com/smarco/WFA2-lib
- Version: v2.3.5
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

Read-level QC algorithms (adapter trimming, polyG/polyX tail trimming,
sliding-window quality trimming, per-read quality/length/N-base
filtering) used by the `qc_*` scalar function family
(`trim_adapters`, `trim_polyg`, `trim_polyx`, `trim_quality_3p`,
`trim_quality_5p`, `trim_quality_sliding`, `filter_read`). fastp is
**not** a build dependency; specific algorithms have been ported into
`src/qc_algorithms.cpp` with per-file SPDX attribution. Adapter
auto-detection is **not** ported — adapters are user-supplied.

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
