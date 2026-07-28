# MIINT Documentation

MIINT provides a range of bioinformatic capabilities. The materials included here are coupled where possible with detailed explanation and examples.

## Reading & writing

- [Reading files](reading.md) - File formats with read support in MIINT (FASTA/FASTQ, SAM/BAM, SFF, BIOM, mzML/mzXML, GFF, jplace, Newick).
- [Writing files](writing.md) - File formats with write support in MIINT.

## Public data repositories

- [EBI/ENA](insdc_ena.md) - Interaction with EBI/ENA including submission and fetching of data.
- [NCBI](insdc_ncbi.md) - Interaction with NCBI including fetching of data and the NCBI taxonomy (tree, lineages, retired-taxid remapping).

## Alignment

- [Reference based alignment](alignment_reference.md) - High throughput reference based read alignment.
- [Pairwise alignment](alignment_pairwise.md) - Pairwise alignment methods.
- [Multiple sequence alignment](alignment_multiple.md) - Multiple sequence alignment methods.
- [Alignment analysis](alignment_analysis.md) - Operations for the analysis of alignment data.

## Amplicon & community analysis

- [Quality control](qc.md) - Adapter / polyG / polyX / quality trimming and per-read filtering.
- [Denoising](denoising.md) - Denoise sequencing data into sub-OTU variants.
- [Chimera checking](chimera.md) - Reference-based and de novo chimera detection.
- [Sequence clustering](clustering.md) - Cluster sequences into OTUs at an identity threshold.
- [Sequence search](search.md) - Search sequences against a reference database.
- [Taxonomic classification](classification.md) - Classify sequences by minimizer/k-mer content (RYpe).
- [Profiling & feature tables](profiling.md) - Abundance profiling (sylph) and OGU feature tables (Woltka).

## Phylogeny & diversity

- [Phylogeny estimation](phylogeny.md) - Phylogenetic estimation, tree manipulation (shear, resolve multifurcations/placements), and comparative methods (independent contrasts, ancestral state reconstruction: Brownian-motion, parsimony, and Mk maximum likelihood).
- [Alpha and beta diversity](diversity.md) - Methods to compute and analyze alpha and beta diversity, including non-phylogenetic community distances, ordination, and sample clustering.
- [Community simulation](simulation.md) - Generate synthetic gradient / clustered OTU tables with known ground truth for benchmarking resemblance and ordination methods.

## Mass spectrometry

- [MassQL](massql.md) - Operating on mass spectrometry data through MassQL.

## Other

- [Utilities](utilities.md) - Sequence helpers, versions/diagnostics, and optional-tool installation.
