# scigantic-ucsc

Python readers for four kinds of data in UCSC's public S3 bucket (`s3://genome-browser`, us-east-1, anonymous access): per-base conservation scores, hg38 multiz alignment blocks, TOGA gene annotations, and the assembly lists of UCSC's genome hubs. Each function reads only the byte ranges it needs over HTTPS, so no file is downloaded in full.

```python
import scigantic_ucsc as ucsc

# phyloP scores, hg38, 0-based half-open coordinates
scores = ucsc.conservation("chr17", 43_044_295, 43_044_395, track="phyloP447way")

# mean score in 50 equal windows
bins = ucsc.conservation_bins("chr17", 43_044_295, 43_054_295, bins=50)

# multiz alignment blocks overlapping a region
blocks = ucsc.multiz_blocks("chr1", 1_000_000, 1_000_300, species=["mm10", "canFam3"])

# the same region as one string per species, one column per hg38 base
aln = ucsc.multiz_alignment("chr1", 1_000_000, 1_000_100, species=["mm10", "canFam3"])

# assemblies listed in the mammals hub
ucsc.find("marmoset")

# TOGA genes classified "Lost" in a region of a query assembly
genes = ucsc.toga_genes("GCF_001458135.1", "CZRN01000001", 0, 5_000_000, status="Lost")
```

## Installation

```console
$ pip install scigantic-ucsc
```

Requires Python 3.10 or later, numpy and pyBigWig.

## Functions and sources

| Function | Reads | Method |
|---|---|---|
| `conservation`, `conservation_bins` | bigWig files under `goldenPath/hg38/` | pyBigWig remote range reads |
| `multiz_blocks`, `multiz_alignment` | `gbdb/hg38/multiz100way/defaultMaf/<chrom>.maf` | binary search by position over HTTP byte ranges; the files have no index |
| `assemblies`, `find` | `hubs/<group>/assemblyList.json` | one JSON file per group, cached in memory |
| `toga_genes`, `toga_chromosomes` | the `*TOGA*.bb` bigBed in `hubs/<GCF or GCA>/.../<accession>/bbi/` | pyBigWig remote range reads; the file is found by listing the assembly's `bbi/` prefix |

Conservation tracks: `phyloP100way`, `phastCons100way`, `phyloP447way` (default), `phyloP447wayPrimates`, `phyloP447wayLRT`, and `phyloP241way`. All are scored on hg38 coordinates. `phyloP241way` is the Zoonomia 241-way phyloP file (21.9 GB), hosted on UCSC's `cgl.gi.ucsc.edu` server, not in the S3 bucket; the others are in the bucket.

`toga_genes` returns the reference transcript, the region on the query assembly, UCSC's gene classification (`status`) and the chain and intactness scores stored in the file. It omits the four columns that hold HTML and SVG for the browser popup.

## Measurements

Taken on 2026-10-13 from a single machine, one request at a time.

- Conservation: opening a bigWig took about 1 s the first time. A 20-base query took 0.13 s, and a mean over 10 kb took 0.54 s.
- Multiz on chr1 (a 4.9 GB file): a 300 bp window took 17 requests, 1.2 MB and 3.6 s. A 2 Mb window (80,022 blocks) took 62 requests, 47 MB and 16.4 s.
- TOGA on one marmoset assembly: 239 genes from a 5 Mb region took 1.4 s.
- Multiz search correctness: the position search was compared with a full parse of `chrY.maf` (558,379 blocks) and `chr21.maf` (1,404,006 blocks), 400 random windows each, widths from 1 bp to 200 kb. The block lists were identical in all 800 comparisons.
- `phyloP241way` (read from `cgl.gi.ucsc.edu`): opens in 1.3 s and takes 0.2 s for a 10-base query. Over three 20 kb windows (chr17, chr1, chrX) its scores correlated with `phyloP447way` at r = 0.55 to 0.85 and with `phyloP100way` at r = 0.62 to 0.77, against r near 0 for the 447-way track shifted by 5 kb, so the coordinates match.
- The five bucket conservation tracks returned data for 40 random 1 kb windows each. Between 8% and 17% of bases had no score (NaN).
- The mammals hub listed 1,006 assemblies on 2026-10-13.

## Limits

- The multiz files that can be read this way contain nine species: hg38, rheMac3, canFam3, loxAfr3, galGal4, mm10, xenTro7, danRer10 and petMar2 (the species the UCSC browser shows by default). The full 100-way MAFs under `goldenPath/hg38/multiz100way/maf/` are plain gzip, which cannot be read from the middle, so this package does not use them.
- The Cactus 241-way and 447-way alignments are HAL files of about 0.8 and 1.2 TB on UCSC's own server, not in this bucket. This package does not read them. Only their conservation scores are read, from the bigWig files UCSC publishes.
- Conservation and multiz cover hg38 only. `phyloP241way` is read from `cgl.gi.ucsc.edu`, a single academic server, so it is the track most likely to be slow or unavailable.
- Base-level conservation reads are slow for long windows: a 10 Mb window took 97.5 s (about 10 s per megabase) and returns an 80 MB array. `conservation` raises above `max_bases` (10 million by default). `conservation_bins` summarizes a long region and is the intended call for that.
- Parallel threads do not speed up `conservation`: 300 queries over 32 threads took 108 s, about the same as running them one at a time. Use processes for parallelism.
- TOGA annotations exist for a minority of assemblies. In a random sample of 20 mammal and 20 bird assemblies, 2 of each had one. `toga_genes` raises `ValueError` if an assembly has none. Coordinates are on the query assembly and use its own sequence names, listed by `toga_chromosomes(accession)`.
- One TOGA file in the sampled mammal assemblies (GCA_001685075.1) did not open within 150 s, twice. Its sequence count (58,423) and file size (203 MB) are similar to assemblies that open in under a second, and the cause has not been identified. The first open of each remote bigWig or bigBed runs in a subprocess and raises `TimeoutError` after 60 s.
- `multiz_alignment` projects onto hg38. Insertions relative to hg38 are dropped, a species with a gap gets `-`, and a base that no block covers is `.` in every row.
- `multiz_blocks` raises if a window needs more than `max_bytes` (64 MB by default) of a file.
- Chromosome names `1`, `chr1`, `CHR1` and `MT` are normalized to UCSC's spelling. A name that is not in hg38, a negative start, or a window past the end of the chromosome raises `ValueError`.
- The tests that read the live bucket are skipped unless `SCIGANTIC_UCSC_LIVE=1` is set, so CI does not depend on UCSC or AWS.

## Data and citation

The data is produced by UCSC and by the groups behind each track. The track description pages on genome.ucsc.edu give the citation for each. This package is not affiliated with UCSC.
