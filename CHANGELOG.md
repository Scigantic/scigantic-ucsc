# Changelog

All notable changes to this project are documented here. Versions correspond to PyPI releases.

## 0.2.0 - unreleased

- New conservation track `phyloP241way`: the Zoonomia 241-way phyloP scores on the human assembly, read by range from UCSC's cgl server (`cgl.gi.ucsc.edu`), not the S3 bucket. A track can now be an absolute URL.

## 0.1.1 - 2026-10-13

- `conservation` reads into numpy directly. A 50 Mb window peaked at about 2 GB of Python memory before; it also now raises above `max_bases` (10 million by default).
- The first open of each remote bigWig or bigBed runs in a subprocess and gives up after 60 s with `TimeoutError`. pyBigWig has no timeout and holds the GIL during the open, so a thread cannot interrupt it. One TOGA file (GCA_001685075.1) does not open within 150 s.
- `multiz_blocks` starts with a 128 KB read after the search and grows it, instead of a fixed 1 MB.
- The PyPI summary now mentions TOGA.

## 0.1.0 - 2026-10-13

- `conservation` and `conservation_bins`: phyloP and phastCons scores for hg38 from UCSC's bigWig files, read by range.
- `multiz_blocks`: slice the nine-species hg38 multiz alignment by position without downloading a chromosome.
- `multiz_alignment`: the same blocks projected onto hg38, one column per reference base, `.` where no block covers.
- `toga_genes` and `toga_chromosomes`: TOGA gene models and loss classification for any assembly in a UCSC hub, read by region.
- `assemblies` and `find`: the assembly list for any UCSC hub group.
