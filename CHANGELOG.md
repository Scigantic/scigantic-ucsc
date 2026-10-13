# Changelog

All notable changes to this project are documented here. Versions correspond to PyPI releases.

## 0.1.0 - unreleased

- `conservation` and `conservation_bins`: phyloP and phastCons scores for hg38 from UCSC's bigWig files, read by range.
- `multiz_blocks`: slice the nine-species hg38 multiz alignment by position without downloading a chromosome.
- `multiz_alignment`: the same blocks projected onto hg38, one column per reference base, `.` where no block covers.
- `toga_genes` and `toga_chromosomes`: TOGA gene models and loss classification for any assembly in a UCSC hub, read by region.
- `assemblies` and `find`: the assembly list for any UCSC hub group.
