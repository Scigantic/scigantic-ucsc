# scigantic-ucsc

Read UCSC Genome Browser conservation scores, multiz alignments and assembly lists straight from UCSC's public S3 bucket (`s3://genome-browser`). Nothing is downloaded first.

```python
import scigantic_ucsc as ucsc

# per-base phyloP from the 447-way mammal alignment, hg38
scores = ucsc.conservation("chr17", 43_044_295, 43_044_395)

# mean score in 50 windows across a 10 kb region
bins = ucsc.conservation_bins("chr17", 43_044_295, 43_054_295, bins=50)

# alignment blocks overlapping a region, human plus selected species
blocks = ucsc.multiz_blocks("chr1", 1_000_000, 1_000_300, species=["mm10", "canFam3"])

# the same region projected onto hg38 coordinates: one string per species
aln = ucsc.multiz_alignment("chr1", 1_000_000, 1_000_100, species=["mm10", "canFam3"])

# the 1,006 mammal assemblies UCSC lists
ucsc.find("marmoset")

# TOGA gene models and loss status on any assembly in a hub
genes = ucsc.toga_genes("GCF_001458135.1", "CZRN01000001", 0, 5_000_000, status="Lost")
```

## Installation

```console
$ pip install scigantic-ucsc
```

## What it reads

| Function | Source in the bucket | How |
|---|---|---|
| `conservation`, `conservation_bins` | `goldenPath/hg38/phyloP*way`, `phastCons100way` bigWig files | pyBigWig range reads, about 0.1 to 0.5 s per query |
| `multiz_blocks` | `gbdb/hg38/multiz100way/defaultMaf/<chrom>.maf` | binary search over byte-range reads, about 16 requests and a few seconds |
| `multiz_alignment` | same | the blocks above, projected onto hg38 so every row has one column per base |
| `assemblies`, `find` | `hubs/<group>/assemblyList.json` | one small JSON file, cached |
| `toga_genes`, `toga_chromosomes` | `hubs/<GCF or GCA>/.../<accession>/bbi/*TOGA*.bb` | bigBed range reads, about a second per region |

Tracks: `phyloP100way`, `phastCons100way`, `phyloP447way` (default), `phyloP447wayPrimates`, `phyloP447wayLRT`. All are hg38.

## Limits

- `multiz_blocks` covers nine species, not the full 100-way. UCSC's `defaultMaf` files hold only the species the browser shows by default: human, rhesus, dog, elephant, chicken, mouse, frog, zebrafish and sea lamprey. The full 100-way MAFs are plain gzip, which cannot be sliced without reading from the start, so they are not used.
- `toga_genes` coordinates are on the query assembly, not on hg38, and use that assembly's own sequence names (scaffold ids for most). `toga_chromosomes(accession)` lists them. The columns that hold HTML and SVG for the browser popup are left out. TOGA exists for only a minority of assemblies: in a random sample of 20 mammal and 20 bird assemblies, 2 of each had it. `toga_genes` raises if an assembly has none.
- Chromosome names are forgiving (`1`, `chr1`, `MT` all work), but a name that is not in hg38, or a window outside the chromosome, raises a `ValueError` that says so.
- The 241-way and 447-way Cactus alignments are distributed as HAL files (about 0.8 and 1.2 TB) on UCSC's own server, not in this bucket. This package does not read them. The 447-way conservation scores are here because UCSC publishes them as bigWig.
- hg38 only for conservation and alignments. `assemblies` works for any hub group UCSC publishes (`mammals`, `birds`, `fish`, `VGP`, and others).
- A window that needs more than `max_bytes` (64 MB by default) of a MAF raises instead of reading on. Narrow the window.

## Data

The data belongs to UCSC and the groups that produced each track. Cite the Zoonomia consortium for the 241-way and 447-way scores, and the UCSC Genome Browser for the rest.
