"""Per-base conservation scores (phyloP, phastCons) from UCSC bigWig files."""
from __future__ import annotations

import threading
from typing import Any

import numpy as np
import numpy.typing as npt

from ._http import open_remote
from ._names import chrom_name

# track name -> published path in the bucket, or an absolute URL for a file
# hosted elsewhere. hg38 only: the bucket paths are the files UCSC publishes
# under goldenPath/hg38, and the 241-way file is scored on the human assembly
# of the Zoonomia alignment (its chr17 length is hg38's, 83,257,441).
TRACKS: dict[str, str] = {
    # Hosted on UCSC's cgl server, not in the S3 bucket. 21.9 GB; the server
    # accepts range requests, so it is read in place like the others.
    "phyloP241way": "https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.phylop-Homo_sapiens.bigWig",
    "phyloP100way": "goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw",
    "phastCons100way": "goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw",
    "phyloP447way": "goldenPath/hg38/phyloP447way/hg38.phyloP447way.bw",
    "phyloP447wayPrimates": "goldenPath/hg38/phyloP447way/hg38.phyloP447wayPrimates.bw",
    "phyloP447wayLRT": "goldenPath/hg38/phyloP447way/hg38.phyloP447wayLRT.bw",
}


# pyBigWig does not document its handles as thread-safe, so each thread opens
# its own rather than sharing one behind a lock.
_local = threading.local()


def _open(track: str) -> Any:
    try:
        key = TRACKS[track]
    except KeyError:
        raise ValueError(f"unknown track {track!r}; choose from {sorted(TRACKS)}") from None
    handles = _local.__dict__.setdefault("handles", {})
    if track not in handles:
        handles[track] = open_remote(key)
    return handles[track]


def _check(bw: Any, chrom: str, start: int, end: int) -> str:
    if end <= start:
        raise ValueError("end must be greater than start")
    chrom = chrom_name(chrom)
    sizes = bw.chroms()
    if chrom not in sizes:
        raise ValueError(f"{chrom!r} is not in this track's assembly (hg38)")
    if start < 0 or end > sizes[chrom]:
        raise ValueError(f"{chrom}:{start}-{end} is outside 0-{sizes[chrom]}")
    return chrom


MAX_BASES = 10_000_000


def conservation(
    chrom: str, start: int, end: int, track: str = "phyloP447way", max_bases: int = MAX_BASES
) -> npt.NDArray[np.float64]:
    """Score for every base in chrom:start-end (0-based, end exclusive).

    Bases the track has no value for come back as NaN. Windows longer than
    `max_bases` raise; use `conservation_bins` for a summary of a long region.
    """
    if end <= start:
        raise ValueError("end must be greater than start")
    if end - start > max_bases:
        raise ValueError(
            f"{end - start:,} bases is over max_bases={max_bases:,}; "
            "raise max_bases or use conservation_bins"
        )
    bw = _open(track)
    chrom = _check(bw, chrom, start, end)
    # numpy=True avoids building a Python list of floats first, which held
    # about five times the array's size in memory at peak.
    return np.asarray(bw.values(chrom, start, end, numpy=True), dtype=np.float64)


def conservation_bins(
    chrom: str, start: int, end: int, bins: int = 100, track: str = "phyloP447way"
) -> npt.NDArray[np.float64]:
    """Mean score in `bins` equal windows across chrom:start-end. NaN where empty."""
    if end <= start:
        raise ValueError("end must be greater than start")
    bw = _open(track)
    chrom = _check(bw, chrom, start, end)
    out = bw.stats(chrom, start, end, type="mean", nBins=bins)
    return np.array([np.nan if v is None else v for v in out], dtype=np.float64)


def chromosomes(track: str = "phyloP447way") -> dict[str, int]:
    """Chromosome name -> length for the track's assembly."""
    return dict(_open(track).chroms())
