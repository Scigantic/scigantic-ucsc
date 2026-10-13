"""Per-base conservation scores (phyloP, phastCons) from UCSC bigWig files."""
from __future__ import annotations

import threading
from typing import Any

import numpy as np
import numpy.typing as npt

from ._http import url
from ._names import chrom_name

# track name -> published path in the bucket. hg38 only: these are the files
# UCSC publishes under goldenPath/hg38.
TRACKS: dict[str, str] = {
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
    import pyBigWig

    try:
        key = TRACKS[track]
    except KeyError:
        raise ValueError(f"unknown track {track!r}; choose from {sorted(TRACKS)}") from None
    handles = _local.__dict__.setdefault("handles", {})
    if track not in handles:
        bw = pyBigWig.open(url(key))
        if bw is None:
            raise OSError(f"could not open {url(key)}")
        handles[track] = bw
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


def conservation(
    chrom: str, start: int, end: int, track: str = "phyloP447way"
) -> npt.NDArray[np.float64]:
    """Score for every base in chrom:start-end (0-based, end exclusive).

    Bases the track has no value for come back as NaN.
    """
    if end <= start:
        raise ValueError("end must be greater than start")
    bw = _open(track)
    chrom = _check(bw, chrom, start, end)
    values = bw.values(chrom, start, end)
    return np.asarray(values, dtype=np.float64)


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
