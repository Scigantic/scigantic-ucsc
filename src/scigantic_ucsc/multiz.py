"""Slice the hg38 multiz alignment by position without downloading a chromosome.

UCSC publishes the alignment as one uncompressed MAF per chromosome in two
places under gbdb/hg38/multiz100way: defaultMaf (the nine species the browser
shows by default: human, rhesus, dog, elephant, chicken, mouse, frog, zebrafish,
sea lamprey) and maf (the full 100-way, about nine times larger, for example
70 GB for chr1). A MAF is sorted by reference position and has no index, so a
position is found by binary search over byte-range reads: about 16 requests for
the chr1 defaultMaf file, a few seconds.

The copies under goldenPath/hg38/multiz100way/maf are plain gzip and cannot be
seeked, so they are not used here.
"""
from __future__ import annotations

import re
import urllib.error
from dataclasses import dataclass
from typing import Callable, Iterable

from . import _http
from ._names import chrom_name

SOURCES = {
    "default": "gbdb/hg38/multiz100way/defaultMaf",
    "full": "gbdb/hg38/multiz100way/maf",
}
PREFIX = SOURCES["default"]
_PROBE = 1 << 16
_READ_FIRST = 1 << 17
_READ_MAX = 1 << 20

_HEADER = re.compile(r"(?m)^a score")
_S_LINE = re.compile(r"\ns (\S+)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\S+)")


@dataclass(frozen=True)
class Row:
    src: str
    start: int
    size: int
    strand: str
    src_size: int
    text: str

    @property
    def species(self) -> str:
        return self.src.split(".", 1)[0]


@dataclass(frozen=True)
class Block:
    score: float
    rows: tuple[Row, ...]

    @property
    def ref(self) -> Row:
        return self.rows[0]

    def by_species(self) -> dict[str, Row]:
        return {r.species: r for r in self.rows}


Reader = Callable[[str, int, int], bytes]


def parse_block(text: str) -> Block:
    lines = text.strip("\n").split("\n")
    score = 0.0
    if lines[0].startswith("a"):
        m = re.search(r"score=(\S+)", lines[0])
        score = float(m.group(1)) if m else 0.0
    rows = []
    for ln in lines[1:]:
        if ln.startswith("s "):
            f = ln.split()
            rows.append(Row(f[1], int(f[2]), int(f[3]), f[4], int(f[5]), f[6]))
    return Block(score, tuple(rows))


def _next_block(data: str, off: int) -> tuple[int, int] | None:
    """(offset in data of the first block start at/after off, its ref start)."""
    i = data.find("\na score", off)
    if i < 0:
        return None
    i += 1
    m = _S_LINE.search(data, i)
    if m is None:
        return None
    return i, int(m.group(2))


def _seek(key: str, pos: int, n: int, read: Reader) -> int:
    """Byte offset from which reading forward reaches every block overlapping pos."""
    lo, hi = 0, n
    while hi - lo > 2 * _PROBE:
        mid = (lo + hi) // 2
        data = read(key, mid, min(mid + _PROBE, n - 1)).decode("latin1")
        nb = _next_block(data, 0)
        if nb is None or nb[1] >= pos:
            hi = mid
        else:
            lo = mid
    return lo


def blocks(
    chrom: str,
    start: int,
    end: int,
    species: Iterable[str] | None = None,
    max_bytes: int = 64 << 20,
    source: str = "default",
    *,
    _read: Reader | None = None,
    _size: int | None = None,
) -> list[Block]:
    """Alignment blocks overlapping chrom:start-end (0-based, end exclusive).

    Blocks are returned whole, so they can extend past the requested window.
    `species` keeps only those rows (the human reference row is always kept).
    `source` is "default" (nine species) or "full" (the 100-way, with many more
    rows per block, so windows cost roughly nine times the bytes).
    Raises if the window needs more than `max_bytes` of the file.
    """
    if end <= start:
        raise ValueError("end must be greater than start")
    if source not in SOURCES:
        raise ValueError(f"unknown source {source!r}; choose from {sorted(SOURCES)}")
    chrom = chrom_name(chrom)
    key = f"{SOURCES[source]}/{chrom}.maf"
    read = _read or _http.get_range
    try:
        n = _size if _size is not None else _http.size(key)
    except urllib.error.HTTPError as exc:
        if exc.code in (403, 404):
            raise ValueError(f"no hg38 alignment file for {chrom!r}") from None
        raise
    keep = set(species) if species is not None else None

    off = _seek(key, start, n, read)
    out: list[Block] = []
    carry = ""
    fetched = 0
    step = _READ_FIRST
    done = False
    while off < n and not done:
        if fetched >= max_bytes:
            raise ValueError(f"window needs more than {max_bytes} bytes of {key}; narrow it")
        chunk = read(key, off, min(off + step, n) - 1).decode("latin1")
        step = min(step * 2, _READ_MAX)
        fetched += len(chunk)
        off += len(chunk)
        text = carry + chunk
        parts = text.split("\n\n")
        carry = parts.pop() if off < n else ""
        if off >= n:
            parts = text.split("\n\n")
        for part in parts:
            m = _HEADER.search(part)
            if m is None:
                continue
            b = parse_block(part[m.start():])
            if not b.rows:
                continue
            r = b.ref
            if r.start >= end:
                done = True
                break
            if r.start + r.size <= start:
                continue
            if keep is not None:
                b = Block(b.score, tuple(x for x in b.rows if x is b.rows[0] or x.species in keep))
            out.append(b)
    return out


def alignment(
    chrom: str,
    start: int,
    end: int,
    species: Iterable[str] | None = None,
    max_bytes: int = 64 << 20,
    source: str = "default",
    **kw: object,
) -> dict[str, str]:
    """Species -> aligned string of exactly end-start columns, one per hg38 base.

    The alignment is projected onto the human reference: insertions relative to
    hg38 are dropped, a species with a gap gets '-', and a base that no block
    covers is '.' in every row. Soft-masking (lowercase) is kept.
    """
    bl = blocks(chrom, start, end, species=species, max_bytes=max_bytes, source=source, **kw)  # type: ignore[arg-type]
    width = end - start
    names: list[str] = ["hg38"]
    for b in bl:
        for r in b.rows[1:]:
            if r.species not in names:
                names.append(r.species)
    grid = {n: ["."] * width for n in names}
    for b in bl:
        pos = b.ref.start
        present = {r.species: r.text for r in b.rows}
        for col, base in enumerate(b.ref.text):
            if base == "-":
                continue  # insertion relative to hg38
            if start <= pos < end:
                i = pos - start
                for n in names:
                    t = present.get(n)
                    grid[n][i] = "-" if t is None else t[col]
            pos += 1
    return {n: "".join(v) for n, v in grid.items()}
