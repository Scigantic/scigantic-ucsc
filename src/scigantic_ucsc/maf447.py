"""Slice the Cactus 447-way mammal alignment by hg38 position.

UCSC publishes the 447-way alignment as one bgzip-compressed MAF of about
1.2 TB on its cgl server (cgl.gi.ucsc.edu/data/cactus), referenced on hg38,
with a taffy `.tai` index of about 6 MB. Only the index and the few bgzip
blocks a query needs are fetched, by HTTP range request.

The `.tai` format, as written by `taffy index` and checked against the file:

    <sequence>\\t<ref start>\\t<virtual offset>

one row per index interval of about 10 kb of reference (`hg38.chr1`, ...).
A row whose first column is `*` is on the same sequence as the row above and
holds the DELTA of both numbers from that row. The offset is a bgzip virtual
offset: `compressed_offset << 16 | offset_in_decompressed_block`. A query seeks
to the last row at or before its start and reads blocks forward until the
reference start passes its end, so the cost is the data between that row and
the end of the window, not only the window itself.
"""
from __future__ import annotations

import bisect
import re
import struct
import threading
import zlib
from dataclasses import dataclass
from typing import Callable, Iterable, Iterator

from . import _http
from ._names import chrom_name
from .multiz import Block, Row

HOST = "https://cgl.gi.ucsc.edu/data/cactus"
# variant -> file. "single" has at most one row per species per block.
VARIANTS: dict[str, str] = {
    "single": "447-mammalian-2022v1.fix2.single.maf.gz",
    "all": "447-mammalian-2022v1.fix2.maf.gz",
}
REF = "hg38"

_FIRST_FETCH = 1 << 18
_MAX_FETCH = 1 << 22
_HEADER = re.compile(r"(?m)^a(?: .*)?$")
_REF_S = re.compile(rf"\ns\s+({REF}\.\S+)\s+(\d+)\s+(\d+)")

Reader = Callable[[str, int, int], bytes]


def _parse_block(text: str) -> Block:
    """One MAF block. The 447-way file separates fields with tabs; multiz.py's parser expects spaces."""
    lines = text.strip("\n").split("\n")
    score = 0.0
    m = re.search(r"score=(\S+)", lines[0])
    if m:
        score = float(m.group(1))
    rows = []
    for ln in lines[1:]:
        if ln[:2] in ("s ", "s\t"):
            f = ln.split()
            rows.append(Row(f[1], int(f[2]), int(f[3]), f[4], int(f[5]), f[6]))
    return Block(score, tuple(rows))


def _project(bl: list[Block], start: int, end: int) -> dict[str, str]:
    """Blocks as one string per species, one column per hg38 base in start..end.

    Insertions relative to hg38 are dropped, a species missing from a block gets
    '-', and a base no block covers is '.' in every row.
    """
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
                continue
            if start <= pos < end:
                i = pos - start
                for n in names:
                    t = present.get(n)
                    grid[n][i] = "-" if t is None else t[col]
            pos += 1
    return {n: "".join(v) for n, v in grid.items()}


@dataclass(frozen=True)
class Tai:
    """Parsed `.tai`: per sequence, ascending reference starts and virtual offsets."""

    starts: dict[str, list[int]]
    offsets: dict[str, list[int]]

    def entry(self, seq: str, pos: int) -> tuple[int, int]:
        """(reference start, virtual offset) of the last row at or before `pos`."""
        starts = self.starts[seq]
        k = max(bisect.bisect_right(starts, pos) - 1, 0)
        return starts[k], self.offsets[seq][k]

    def compressed_bytes(self, seq: str, start: int, end: int) -> int:
        """Compressed bytes between the row serving `start` and the row at/after `end`."""
        starts = self.starts[seq]
        a = max(bisect.bisect_right(starts, start) - 1, 0)
        b = min(bisect.bisect_left(starts, end), len(starts) - 1)
        return (self.offsets[seq][b] >> 16) - (self.offsets[seq][a] >> 16)


def parse_tai(text: str) -> Tai:
    starts: dict[str, list[int]] = {}
    offsets: dict[str, list[int]] = {}
    seq = ""
    pos = off = 0
    for line in text.splitlines():
        if not line:
            continue
        name, p, o = line.split("\t")
        if name == "*":
            if not seq:
                raise ValueError("malformed .tai: relative row before any sequence")
            pos += int(p)
            off += int(o)
        else:
            seq, pos, off = name, int(p), int(o)
        starts.setdefault(seq, []).append(pos)
        offsets.setdefault(seq, []).append(off)
    return Tai(starts, offsets)


_tai_cache: dict[str, Tai] = {}
_tai_lock = threading.Lock()


def _key(variant: str) -> str:
    try:
        return f"{HOST}/{VARIANTS[variant]}"
    except KeyError:
        raise ValueError(f"unknown variant {variant!r}; choose from {sorted(VARIANTS)}") from None


def load_tai(variant: str = "single") -> Tai:
    """The index for a variant, fetched once per process (about 6 MB)."""
    key = _key(variant)
    with _tai_lock:
        if key not in _tai_cache:
            _tai_cache[key] = parse_tai(_http.get(key + ".tai").decode("ascii"))
        return _tai_cache[key]


def _bgzf_block(buf: bytes, at: int) -> tuple[bytes, int] | None:
    """Decompress the bgzip block starting at buf[at]; None if it is incomplete."""
    if len(buf) - at < 18:
        return None
    if buf[at : at + 4] != b"\x1f\x8b\x08\x04":
        raise ValueError("not a bgzip block at this offset; the file or index does not match")
    xlen = struct.unpack_from("<H", buf, at + 10)[0]
    i = at + 12
    bsize = 0
    while i < at + 12 + xlen:
        slen = struct.unpack_from("<H", buf, i + 2)[0]
        if buf[i] == 66 and buf[i + 1] == 67:  # 'B', 'C'
            bsize = struct.unpack_from("<H", buf, i + 4)[0] + 1
        i += 4 + slen
    if bsize == 0:
        raise ValueError("bgzip block without a BC field")
    if len(buf) - at < bsize:
        return None
    data = zlib.decompress(buf[at + 12 + xlen : at + bsize - 8], -15)
    crc, isize = struct.unpack_from("<II", buf, at + bsize - 8)
    if isize != len(data) or crc != zlib.crc32(data):
        raise ValueError("bgzip block failed its checksum")
    return data, bsize


def _stream(read: Reader, key: str, size: int, coffset: int, max_bytes: int) -> Iterator[bytes]:
    """Decompressed bgzip blocks from `coffset` onward, fetching a window at a time."""
    buf = b""
    at = 0
    next_fetch = coffset
    fetched = 0
    window = _FIRST_FETCH
    while True:
        got = _bgzf_block(buf, at)
        if got is not None:
            data, bsize = got
            at += bsize
            if data:
                yield data
            continue
        if next_fetch >= size:
            return
        remaining = max_bytes - fetched
        if remaining <= 0:
            raise ValueError(
                f"the window needs more than {max_bytes:,} compressed bytes of {key}; narrow it"
            )
        end = min(next_fetch + min(window, remaining), size) - 1
        chunk = read(key, next_fetch, end)
        fetched += len(chunk)
        next_fetch += len(chunk)
        buf = buf[at:] + chunk
        at = 0
        window = min(window * 2, _MAX_FETCH)


def blocks(
    chrom: str,
    start: int,
    end: int,
    species: Iterable[str] | None = None,
    variant: str = "single",
    max_bytes: int = 64 << 20,
    *,
    _read: Reader | None = None,
    _tai: Tai | None = None,
    _size: int | None = None,
) -> list[Block]:
    """447-way alignment blocks overlapping chrom:start-end (hg38, 0-based, end exclusive).

    Blocks are returned whole, so they can extend past the window. `species`
    keeps only those species' rows (names like "Pan_troglodytes"; the hg38 row
    is always kept). `variant` is "single" (at most one row per species per
    block) or "all". Raises if more than `max_bytes` of COMPRESSED data would
    be read; see the module note on what a query costs.
    """
    if end <= start:
        raise ValueError("end must be greater than start")
    chrom = chrom_name(chrom)
    seq = f"{REF}.{chrom}"
    tai = _tai if _tai is not None else load_tai(variant)
    if seq not in tai.starts:
        raise ValueError(f"no hg38 sequence {chrom!r} in the 447-way index")
    key = _key(variant)
    read = _read or _http.get_range
    size = _size if _size is not None else _http.size(key)
    keep = set(species) if species is not None else None

    _, voff = tai.entry(seq, start)
    coffset, uoffset = voff >> 16, voff & 0xFFFF
    out: list[Block] = []
    tail = ""
    first = True
    done = False
    for data in _stream(read, key, size, coffset, max_bytes):
        text = data.decode("latin1")
        if first:
            text, first = text[uoffset:], False
        pieces = (tail + text).split("\n\n")
        tail = pieces.pop()
        done = _consume(pieces, seq, start, end, keep, out)
        if done:
            break
    if not done and tail:
        _consume([tail], seq, start, end, keep, out)
    return out


def _consume(
    pieces: list[str], seq: str, start: int, end: int, keep: set[str] | None, out: list[Block]
) -> bool:
    """Add overlapping blocks from `pieces` to `out`; True once the window is passed."""
    for piece in pieces:
        m = _HEADER.search(piece)
        if m is None:
            continue
        body = piece[m.start() :]
        s = _REF_S.search(body)
        if s is None:
            continue  # no hg38 row, nothing to place it by
        if s.group(1) != seq:
            return True  # the next sequence: past the end of this one
        rstart, rsize = int(s.group(2)), int(s.group(3))
        if rstart >= end:
            return True
        if rstart + rsize <= start:
            continue
        b = _parse_block(body)
        rows = list(b.rows)
        if not rows[0].src.startswith(REF + "."):
            ref = [r for r in rows if r.src.startswith(REF + ".")]
            rows = ref + [r for r in rows if r not in ref]
        if keep is not None:
            rows = [r for i, r in enumerate(rows) if i == 0 or r.species in keep]
        out.append(Block(b.score, tuple(rows)))
    return False


def alignment(
    chrom: str,
    start: int,
    end: int,
    species: Iterable[str] | None = None,
    variant: str = "single",
    max_bytes: int = 64 << 20,
    **kw: object,
) -> dict[str, str]:
    """Species -> aligned string of exactly end-start columns, one per hg38 base.

    Same projection as `multiz_alignment`: insertions relative to hg38 are
    dropped, a gap in a species is '-', and a base no block covers is '.'.
    """
    bl = blocks(chrom, start, end, species=species, variant=variant, max_bytes=max_bytes, **kw)  # type: ignore[arg-type]
    return _project(bl, start, end)
