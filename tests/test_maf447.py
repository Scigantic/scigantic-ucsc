import os
import random
import struct
import zlib

import pytest

from scigantic_ucsc import maf447

BLOCK = 4096  # small bgzip blocks so blocks of text straddle them


def _bgzf(data: bytes, block: int = BLOCK) -> tuple[bytes, list[tuple[int, int]]]:
    """bgzip-compress; returns the file and (compressed offset, uncompressed offset) per block."""
    out = bytearray()
    spans = []
    for i in range(0, len(data), block):
        chunk = data[i : i + block]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        raw = c.compress(chunk) + c.flush()
        bsize = 12 + 6 + len(raw) + 8
        out += b"\x1f\x8b\x08\x04" + b"\0\0\0\0" + b"\0\xff" + struct.pack("<H", 6)
        out += b"BC" + struct.pack("<HH", 2, bsize - 1) + raw
        out += struct.pack("<II", zlib.crc32(chunk), len(chunk))
        spans.append((len(out) - bsize, i))
    eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    return bytes(out) + eof, spans


def _voff(spans: list[tuple[int, int]], upos: int) -> int:
    coff, ustart = [s for s in spans if s[1] <= upos][-1]
    return (coff << 16) | (upos - ustart)


def _build(n_blocks=900, seed=3, interval=10_000):
    """A reference-based MAF for hg38.chr1 and hg38.chr2, its bgzip file and a .tai text."""
    rnd = random.Random(seed)
    species = ["Pan_troglodytes", "Mus_musculus", "Canis_lupus", "Bos_taurus"]
    text = "##maf version=1 scoring=N/A\n"
    marks = []  # (seq, ref start, uncompressed position of the block)
    spans_ref = {}
    for seq in ("chr1", "chr2"):
        pos = 100
        for _ in range(n_blocks // 2):
            size = rnd.randint(20, 300)
            marks.append((f"hg38.{seq}", pos, len(text.encode())))
            ref = "".join(rnd.choice("acgtACGT") for _ in range(size))
            block = f"\na\ns\thg38.{seq}\t{pos}\t{size}\t+\t999999\t{ref}\n"
            for sp in rnd.sample(species, rnd.randint(1, 4)):
                block += f"s\t{sp}.c{rnd.randint(1, 9)}\t{rnd.randint(0, 9999)}\t{size}\t-\t50000\t{ref}\n"
            text += block + "\n"
            spans_ref.setdefault(seq, []).append((pos, size))
            pos += size + rnd.randint(1, 200)
    data = text.encode("latin1")
    file, spans = _bgzf(data)
    # the index writer, written from the taffy description: one row per interval,
    # later rows on the same sequence as deltas of both numbers
    rows, prev = [], {}
    last_pos = {}
    for seq, pos, upos in marks:
        if seq in last_pos and pos - last_pos[seq] < interval:
            continue
        v = _voff(spans, upos)
        if seq in prev:
            rows.append(f"*\t{pos - prev[seq][0]}\t{v - prev[seq][1]}")
        else:
            rows.append(f"{seq}\t{pos}\t{v}")
        prev[seq] = (pos, v)
        last_pos[seq] = pos
    return file, "\n".join(rows) + "\n", spans_ref


FILE, TAI_TEXT, REFS = _build()
TAI = maf447.parse_tai(TAI_TEXT)


def _read(key, a, b):
    return FILE[a : b + 1]


def _blocks(chrom, start, end, **kw):
    return maf447.blocks(chrom, start, end, _read=_read, _tai=TAI, _size=len(FILE), **kw)


def test_tai_deltas_are_decoded_per_sequence():
    t = maf447.parse_tai("hg38.chr1\t0\t0\n*\t10078\t100\n*\t10001\t50\nhg38.chr2\t5\t7\n*\t10\t1\n")
    assert t.starts == {"hg38.chr1": [0, 10078, 20079], "hg38.chr2": [5, 15]}
    assert t.offsets == {"hg38.chr1": [0, 100, 150], "hg38.chr2": [7, 8]}
    with pytest.raises(ValueError):
        maf447.parse_tai("*\t1\t1\n")


def test_blocks_match_the_ground_truth_for_many_windows():
    rnd = random.Random(9)
    for seq in ("chr1", "chr2"):
        refs = REFS[seq]
        for _ in range(60):
            s = rnd.randint(0, refs[-1][0] + 300)
            e = s + rnd.choice([1, 5, 80, 1000, 12_000])
            want = [p for p, n in refs if p < e and p + n > s]
            got = [b.ref.start for b in _blocks(seq, s, e)]
            assert got == want, (seq, s, e)


def test_ref_row_first_species_filter_and_names():
    s, n = REFS["chr1"][200]
    full = _blocks("chr1", s, s + n)[0]
    assert full.ref.src == "hg38.chr1" and full.ref.species == "hg38"
    only = _blocks("chr1", s, s + n, species=["Mus_musculus"])[0]
    assert {r.species for r in only.rows} <= {"hg38", "Mus_musculus"}
    assert only.ref.src == "hg38.chr1"


def test_sequence_boundary_stops_the_scan():
    last_pos, last_size = REFS["chr1"][-1]
    got = _blocks("chr1", last_pos, last_pos + last_size + 10_000_000)
    assert [b.ref.start for b in got][-1] == last_pos
    assert all(b.ref.src == "hg38.chr1" for b in got)


def test_alignment_projects_onto_hg38():
    s, n = REFS["chr1"][50]
    a = maf447.alignment("chr1", s, s + n, _read=_read, _tai=TAI, _size=len(FILE))
    assert set(a) >= {"hg38"} and all(len(v) == n for v in a.values())
    assert "." not in a["hg38"]


def test_budget_checksum_and_bad_inputs():
    with pytest.raises(ValueError, match="narrow it"):
        _blocks("chr1", 0, 10**9, max_bytes=1)
    bad = bytearray(FILE)
    bad[40] ^= 0xFF  # corrupt the first block's compressed data
    with pytest.raises((ValueError, zlib.error)):
        maf447.blocks("chr1", 0, 50, _read=lambda k, a, b: bytes(bad[a : b + 1]), _tai=TAI, _size=len(FILE))
    with pytest.raises(ValueError, match="no hg38 sequence"):
        _blocks("chrZ", 0, 10)
    with pytest.raises(ValueError):
        _blocks("chr1", 10, 10)
    with pytest.raises(ValueError, match="unknown variant"):
        maf447.blocks("chr1", 0, 10, variant="nope", _tai=TAI)


def test_compressed_bytes_estimate_is_monotonic():
    a = TAI.compressed_bytes("hg38.chr1", 0, 20_000)
    b = TAI.compressed_bytes("hg38.chr1", 0, 200_000)
    assert 0 <= a <= b


@pytest.mark.skipif(not os.environ.get("SCIGANTIC_UCSC_LIVE"), reason="live server")
def test_live_447_small_window():
    bl = maf447.blocks("chr17", 43_044_295, 43_044_395)
    assert bl and bl[0].ref.src == "hg38.chr17"
    assert len({r.species for b in bl for r in b.rows}) > 20
