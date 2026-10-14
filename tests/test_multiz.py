import os, random

import pytest

from scigantic_ucsc import multiz


def _maf(n_blocks=4000, seed=1):
    rnd = random.Random(seed)
    out = ["##maf version=1 scoring=roast.v3.3\n"]
    pos = 100
    spans = []
    for i in range(n_blocks):
        size = rnd.randint(20, 400)
        seq = "".join(rnd.choice("acgt") for _ in range(size))
        out.append(f"a score={i}.0\n"
                   f"s hg38.chr1 {pos} {size} + 1000000 {seq}\n"
                   f"s mm10.chr2 {pos + 5} {size} - 900000 {seq}\n"
                   f"s canFam4.chr3 {pos + 9} {size} + 800000 {seq}\n\n")
        spans.append((pos, pos + size))
        pos += size + rnd.randint(1, 300)
    return "".join(out).encode(), spans


DATA, SPANS = _maf()


def _read(key, a, b):
    return DATA[a:b + 1]


def _blocks(start, end, **kw):
    return multiz.blocks("chr1", start, end, _read=_read, _size=len(DATA), **kw)


def test_finds_exactly_the_overlapping_blocks():
    start, end = SPANS[2500][0] + 3, SPANS[2510][1] - 3
    got = [b.ref.start for b in _blocks(start, end)]
    want = [s for s, e in SPANS if s < end and e > start]
    assert got == want


def test_window_inside_one_block_and_in_a_gap():
    s, e = SPANS[1000]
    assert [b.ref.start for b in _blocks(s + 1, s + 2)] == [s]
    gap_a, gap_b = SPANS[1000][1], SPANS[1001][0]
    if gap_b > gap_a:
        assert _blocks(gap_a, gap_b) == []


def test_window_at_the_start_and_end_of_the_file():
    assert _blocks(0, SPANS[0][1])[0].ref.start == SPANS[0][0]
    assert _blocks(SPANS[-1][0], SPANS[-1][1] + 10)[-1].ref.start == SPANS[-1][0]


def test_species_filter_keeps_reference():
    b = _blocks(SPANS[10][0], SPANS[10][1], species=["mm10"])[0]
    assert sorted(b.by_species()) == ["hg38", "mm10"]


def test_row_fields():
    r = _blocks(SPANS[5][0], SPANS[5][1])[0].by_species()["mm10"]
    assert (r.strand, r.src_size, r.start) == ("-", 900000, SPANS[5][0] + 5)


def test_budget_is_enforced():
    with pytest.raises(ValueError, match="narrow it"):
        _blocks(0, 10**9, max_bytes=1 << 16)


def test_empty_window_rejected():
    with pytest.raises(ValueError):
        _blocks(10, 10)


@pytest.mark.live
@pytest.mark.skipif(not os.environ.get("SCIGANTIC_UCSC_LIVE"), reason="live bucket")
def test_live_chr1_region():
    bl = multiz.blocks("chr1", 1_000_000, 1_000_300)
    assert bl and bl[0].ref.src.startswith("hg38.chr1")
    assert any(b.ref.start <= 1_000_000 < b.ref.start + b.ref.size or b.ref.start > 1_000_000 for b in bl)


def test_alignment_projects_onto_reference():
    s, e = SPANS[40]
    a = multiz.alignment("chr1", s + 2, e - 2, _read=_read, _size=len(DATA))
    assert set(a) == {"hg38", "mm10", "canFam4"}
    assert all(len(v) == e - s - 4 for v in a.values())
    # mm10 starts 5 bases later than hg38 in the synthetic block, canFam4 9 later
    # but the synthetic rows reuse one sequence, so columns line up one to one
    assert a["hg38"] == a["mm10"] == a["canFam4"]


def test_alignment_marks_uncovered_bases():
    gap_a, gap_b = SPANS[1000][1], SPANS[1001][0]
    a = multiz.alignment("chr1", gap_a - 3, gap_a + (gap_b - gap_a) + 3, _read=_read, _size=len(DATA))
    assert "." in a["hg38"] and a["hg38"].count(".") == gap_b - gap_a


def test_source_selects_the_file_and_is_validated():
    seen = []

    def rd(key, a, b):
        seen.append(key)
        return DATA[a:b + 1]

    multiz.blocks("chr1", SPANS[3][0], SPANS[3][1], _read=rd, _size=len(DATA))
    multiz.blocks("chr1", SPANS[3][0], SPANS[3][1], source="full", _read=rd, _size=len(DATA))
    assert seen[0].startswith("gbdb/hg38/multiz100way/defaultMaf/")
    assert any(k.startswith("gbdb/hg38/multiz100way/maf/") for k in seen)
    with pytest.raises(ValueError, match="unknown source"):
        multiz.blocks("chr1", 0, 10, source="nope", _read=rd, _size=len(DATA))


@pytest.mark.live
@pytest.mark.skipif(not os.environ.get("SCIGANTIC_UCSC_LIVE"), reason="live bucket")
def test_live_full_alignment_has_many_species():
    bl = multiz.blocks("chrM", 0, 16_000, source="full")
    species = {r.species for b in bl for r in b.rows}
    assert "hg38" in species and len(species) > 50
