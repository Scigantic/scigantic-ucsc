import os

import numpy as np
import pytest

import scigantic_ucsc as u

live = pytest.mark.skipif(not os.environ.get("SCIGANTIC_UCSC_LIVE"), reason="live bucket")


def test_unknown_track():
    with pytest.raises(ValueError, match="unknown track"):
        u.conservation("chr1", 0, 10, track="nope")


def test_empty_window():
    with pytest.raises(ValueError):
        u.conservation("chr1", 5, 5)


@live
def test_live_conservation():
    v = u.conservation("chr1", 1_000_000, 1_000_020)
    assert v.shape == (20,) and np.isfinite(v).any()
    assert u.conservation_bins("chr17", 43_044_295, 43_054_295, bins=10).shape == (10,)


@live
def test_live_assemblies():
    a = u.assemblies("mammals")
    assert len(a) > 500
    assert any("aardvark" in x.common_name for x in u.find("aardvark"))


def test_chrom_names():
    from scigantic_ucsc._names import chrom_name
    assert [chrom_name(x) for x in ("1", "chr1", "CHR1", "x", "MT", "chrM", "chrUn_KI270302v1")] == [
        "chr1", "chr1", "chr1", "chrX", "chrM", "chrM", "chrUn_KI270302v1"]
    with pytest.raises(ValueError):
        chrom_name(" ")


@live
def test_live_bad_inputs_are_readable_errors():
    with pytest.raises(ValueError, match="not in this track"):
        u.conservation("chrZ", 0, 10)
    with pytest.raises(ValueError, match="outside"):
        u.conservation("chr1", -5, 10)
    with pytest.raises(ValueError, match="outside"):
        u.conservation("chr1", 248956400, 248956500)
    with pytest.raises(ValueError, match="no hg38 alignment"):
        u.multiz_blocks("chrZ", 0, 10)
    assert u.conservation("1", 1_000_000, 1_000_010).shape == (10,)


def test_max_bases_checked_before_any_network():
    with pytest.raises(ValueError, match="max_bases"):
        u.conservation("chr1", 0, 20_000_000)
