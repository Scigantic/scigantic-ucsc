"""Chromosome name handling shared by the readers."""
from __future__ import annotations


def chrom_name(chrom: str) -> str:
    """'1', 'chr1', 'CHR1' -> 'chr1'; 'MT', 'M', 'chrMT' -> 'chrM'."""
    c = chrom.strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    if c.upper() in ("M", "MT"):
        return "chrM"
    if not c:
        raise ValueError("empty chromosome name")
    # UCSC spells X, Y and alt contigs with the chr prefix and keeps the rest.
    return "chr" + (c.upper() if c.upper() in ("X", "Y") or c.isdigit() else c)
