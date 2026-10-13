"""TOGA gene annotations for any assembly in a UCSC hub, read by region.

TOGA (Tool to infer Orthologs from Genome Alignments) projects a reference
annotation (human for mammals, chicken for birds) onto a query assembly and
classifies each gene as intact, lost and so on. UCSC mirrors the result as one
bigBed per assembly. Coordinates are on the query assembly, not on hg38.
"""
from __future__ import annotations

import re
import threading
import urllib.request
from dataclasses import dataclass
from functools import lru_cache
from typing import Any

from . import _http

_ACC = re.compile(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.(\d+)$")

# bigBed columns after chrom/start/end, in file order, up to the last scalar.
# The remaining four are HTML and SVG blobs that are tens of kB per gene.
_FIELDS = (
    "name", "score", "strand", "thick_start", "thick_end", "item_rgb", "block_count",
    "block_sizes", "block_starts", "ref_transcript", "ref_region", "query_region",
    "chain_score", "chain_synteny", "chain_flank", "chain_gl_cds_fract",
    "chain_loc_cds_fract", "chain_exon_cov", "chain_intron_cov", "status",
    "perc_intact_ign_missing", "perc_intact_int_missing", "intact_codon_prop",
    "ouf_prop", "mid_intact", "mid_present",
)
_FLOATS = {
    "chain_score", "chain_synteny", "chain_flank", "chain_gl_cds_fract",
    "chain_loc_cds_fract", "chain_exon_cov", "chain_intron_cov",
    "perc_intact_ign_missing", "perc_intact_int_missing", "intact_codon_prop", "ouf_prop",
}


@dataclass(frozen=True)
class Gene:
    chrom: str
    start: int
    end: int
    name: str
    strand: str
    ref_transcript: str
    ref_region: str
    query_region: str
    status: str
    fields: dict[str, Any]

    @property
    def gene_symbol(self) -> str:
        """The reference gene symbol, from names like ENST00000308824.PCMTD2.47."""
        parts = self.name.split(".")
        return parts[1] if len(parts) >= 3 else self.name


def _prefix(accession: str) -> str:
    m = _ACC.match(accession.strip())
    if m is None:
        raise ValueError(f"{accession!r} is not a GenBank or RefSeq accession like GCF_001458135.1")
    return f"hubs/{m[1]}/{m[2]}/{m[3]}/{m[4]}/{accession.strip()}/bbi/"


@lru_cache(maxsize=256)
def toga_key(accession: str) -> str:
    """Bucket key of the TOGA bigBed for an assembly accession."""
    prefix = _prefix(accession)
    xml = urllib.request.urlopen(
        f"{_http.BASE}/?list-type=2&prefix={prefix}", timeout=30
    ).read().decode()
    keys = [k for k in re.findall(r"<Key>([^<]+)</Key>", xml) if "TOGA" in k and k.endswith(".bb")]
    if not keys:
        raise ValueError(f"no TOGA annotation for {accession}")
    return str(keys[0])


_local = threading.local()


def _open(accession: str) -> Any:
    import pyBigWig

    handles = _local.__dict__.setdefault("handles", {})
    if accession not in handles:
        bb = pyBigWig.open(_http.url(toga_key(accession)))
        if bb is None:
            raise OSError(f"could not open the TOGA track for {accession}")
        handles[accession] = bb
    return handles[accession]


def toga_chromosomes(accession: str) -> dict[str, int]:
    """Sequence name -> length for the query assembly."""
    return dict(_open(accession).chroms())


def parse_entry(chrom: str, start: int, end: int, rest: str) -> Gene:
    vals = rest.split("\t")
    f: dict[str, Any] = {}
    for key, raw in zip(_FIELDS, vals):
        f[key] = float(raw) if key in _FLOATS else raw
    return Gene(
        chrom=chrom, start=start, end=end, name=f["name"], strand=f["strand"],
        ref_transcript=f.get("ref_transcript", ""), ref_region=f.get("ref_region", ""),
        query_region=f.get("query_region", ""), status=f.get("status", ""), fields=f,
    )


def toga_genes(
    accession: str, chrom: str, start: int, end: int, status: str | None = None
) -> list[Gene]:
    """TOGA genes overlapping chrom:start-end on the query assembly.

    `status` keeps one classification, for example "Intact" or "Lost" (compared
    case-insensitively). The sequence name is the assembly's own, such as a
    scaffold id; see `toga_chromosomes`.
    """
    if end <= start:
        raise ValueError("end must be greater than start")
    bb = _open(accession)
    sizes = bb.chroms()
    if chrom not in sizes:
        raise ValueError(f"{chrom!r} is not a sequence in {accession}")
    if start < 0 or end > sizes[chrom]:
        raise ValueError(f"{chrom}:{start}-{end} is outside 0-{sizes[chrom]}")
    out = []
    for s, e, rest in bb.entries(chrom, start, end) or []:
        g = parse_entry(chrom, s, e, rest)
        if status is None or g.status.lower() == status.lower():
            out.append(g)
    return out
