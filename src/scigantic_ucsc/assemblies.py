"""The genome assemblies UCSC lists for a hub group, such as mammals or birds."""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache

from . import _http


@dataclass(frozen=True)
class Assembly:
    tax_id: int
    asm_id: str
    genbank: str | None
    refseq: str | None
    scientific_name: str
    common_name: str
    browser_url: str


@lru_cache(maxsize=8)
def assemblies(group: str = "mammals") -> tuple[Assembly, ...]:
    """Every assembly in hubs/<group> (for example 'mammals', 'birds', 'fish', 'VGP')."""
    raw = json.loads(_http.get(f"hubs/{group}/assemblyList.json"))
    cols = [c["name"] for c in raw["columns"]]
    out = []
    for row in raw["data"]:
        # UCSC writes each row as an object; fall back to positional for a list.
        r = row if isinstance(row, dict) else dict(zip(cols, row))
        out.append(
            Assembly(
                tax_id=r["taxId"],
                asm_id=r["asmId"],
                genbank=r.get("genBank"),
                refseq=r.get("refSeq"),
                scientific_name=r["sciName"],
                common_name=r["comName"],
                browser_url=r["ucscBrowser"],
            )
        )
    return tuple(out)


def find(text: str, group: str = "mammals") -> list[Assembly]:
    """Assemblies whose scientific or common name contains `text` (case-insensitive)."""
    t = text.lower()
    return [a for a in assemblies(group) if t in a.scientific_name.lower() or t in a.common_name.lower()]
