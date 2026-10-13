"""Byte-range reads against the public UCSC bucket."""
from __future__ import annotations

import time
import urllib.error
import urllib.request
from typing import Any

BASE = "https://genome-browser.s3.amazonaws.com"
_RETRIES = 3


def url(key: str) -> str:
    return f"{BASE}/{key.lstrip('/')}"


def _open(req: urllib.request.Request) -> Any:
    last: Exception | None = None
    for attempt in range(_RETRIES):
        try:
            return urllib.request.urlopen(req, timeout=30)
        except urllib.error.HTTPError as exc:
            # 4xx other than throttling will not get better by asking again.
            if exc.code < 500 and exc.code != 429:
                raise
            last = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
        time.sleep(0.5 * 2**attempt)
    assert last is not None
    raise last


def size(key: str) -> int:
    resp = _open(urllib.request.Request(url(key), method="HEAD"))
    return int(resp.headers["Content-Length"])


def get_range(key: str, start: int, end: int) -> bytes:
    """Bytes start..end inclusive."""
    req = urllib.request.Request(url(key), headers={"Range": f"bytes={start}-{end}"})
    return _open(req).read()  # type: ignore[no-any-return]


def get(key: str) -> bytes:
    return _open(urllib.request.Request(url(key))).read()  # type: ignore[no-any-return]
