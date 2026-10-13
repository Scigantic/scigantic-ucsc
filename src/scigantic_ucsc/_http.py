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


_PROBE = (
    "import sys, pyBigWig\n"
    "b = pyBigWig.open(sys.argv[1])\n"
    "sys.exit(0 if b is not None else 1)\n"
)
_probed: set[str] = set()


def open_remote(key: str, timeout: float = 60.0) -> Any:
    """pyBigWig.open on a bucket key, giving up after `timeout` seconds.

    pyBigWig has no timeout and holds the GIL for the whole open, so a watchdog
    thread cannot interrupt it. The first open of each key is therefore tried in
    a short-lived subprocess that can be killed. At least one TOGA file in the
    bucket (GCA_001685075.1) does not open within minutes. A key that opened
    once is not probed again.
    """
    import subprocess
    import sys

    import pyBigWig

    if key not in _probed:
        try:
            r = subprocess.run(
                [sys.executable, "-c", _PROBE, url(key)],
                timeout=timeout, capture_output=True,
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(f"opening {key} took longer than {timeout:.0f}s") from None
        if r.returncode != 0:
            raise OSError(f"could not open {url(key)}")
        _probed.add(key)
    bw = pyBigWig.open(url(key))
    if bw is None:
        raise OSError(f"could not open {url(key)}")
    return bw
