"""MediaDive (DSMZ cultivation media database) REST client with an on-disk JSON cache.

API documentation: https://mediadive.dsmz.de/doc/index.html (endpoint list read from the page's
apidoc bundle, 2026-10-23). Every response is {"status", "count", "data", ...}; list endpoints
take limit and offset.
"""

import json
import time
from pathlib import Path

from plate_truth.http import get

BASE = "https://mediadive.dsmz.de/rest"


def cached_json(url: str, dest: Path, params: dict | None = None, pause: float = 0.2,
                missing_ok: bool = False) -> dict:
    """Return the JSON body of url, reading dest if it exists, else fetching and saving it.

    A response whose status is not 200 is raised, never cached, so a maintenance page cannot
    be mistaken for data. With missing_ok, a 404 (the record does not exist) is cached as
    {"status": 404, "data": None}.
    """
    dest = Path(dest)
    if dest.exists():
        return json.loads(dest.read_text(encoding="utf-8"))
    resp = get(url, params=params)
    if missing_ok and resp.status_code == 404:
        body = {"status": 404, "data": None}
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(body), encoding="utf-8")
        return body
    resp.raise_for_status()
    body = resp.json()
    if body.get("status") != 200:
        raise RuntimeError(f"{url}: status {body.get('status')}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    part.write_text(json.dumps(body), encoding="utf-8")
    part.replace(dest)
    time.sleep(pause)
    return body


def paged(endpoint: str, out_dir: Path, page_size: int = 500) -> list[dict]:
    """Every record of a list endpoint ('media', 'ingredients'), one cached file per page."""
    records, offset = [], 0
    while True:
        body = cached_json(f"{BASE}/{endpoint}", out_dir / f"{endpoint}_{offset:05d}.json",
                           params={"limit": page_size, "offset": offset})
        records += body["data"]
        offset += page_size
        if offset >= body["count"]:
            return records
