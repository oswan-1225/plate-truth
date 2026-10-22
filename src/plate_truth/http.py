"""Polite HTTP downloads with retries and an on-disk cache."""

import time
from pathlib import Path

import requests

USER_AGENT = "plate-truth/0.1 (research pipeline; https://github.com/oswan-1225)"
RETRY_STATUS = (429, 500, 502, 503, 504)


def get(url: str, params: dict | None = None, retries: int = 5, timeout: int = 120,
        stream: bool = False) -> requests.Response:
    """GET a URL, retrying with exponential backoff on network errors and 429/5xx."""
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=timeout, stream=stream,
                                headers={"User-Agent": USER_AGENT})
            if resp.status_code not in RETRY_STATUS:
                return resp
        except requests.RequestException:
            if attempt == retries - 1:
                raise
        time.sleep(2 ** attempt)
    return resp


def download(url: str, dest: Path, params: dict | None = None, retries: int = 5) -> Path:
    """Stream a URL to dest unless dest already exists; write via a .part file so an
    interrupted download never looks complete."""
    dest = Path(dest)
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    for attempt in range(retries):
        resp = get(url, params=params, stream=True)
        resp.raise_for_status()
        try:
            with part.open("wb") as handle:
                for chunk in resp.iter_content(chunk_size=1 << 20):
                    handle.write(chunk)
            break
        except requests.RequestException:          # connection dropped mid-stream: start over
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)
    part.replace(dest)
    return dest
