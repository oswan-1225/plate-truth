"""MGnify API access for the Phase 5 survey census (cached JSON) and streamed read screens.

API documentation: https://www.ebi.ac.uk/metagenomics/api/v2/openapi.json (v2: studies by
biome lineage, a study's analyses) and the v1 API (per-analysis summary and download list).
Every response is cached under data/raw/mgnify/api/, so a rerun repeats the same census.
"""

import hashlib
import json
import zlib
from pathlib import Path

from plate_truth.http import get

V1 = "https://www.ebi.ac.uk/metagenomics/api/v1"
V2 = "https://www.ebi.ac.uk/metagenomics/api/v2"


def cached_json(url: str, params: dict | None, cache: Path) -> dict:
    """GET JSON, cached on disk by URL and parameters."""
    key = hashlib.sha1((url + json.dumps(params or {}, sort_keys=True)).encode()).hexdigest()
    path = Path(cache) / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    resp = get(url, params=params)
    resp.raise_for_status()
    data = resp.json()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")
    return data


def paged_items(url: str, params: dict, cache: Path, max_pages: int = 1000) -> list[dict]:
    """All items of a paged v2 list endpoint ({count, items})."""
    items, page = [], 1
    while page <= max_pages:
        d = cached_json(url, {**params, "page": page, "page_size": 100}, cache)
        items += d["items"]
        if len(items) >= d["count"] or not d["items"]:
            break
        page += 1
    return items


def studies_by_lineage(lineage: str, cache: Path) -> list[dict]:
    """MGnify studies whose biome is this lineage or a descendant (v2)."""
    return paged_items(f"{V2}/studies/", {"biome_lineage": lineage}, cache)


def study_analyses(accession: str, cache: Path, max_pages: int = 30) -> list[dict]:
    """A study's analyses (v2): experiment type, pipeline version, run and sample."""
    return paged_items(f"{V2}/studies/{accession}/analyses/", {}, cache, max_pages)


def analysis_rows(items: list[dict]) -> list[dict]:
    """Flatten v2 analysis items to accession, experiment type, pipeline version, run, sample."""
    out = []
    for a in items:
        run = a.get("run") or {}
        sample = a.get("sample") or {}
        out.append({"analysis": a["accession"], "experiment_type": a.get("experiment_type"),
                    "pipeline_version": a.get("pipeline_version"),
                    "run_accession": run.get("accession"),
                    "sample_accession": sample.get("accession") or run.get("sample_accession")})
    return out


def analysis_summary(accession: str, cache: Path) -> dict:
    """The v1 analysis-summary of one analysis as {key: value}."""
    d = cached_json(f"{V1}/analyses/{accession}", None, cache)
    return {x["key"]: x["value"] for x in d["data"]["attributes"].get("analysis-summary", [])}


def ssu_fasta_url(accession: str, cache: Path) -> str | None:
    """URL of the 'encoding SSU rRNA' FASTA of a v1 analysis, if published."""
    d = cached_json(f"{V1}/analyses/{accession}/downloads", {"page_size": 100}, cache)
    for x in d["data"]:
        alias = x["attributes"]["alias"]
        label = x["attributes"]["description"]["label"]
        if alias.endswith("_SSU.fasta.gz") and "SSU rRNA" in label:
            return x["links"]["self"]
    return None


def fastq_sequences(lines, n: int) -> list[str]:
    """The sequence lines of the first n FASTQ records in an iterable of text lines."""
    out = []
    for i, line in enumerate(lines):
        if i % 4 == 1:
            out.append(line.strip())
            if len(out) >= n:
                break
    return out


def stream_fastq_head(url: str, n: int) -> list[str]:
    """Stream a gzipped FASTQ over HTTP and return the first n sequences (no full download)."""
    resp = get(url, stream=True)
    resp.raise_for_status()
    dec = zlib.decompressobj(zlib.MAX_WBITS | 32)
    buf, lines = "", []
    try:
        for chunk in resp.iter_content(1 << 16):
            buf += dec.decompress(chunk).decode("ascii", "replace")
            *full, buf = buf.split("\n")
            lines += full
            if len(lines) >= 4 * n + 4:
                break
    finally:
        resp.close()
    return fastq_sequences(lines, n)
