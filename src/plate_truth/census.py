"""Phase 1 dataset census: literature search (Europe PMC), keyword screening, accession lookup.

Europe PMC REST API (https://www.ebi.ac.uk/europepmc/webservices/rest/): /search with
cursorMark paging (resultType=core gives abstracts); /{PMCID}/fullTextXML for open-access
full text, from which project accessions are read.
"""

import re

import pandas as pd

from plate_truth.http import get

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Screening keywords: evidence that plates were sequenced as communities, and that a matching
# environmental sample was sequenced too. Used only to ORDER candidates for manual screening.
PLATE_TERMS = [r"plate scrap", r"scraped", r"plate[- ]?wash", r"plate sweep", r"pooled colon",
               r"culture[- ]enriched", r"plate[- ]?based amplicon", r"cultivable taxa",
               r"culturomic", r"agar plate"]
ENV_TERMS = [r"culture[- ]independent", r"original sample", r"direct(ly)? extract",
             r"environmental dna", r"edna", r"inocul", r"source community", r"bulk soil"]


def search_europepmc(query: str, page_size: int = 1000, max_records: int = 5000) -> list[dict]:
    """All records (core: with abstract) for a query, following cursorMark pages."""
    out, cursor = [], "*"
    while len(out) < max_records:
        resp = get(EPMC, params={"query": query, "format": "json", "resultType": "core",
                                 "pageSize": page_size, "cursorMark": cursor})
        resp.raise_for_status()
        data = resp.json()
        batch = data.get("resultList", {}).get("result", [])
        out += batch
        nxt = data.get("nextCursorMark")
        if not batch or not nxt or nxt == cursor:
            break
        cursor = nxt
    return out


def record_row(rec: dict) -> dict:
    """Flatten one Europe PMC core record to the census candidate fields."""
    acc_types = (rec.get("tmAccessionTypeList") or {}).get("accessionType", [])
    return {"id": f"{rec.get('source')}:{rec.get('id')}", "pmid": rec.get("pmid"),
            "pmcid": rec.get("pmcid"), "doi": rec.get("doi"), "title": rec.get("title"),
            "year": rec.get("pubYear"), "journal": (rec.get("journalInfo") or {}).get(
                "journal", {}).get("title"), "pub_type": rec.get("pubType"),
            "abstract": rec.get("abstractText") or "",
            "tm_accession_types": ";".join(sorted(set(acc_types)))}


def keyword_score(text: str) -> dict:
    """Count plate-community and environment-pairing keyword hits in a text (case-insensitive)."""
    t = (text or "").lower()
    plate = sum(bool(re.search(p, t)) for p in PLATE_TERMS)
    env = sum(bool(re.search(p, t)) for p in ENV_TERMS)
    return {"plate_terms": plate, "env_terms": env, "score": plate * 2 + env}


def bioproject_accessions(text: str) -> list[str]:
    """INSDC project/study accessions mentioned in a text (PRJNA/PRJEB/PRJDB, SRP/ERP/DRP)."""
    pat = r"\b(PRJ(?:NA|EB|DB)\d+|[SED]RP\d{6,})\b"
    return sorted(set(re.findall(pat, text or "")))


def fulltext_accessions(pmcid: str) -> list[str]:
    """Project/study accessions in an open-access article's full text (Europe PMC XML).
    Returns [] when no full text is available."""
    resp = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
    if resp.status_code != 200:
        return []
    return bioproject_accessions(re.sub(r"<[^>]+>", " ", resp.text))


def dedupe_candidates(rows: list[dict]) -> pd.DataFrame:
    """One row per article, with the list of queries that found it."""
    df = pd.DataFrame(rows)
    agg = df.groupby("id", as_index=False).agg(
        {**{c: "first" for c in df.columns if c not in ("id", "query")},
         "query": lambda q: ";".join(sorted(set(q)))})
    return agg


ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
RUN_FIELDS = ("run_accession,sample_accession,library_name,library_strategy,library_source,"
              "instrument_platform,read_count,sample_title,sample_alias")
SAMPLE_FIELDS = ("sample_accession,sample_title,description,isolation_source,"
                 "environment_material,host,collection_date,scientific_name")


def ena_table(accession: str, result: str, fields: str, dest) -> "pd.DataFrame":
    """Cached ENA table for a study: result='read_run' (filereport) or 'sample' (search).

    Sample metadata cannot be fetched by project accession from filereport, so the search
    endpoint is used with a study_accession query instead.
    """
    from plate_truth.http import download
    if result == "read_run":
        url, params = ("https://www.ebi.ac.uk/ena/portal/api/filereport",
                       {"accession": accession, "result": result, "format": "tsv",
                        "fields": fields})
    else:
        url, params = (ENA_SEARCH, {"result": result, "format": "tsv", "fields": fields,
                                    "query": f"study_accession={accession}", "limit": 0})
    path = download(url, dest, params=params)
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")
