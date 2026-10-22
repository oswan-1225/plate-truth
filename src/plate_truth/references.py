"""Check every study's citation against the publisher record and the paper itself.

A study's DOI, first author, year and accession are facts about a paper, so they are compared
with Crossref (https://api.crossref.org/works/{doi}) and with the open-access full text from
Europe PMC, never typed from memory. The study ID must start with the real first author's
surname and the year, so a wrong name cannot hide in an identifier.
"""

import re
import unicodedata

import pandas as pd


def name_key(surname: str) -> str:
    """Lowercase ASCII letters only: 'Ortiz-Severín' -> 'ortizseverin', 'Krüger' -> 'kruger'."""
    ascii_ = unicodedata.normalize("NFKD", surname).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z]", "", ascii_.lower())


def citation_parts(citation: str) -> tuple[str, int | None]:
    """First-author surname and year from 'Surname et al. 2021, Journal' (surname '' if none)."""
    m = re.match(r"^\s*(?P<name>[^\d,]+?)(?:\s+et al\.)?\s+(?P<year>(?:19|20)\d{2})\b", citation)
    if not m:
        return "", None
    return m["name"].strip(), int(m["year"])


def check_study(study: dict, crossref: dict, fulltext: str | None) -> dict:
    """Compare one study's recorded fields with its Crossref record and full text.

    crossref: the 'message' of a Crossref /works response. fulltext: plain text of the paper,
    or None if no open-access text exists. Returns the facts and a list of problems.
    """
    authors = crossref.get("author") or [{}]
    first = authors[0].get("family", "")
    year = (crossref.get("issued", {}).get("date-parts") or [[None]])[0][0]
    cited_name, cited_year = citation_parts(study.get("citation", ""))
    key = f"{name_key(first)}{year}"
    problems = []
    if name_key(cited_name) != name_key(first):
        problems.append(f"citation names '{cited_name}', first author is '{first}'")
    if cited_year != year:
        problems.append(f"citation year {cited_year}, published {year}")
    if not study["id"].startswith(key):
        problems.append(f"study id should start with '{key}'")
    if not str(study.get("study_group", "")).startswith(key):
        problems.append(f"study group should start with '{key}'")
    if fulltext is None:
        problems.append("no open-access full text to confirm the accession")
    elif study["accession"] not in fulltext:
        problems.append(f"accession {study['accession']} not in the paper's text")
    return {"study_id": study["id"], "doi": study["doi"], "first_author": first, "year": year,
            "title": (crossref.get("title") or [""])[0], "citation": study.get("citation"),
            "accession": study["accession"], "ok": not problems,
            "problems": "; ".join(problems)}


def check_table(rows: list[dict]) -> pd.DataFrame:
    """Results as a table, problems first."""
    return pd.DataFrame(rows).sort_values(["ok", "study_id"]).reset_index(drop=True)
