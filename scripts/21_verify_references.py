"""Verify every census study's citation against Crossref and the paper's own text.

For each study in config/census_studies.yaml (included and excluded): the DOI must resolve on
Crossref; the citation's first author and year must match it; the study ID and group must start
with that author's surname and the year; and the accession must appear in the open-access full
text (Europe PMC). Full texts are cached in data/raw/census/fulltext_xml/{PMCID}.txt.

Output: results/reference_check.csv. Exits with an error if any study fails, so a wrong
citation stops the pipeline instead of reaching a table or the manuscript.
"""

import re
import sys

from plate_truth.config import data_path, load_config, project_root
from plate_truth.http import get
from plate_truth.references import check_study, check_table

CROSSREF = "https://api.crossref.org/works/"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"


def crossref_record(doi: str) -> dict:
    """The Crossref 'message' for a DOI (raises if the DOI does not resolve)."""
    resp = get(CROSSREF + doi)
    resp.raise_for_status()
    return resp.json()["message"]


def full_text(doi: str) -> str | None:
    """Plain text of the open-access paper with this DOI (cached), or None if there is none."""
    resp = get(EPMC + "search", params={"query": f'DOI:"{doi}"', "format": "json",
                                         "resultType": "lite"})
    resp.raise_for_status()
    hits = [r for r in resp.json()["resultList"]["result"]
            if (r.get("doi") or "").lower() == doi.lower() and r.get("pmcid")]
    if not hits:
        return None
    cache = data_path("raw") / "census" / "fulltext_xml" / f"{hits[0]['pmcid']}.txt"
    if not cache.exists():
        xml = get(f"{EPMC}{hits[0]['pmcid']}/fullTextXML")
        xml.raise_for_status()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml.text)), encoding="utf-8")
    return cache.read_text(encoding="utf-8")


def main() -> None:
    rows = [check_study(s, crossref_record(s["doi"]), full_text(s["doi"]))
            for s in load_config("census_studies")["studies"]]
    table = check_table(rows)
    table.to_csv(project_root() / "results" / "reference_check.csv", index=False)
    for r in table.itertuples():
        print(f"{'ok  ' if r.ok else 'FAIL'} {r.study_id:34s} {r.first_author} {r.year} "
              f"{r.problems}")
    if not table["ok"].all():
        sys.exit(f"{(~table['ok']).sum()} studies failed the reference check")


if __name__ == "__main__":
    main()
