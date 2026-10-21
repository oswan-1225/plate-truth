"""Phase 1 census, step 1: find candidate studies that sequenced both an environmental sample
and plate-grown communities from it.

Searches Europe PMC with fixed queries (QUERIES below, logged with hit counts), de-duplicates
articles, fetches open-access full text (cached), and records project accessions and keyword
scores used to order candidates for manual screening. Keyword scores never decide inclusion.

Outputs: data/raw/census/ (query results, full texts), data/interim/census/candidates.csv,
data/interim/census/query_log.csv
"""

import json
import re
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from plate_truth.census import (
    bioproject_accessions,
    dedupe_candidates,
    keyword_score,
    record_row,
    search_europepmc,
)
from plate_truth.config import data_path
from plate_truth.http import get

QUERIES = {
    "plate_scrape_methods": 'METHODS:("plate scraping" OR "plate scrapings" OR "scraped from '
    'the plates" OR "scraped from plates" OR "colonies were scraped" OR "plate wash" OR '
    '"plate washes" OR "plate-wash" OR "plate sweep" OR "plate sweeps" OR "pooled colonies") '
    'AND ("16S rRNA") AND (amplicon OR MiSeq OR NovaSeq OR HiSeq OR Illumina)',
    "culture_enriched": '"culture-enriched molecular profiling" OR "culture enriched molecular '
    'profiling" OR "culture-enriched metagenomics"',
    "culturomics_env": '(culturomics OR "culture-omics") AND ("16S rRNA" AND amplicon) AND '
    '(soil OR sediment OR marine OR seawater OR freshwater OR rhizosphere OR plant OR sludge '
    'OR lake)',
    "dep_indep_amplicon": '("culture-dependent" AND "culture-independent") AND METHODS:(agar '
    'AND ("16S rRNA" AND (MiSeq OR Illumina OR amplicon)) AND (scraped OR scraping OR '
    '"plate wash" OR pooled))',
    "cultivable_fraction_seq": '("proportion of cultivable" OR "cultivable fraction" OR '
    '"culturable fraction" OR "cultured fraction") AND ("16S rRNA" AND amplicon) AND '
    '(plates OR agar)',
}


def fulltext(pmcid: str, cache_dir) -> str:
    """Plain text of an open-access article (cached); '' if not available."""
    path = cache_dir / f"{pmcid}.txt"
    if path.exists():
        return path.read_text(encoding="utf-8")
    resp = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", resp.text)) \
        if resp.status_code == 200 else ""
    path.write_text(text, encoding="utf-8")
    return text


def main() -> None:
    raw = data_path("raw") / "census"
    ft_dir = raw / "fulltext"
    ft_dir.mkdir(parents=True, exist_ok=True)
    out = data_path("interim") / "census"
    out.mkdir(parents=True, exist_ok=True)

    rows, log = [], []
    for name, query in QUERIES.items():
        recs = search_europepmc(query)
        (raw / f"query_{name}.json").write_text(json.dumps(recs), encoding="utf-8")
        rows += [dict(record_row(r), query=name) for r in recs]
        log.append({"query_name": name, "query": query, "hits": len(recs),
                    "date": pd.Timestamp.today().date().isoformat()})
        print(f"{name}: {len(recs)}", flush=True)
    pd.DataFrame(log).to_csv(out / "query_log.csv", index=False)

    cand = dedupe_candidates(rows)
    with ThreadPoolExecutor(4) as pool:
        texts = list(pool.map(lambda p: fulltext(p, ft_dir) if isinstance(p, str) else "",
                              cand["pmcid"]))
    cand["has_fulltext"] = [bool(t) for t in texts]
    cand["accessions"] = [";".join(bioproject_accessions(t)) for t in texts]
    ft_scores = pd.DataFrame([keyword_score(t) for t in texts]).add_prefix("ft_")
    ab_scores = pd.DataFrame([keyword_score(f"{a} {b}") for a, b in
                              zip(cand["title"], cand["abstract"], strict=True)]).add_prefix("ab_")
    cand = pd.concat([cand, ft_scores, ab_scores], axis=1)
    cand = cand.sort_values(["ft_score", "ab_score"], ascending=False)
    cand.to_csv(out / "candidates.csv", index=False)
    print(f"unique articles: {len(cand)}; with full text: {cand['has_fulltext'].sum()}; "
          f"with a project accession: {(cand['accessions'] != '').sum()}")


if __name__ == "__main__":
    main()
