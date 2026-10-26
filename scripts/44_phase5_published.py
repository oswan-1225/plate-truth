"""Phase 5, step 6: published per-habitat estimates, read from the papers by script
(ANALYSIS_PLAN.md section 24), so no number is typed by hand.

- Lloyd et al. 2018 (mSystems, doi:10.1128/mSystems.00055-18), Table 1: metagenome-based shares
  of microbial cells per environment in cultured species-to-genera, uncultured genera-to-classes
  and uncultured phyla (Europe PMC full text PMC6156271).
- Martiny 2019 (ISME J, doi:10.1038/s41396-019-0410-3), abstract: median similarity to cultured
  relatives and the share of sequences with a cultured relative at >97% (PMC6775996 page).
The cached texts are downloaded on the first run into data/raw/phase5_lit/.

Output: results/phase5/published_estimates.csv
"""

import html
import re

import pandas as pd

from plate_truth.config import data_path
from plate_truth.http import download, get

LLOYD = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6156271/fullTextXML"
MARTINY = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6775996/"
# Lloyd et al. environment -> Phase 5 habitat (None: no matching habitat)
LLOYD_HABITAT = {"Marine sediment": "marine_sediment", "Soil": "soil", "Seawater": "marine_water",
                 "Freshwater": "freshwater", "Plant hosts": "plant",
                 "Terrestrial subsurface": None, "Animal hosts": None}


def text(fragment: str) -> str:
    """Plain text of an XML/HTML fragment."""
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", fragment))).strip()


def lloyd_table1(xml: str) -> pd.DataFrame:
    """Rows of Table 1: environment and the three percentages in parentheses."""
    i = xml.index('<table-wrap id="tab1"')
    body = xml[i:xml.index("</tbody>", i)]
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", body, flags=re.S):
        cells = [text(c) for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, flags=re.S)]
        if len(cells) != 5 or not cells[0]:
            continue
        env = re.sub(r"\s*\(\s*\d+\s*\)$", "", cells[0])
        pct = [float(re.search(r"\((\d+)\)", c).group(1)) for c in cells[2:]]
        rows.append({"environment": env, "cultured_species_to_genera": pct[0],
                     "uncultured_genera_to_classes": pct[1], "uncultured_phyla": pct[2]})
    return pd.DataFrame(rows)


def crossref(doi: str) -> dict:
    """First author's surname, year and journal from Crossref (the DOI must resolve)."""
    r = get("https://api.crossref.org/works/" + doi)
    r.raise_for_status()
    m = r.json()["message"]
    year = (m.get("published-print") or m.get("published-online") or m["issued"])
    return {"crossref_first_author": m["author"][0]["family"],
            "crossref_year": year["date-parts"][0][0],
            "crossref_journal": (m.get("container-title") or [""])[0]}


def main() -> None:
    raw = data_path("raw") / "phase5_lit"
    xml = download(LLOYD, raw / "PMC6156271.xml").read_text(encoding="utf-8")
    t1 = lloyd_table1(xml)
    out = []
    for r in t1.itertuples():
        out.append({"source": "Lloyd et al. 2018, Table 1", "doi": "10.1128/mSystems.00055-18",
                    "url": LLOYD, "environment": r.environment,
                    "habitat": LLOYD_HABITAT.get(r.environment),
                    "measure": "percent of cells in cultured species to genera (metagenomes, "
                               "read depths)", "percent": r.cultured_species_to_genera})
    page = text(download(MARTINY, raw / "martiny_pmc.html").read_text(encoding="utf-8",
                                                                    errors="replace"))
    m = re.search(r"median 16S rRNA similarity of bacteria to known cultured relatives was "
                  r"([\d.]+) ± ([\d.]+)%", page)
    s = re.search(r"([\d.]+) ± ([\d.]+)% of sequences and ([\d.]+) ± ([\d.]+)% of taxa", page)
    out.append({"source": "Martiny 2019, abstract", "doi": "10.1038/s41396-019-0410-3",
                "url": MARTINY, "environment": "six biomes pooled", "habitat": None,
                "measure": "percent of sequences with a cultured relative >97% (mean, s.d. "
                           f"{s.group(2)})", "percent": float(s.group(1))})
    out.append({"source": "Martiny 2019, abstract", "doi": "10.1038/s41396-019-0410-3",
                "url": MARTINY, "environment": "six biomes pooled", "habitat": None,
                "measure": f"median similarity to cultured relatives (s.d. {m.group(2)})",
                "percent": float(m.group(1))})
    out = pd.DataFrame(out)
    meta = {d: crossref(d) for d in out["doi"].unique()}
    out = out.join(pd.DataFrame([meta[d] for d in out["doi"]], index=out.index))
    for src, author, year in (("Lloyd", "Lloyd", 2018), ("Martiny", "Martiny", 2019)):
        rows = out[out["source"].str.startswith(src)]
        if not ((rows["crossref_first_author"] == author) & (rows["crossref_year"] == year)).all():
            raise SystemExit(f"Crossref record does not match {author} {year}")
    res = data_path("results") / "phase5"
    res.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(out).to_csv(res / "published_estimates.csv", index=False)
    print(pd.DataFrame(out).to_string(index=False))


if __name__ == "__main__":
    main()
