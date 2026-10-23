"""Download MediaDive media, compositions and ingredients (Phase 4, collection media classes).

Every response is cached under data/raw/mediadive/, so rerunning resumes where it stopped.
  media_*.json, ingredients_*.json     the full media and ingredient lists (paged)
  composition/<medium id>.json         /medium-composition/<id>: ingredients with g/L, flattened
                                       over the medium's solutions ({"status": 404} if none)
  ingredient/<ingredient id>.json      /ingredient/<id>: complex_compound flag, formula, roles

Usage: python scripts/32_mediadive_download.py
"""

from concurrent.futures import ThreadPoolExecutor

from plate_truth.config import data_path
from plate_truth.mediadive import BASE, cached_json, paged

WORKERS = 4      # a few parallel requests; each response takes ~3 s


def fetch_all(urls_dests: list[tuple[str, object]], label: str) -> None:
    """Fetch (url, dest) pairs in parallel, reporting progress every 200. A record that does
    not exist (404) is cached as missing."""
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for i, _ in enumerate(pool.map(lambda ud: cached_json(*ud, missing_ok=True),
                                       urls_dests), 1):
            if i % 200 == 0:
                print(f"{label} {i}/{len(urls_dests)}", flush=True)


def main() -> None:
    out = data_path("raw") / "mediadive"
    media = paged("media", out)
    ingredients = paged("ingredients", out)
    print(f"media {len(media)}, ingredients {len(ingredients)}", flush=True)
    fetch_all([(f"{BASE}/ingredient/{i['id']}", out / "ingredient" / f"{i['id']}.json")
               for i in ingredients], "ingredients")
    fetch_all([(f"{BASE}/medium-composition/{m['id']}", out / "composition" / f"{m['id']}.json")
               for m in media], "compositions")
    print("done", flush=True)


if __name__ == "__main__":
    main()
