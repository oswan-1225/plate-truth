"""Manuscript assembly: fill {{placeholders}} from computed values and turn [@key] citations
into author-year labels and a reference list built from Crossref records.

A Crossref record is the 'message' of https://api.crossref.org/works/{doi}. Nothing about a
reference (authors, year, title, journal, pages) is typed by hand; it all comes from the record.
"""

import html
import re

PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z0-9_]+)\s*\}\}")
CITATION = re.compile(r"\[(@[A-Za-z0-9_]+(?:\s*;\s*@[A-Za-z0-9_]+)*)\]")


def fill_placeholders(template: str, values: dict[str, str]) -> str:
    """Replace every {{key}} with values[key]; raise listing every key with no value."""
    missing = sorted({m for m in PLACEHOLDER.findall(template) if m not in values})
    if missing:
        raise KeyError(f"no value for placeholders: {', '.join(missing)}")
    return PLACEHOLDER.sub(lambda m: str(values[m.group(1)]), template)


def cited_keys(text: str) -> list[str]:
    """Citation keys in order of first appearance: '[@a; @b] ... [@a]' -> ['a', 'b']."""
    out: list[str] = []
    for group in CITATION.findall(text):
        for k in re.findall(r"@([A-Za-z0-9_]+)", group):
            if k not in out:
                out.append(k)
    return out


def clean(text: str) -> str:
    """Crossref text without markup or repeated whitespace: 'Bac<i>Dive</i> in' -> 'BacDive in'.

    Crossref wraps some titles so that a line break sits before an italic tag inside a word
    ('Bac\\n   <i>Dive</i>'); whitespace containing a line break before an opening tag is
    dropped, while an ordinary space before a tag ('the <i>E. coli</i>') is kept."""
    text = re.sub(r"[ \t]*\n\s*(?=<[a-zA-Z])", "", text)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", text))).strip()


def record_year(rec: dict) -> int:
    """Publication year: print date, else online date, else posted date, else issued."""
    for field in ("published-print", "published-online", "posted", "issued"):
        parts = (rec.get(field) or {}).get("date-parts") or [[None]]
        if parts[0][0]:
            return int(parts[0][0])
    raise ValueError(f"no date in Crossref record {rec.get('DOI')}")


def surnames(rec: dict) -> list[str]:
    """Author surnames in order (an organisation author contributes its name)."""
    out = []
    for a in rec.get("author") or []:
        out.append(a.get("family") or a.get("name") or "")
    if not out or not all(out):
        raise ValueError(f"author list incomplete in Crossref record {rec.get('DOI')}")
    return out


def author_label(rec: dict) -> str:
    """'Martiny', 'Rognes and Flouri' or 'Steen et al.' from the author list."""
    names = surnames(rec)
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return f"{names[0]} and {names[1]}"
    return f"{names[0]} et al."


def citation_labels(records: dict[str, dict]) -> dict[str, str]:
    """Key -> 'Author year' label; works that share author label and year get a, b, ...
    suffixes in order of title, so the labels are unique and stable."""
    base = {k: f"{author_label(r)} {record_year(r)}" for k, r in records.items()}
    out = {}
    for label in set(base.values()):
        keys = [k for k in base if base[k] == label]
        if len(keys) == 1:
            out[keys[0]] = label
            continue
        keys.sort(key=lambda k: clean((records[k].get("title") or [""])[0]).lower())
        for i, k in enumerate(keys):
            out[k] = f"{label}{'abcdefghij'[i]}"
    return out


def replace_citations(text: str, labels: dict[str, str]) -> str:
    """'[@a; @b]' -> '(Label a; Label b)'; raise if a key has no label."""
    def repl(m):
        keys = re.findall(r"@([A-Za-z0-9_]+)", m.group(1))
        unknown = [k for k in keys if k not in labels]
        if unknown:
            raise KeyError(f"no reference for citation keys: {', '.join(unknown)}")
        return "(" + "; ".join(labels[k] for k in keys) + ")"
    return CITATION.sub(repl, text)


def initials(given: str) -> str:
    """'Torbjørn' -> 'T'; 'Jennifer B. H.' -> 'JBH'; 'Jean-Christophe' -> 'J-C'.

    Commas (some Crossref records give 'Koen,J.,van') separate names like spaces, and
    lower-case particles ('van', 'de') are not initials."""
    parts = []
    for word in given.replace(".", " ").replace(",", " ").split():
        if word[0].islower():
            continue
        parts.append("-".join(p[0] for p in word.split("-") if p))
    return "".join(parts)


def format_reference(rec: dict, label_year: str, max_authors: int = 6) -> str:
    """One reference-list entry from a Crossref record.

    'Surname AB, Surname C, ... (2016a) Title. Journal 4:e2584. https://doi.org/...'.
    label_year is the year as it appears in the in-text label (with any a/b suffix)."""
    authors = []
    for a in rec.get("author") or []:
        if a.get("family"):
            ini = initials(a.get("given", ""))
            authors.append(f"{a['family']} {ini}".strip())
        else:
            authors.append(a.get("name", ""))
    if len(authors) > max_authors:
        authors = authors[:max_authors] + ["et al."]
    title = clean((rec.get("title") or [""])[0]).rstrip(".")
    venue = clean((rec.get("container-title") or [""])[0])
    if not venue and rec.get("type") == "posted-content":
        venue = clean(((rec.get("institution") or [{}])[0]).get("name", "")) or "Preprint"
    loc = rec.get("volume", "")
    page = rec.get("page") or rec.get("article-number") or ""
    if page.startswith("10."):  # some records repeat the DOI as the article number
        page = ""
    if loc and page:
        loc = f"{loc}:{page}"
    tail = " ".join(x for x in (venue, loc) if x)
    return (f"{', '.join(authors)} ({label_year}) {title}. {tail + '. ' if tail else ''}"
            f"https://doi.org/{rec['DOI']}")
