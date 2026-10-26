import pytest

from plate_truth.manuscript import (
    author_label,
    citation_labels,
    cited_keys,
    clean,
    fill_placeholders,
    format_reference,
    initials,
    record_year,
    replace_citations,
)


def rec(doi, families, year, title="T", **extra):
    r = {"DOI": doi, "author": [{"family": f, "given": "Ann B."} for f in families],
         "issued": {"date-parts": [[year]]}, "title": [title]}
    r.update(extra)
    return r


def test_fill_placeholders():
    assert fill_placeholders("a {{x}} b {{ y }}", {"x": "1", "y": "2"}) == "a 1 b 2"
    with pytest.raises(KeyError, match="y, z"):
        fill_placeholders("{{x}} {{z}} {{y}}", {"x": "1"})


def test_cited_keys_order_and_dedup():
    assert cited_keys("x [@b; @a] y [@b] z [@c]") == ["b", "a", "c"]


def test_clean_strips_markup():
    assert clean("Bac\n   <i>Dive</i>\n   in 2022 &amp; more") == "BacDive in 2022 & more"
    assert clean("the <i>E. coli</i> genome") == "the E. coli genome"


def test_record_year_prefers_print_date():
    r = rec("d", ["A"], 2012, **{"published-print": {"date-parts": [[2013, 1]]}})
    assert record_year(r) == 2013
    assert record_year(rec("d", ["A"], 2012)) == 2012


def test_author_label():
    assert author_label(rec("d", ["Martiny"], 2019)) == "Martiny"
    assert author_label(rec("d", ["A", "B"], 2019)) == "A and B"
    assert author_label(rec("d", ["Steen", "B", "C"], 2019)) == "Steen et al."


def test_citation_labels_disambiguate_by_title():
    records = {"u": rec("1", ["Edgar"], 2016, "UNOISE2"),
               "s": rec("2", ["Edgar"], 2016, "SINTAX"),
               "m": rec("3", ["Martiny"], 2019)}
    labels = citation_labels(records)
    assert labels == {"s": "Edgar 2016a", "u": "Edgar 2016b", "m": "Martiny 2019"}


def test_replace_citations():
    labels = {"a": "A 2019", "b": "B et al. 2018"}
    assert replace_citations("see [@a; @b].", labels) == "see (A 2019; B et al. 2018)."
    with pytest.raises(KeyError):
        replace_citations("[@zz]", labels)


def test_initials():
    assert initials("Torbjørn") == "T"
    assert initials("Jennifer B. H.") == "JBH"
    assert initials("Jean-Christophe") == "J-C"
    assert initials("Koen,J.,van") == "KJ"


def test_format_reference():
    r = rec("10.7717/peerj.2584", ["Rognes", "Flouri"], 2016, "VSEARCH: a tool.",
            **{"container-title": ["PeerJ"], "volume": "4", "article-number": "e2584"})
    assert format_reference(r, "2016") == ("Rognes AB, Flouri AB (2016) VSEARCH: a tool. "
                                           "PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584")
    many = rec("x", list("ABCDEFG"), 2020, type="posted-content",
               institution=[{"name": "bioRxiv"}])
    assert format_reference(many, "2020").startswith("A AB, B AB, C AB, D AB, E AB, F AB, et al.")
    assert "bioRxiv" in format_reference(many, "2020")
    doi_number = rec("10.1128/x", ["L"], 2018, "T", **{"container-title": ["mSystems"],
                                                       "volume": "3",
                                                       "article-number": "10.1128/x"})
    assert "mSystems 3. https" in format_reference(doi_number, "2018")
