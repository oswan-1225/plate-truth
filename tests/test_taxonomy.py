import pandas as pd

from plate_truth.taxonomy import (
    asv_copy_number,
    clean_taxon,
    copy_number_lookup,
    dada2_to_sintax_header,
    parse_sintax,
)


def test_dada2_to_sintax_header():
    h = dada2_to_sintax_header(3, "Bacteria;Bacillota;Bacilli;;")
    assert h == "r3;tax=d:Bacteria,p:Bacillota,c:Bacilli;"


def test_parse_sintax_keeps_passed_ranks(tmp_path):
    p = tmp_path / "s.tsv"
    p.write_text("a1\td:Bacteria(1.00),p:Foo(0.70)\t+\td:Bacteria\n"
                 "a2\td:Bacteria(1.00),g:Bar(0.95)\t+\td:Bacteria,p:Foo,g:Bar\n")
    t = parse_sintax(p)
    assert pd.isna(t.loc["a1", "phylum"])
    assert t.loc["a2", "genus"] == "Bar"


def test_clean_taxon():
    assert clean_taxon("Bacillus_A") == "Bacillus"
    assert clean_taxon("Burkholderia-Caballeronia-Paraburkholderia") == "Burkholderia"


def test_asv_copy_number_falls_back_up_ranks():
    pan = pd.DataFrame({"rank": ["genus", "family"], "name": ["Bacillus", "Xaceae"],
                        "mean": [9.0, 3.0]})
    lookup = copy_number_lookup(pan)
    tax = pd.DataFrame({"family": ["Bacillaceae", "Xaceae", None],
                        "genus": ["Bacillus_A", "Unknownus", None]}, index=["a", "b", "c"])
    out = asv_copy_number(tax, lookup, default=2.0)
    assert out.loc["a"].tolist() == [9.0, "genus"]
    assert out.loc["b"].tolist() == [3.0, "family"]
    assert out.loc["c"].tolist() == [2.0, "default"]
