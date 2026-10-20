import gzip

from plate_truth.amplicon import reverse_complement
from plate_truth.reference import (
    accession_base,
    ambiguous_fraction,
    dereplicate,
    is_isolate_genome,
    parse_bacdive_record,
    read_fasta,
    sequence_qc,
)


def test_reverse_complement_iupac():
    assert reverse_complement("ACGTN") == "NACGT"
    # 806R GGACTACHVGGGTWTCTAAT: H<->D, V<->B, W stays W
    assert reverse_complement("GGACTACHVGGGTWTCTAAT") == "ATTAGAWACCCBDGTAGTCC"


def test_ambiguous_fraction_and_qc():
    assert ambiguous_fraction("ACGN") == 0.25
    assert sequence_qc(None, 1000, 0.01) == "missing"
    assert sequence_qc("A" * 999, 1000, 0.01) == "too_short"
    assert sequence_qc("A" * 990 + "N" * 11, 1000, 0.01) == "too_ambiguous"  # 11/1001 > 1%
    assert sequence_qc("A" * 1000, 1000, 0.01) == "ok"


def test_is_isolate_genome():
    assert is_isolate_genome("none", "DSM 1234")
    assert not is_isolate_genome("derived from metagenome", "DSM 1234")
    assert not is_isolate_genome("none", "CSMAG_1568")
    assert not is_isolate_genome("none", "SRR123_bin.4")
    assert is_isolate_genome("none", "Robin 7")  # 'bin' inside a word is not a bin name


def test_accession_base():
    assert accession_base("RS_GCF_000331735.1") == "000331735"
    assert accession_base("GCA_000331735") == "000331735"
    assert accession_base(None) == ""


def test_dereplicate_keeps_membership():
    uniq, members = dereplicate([("a", "ACGT"), ("b", "TTTT"), ("c", "ACGT")])
    assert uniq == [("ref_000001", "ACGT"), ("ref_000002", "TTTT")]
    assert members[members.ref_id == "ref_000001"].member_id.tolist() == ["a", "c"]


def test_read_fasta_gz(tmp_path):
    p = tmp_path / "x.fa.gz"
    with gzip.open(p, "wt") as h:
        h.write(">s1 desc\nacg\nt\n>s2\nGG\n")
    assert list(read_fasta(p)) == [("s1 desc", "ACGT"), ("s2", "GG")]


def test_parse_bacdive_record():
    rec = {"Name and taxonomic classification": {"type strain": "yes",
                                                 "LPSN": {"genus": "Foo", "species": "Foo bar"}},
           "Culture and growth conditions": {"culture medium": {"growth": "no"},
                                             "culture temp": [{"growth": "positive"}]},
           "Sequence information": {"16S sequences": {"accession": "X1", "length": 1500},
                                    "Genome sequences": [{"INSDC accession": "GCA_000000001"},
                                                         {"BV-BRC accession": "1.1"}]}}
    out = parse_bacdive_record(7, rec)
    assert out["strain"]["has_growth_record"] is True
    assert out["strain"]["has_medium_and_temp"] is False  # the only medium entry says 'no'
    assert out["strain"]["species"] == "Foo bar"
    assert out["s16"] == [{"bacdive_id": 7, "accession": "X1", "length": 1500}]
    assert out["genomes"] == [{"bacdive_id": 7, "assembly": "GCA_000000001"}]
    empty = parse_bacdive_record(8, {})
    assert empty["strain"]["has_growth_record"] is False and empty["s16"] == []
