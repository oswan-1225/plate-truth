import numpy as np
import pandas as pd
import pytest

from plate_truth.identity import (
    best_hits,
    chunk_records,
    nearest_relative_table,
    read_hits,
    search,
    search_chunked,
)
from plate_truth.reference import write_fasta


def _hits():
    return pd.DataFrame({"query": ["q1", "q1", "q1", "q2"],
                         "target": ["r1", "r2", "r3", "r1"],
                         "id1": [98.0, 98.0, 97.0, 90.0],
                         "id2": [98.0, 99.0, 97.0, 90.0],
                         "mid": [99.0, 99.0, 99.5, 95.0]})


def test_best_hits_keeps_ties():
    b = best_hits(_hits(), "id1")
    assert b.loc["q1", "identity"] == 98.0 and b.loc["q1", "n_tied"] == 2
    assert b.loc["q1", "targets"] == "r1;r2"


def test_best_hit_depends_on_definition():
    b = best_hits(_hits(), "mid")
    assert b.loc["q1", "targets"] == "r3"


def test_nearest_relative_table_no_hit_is_nan():
    t = nearest_relative_table(["q1", "q2", "q3"], _hits())
    assert t.loc["q2", "identity_id1"] == 90.0
    assert np.isnan(t.loc["q3", "identity_id1"])


def test_search_counts_internal_gap(tmp_path):
    """Hand-checkable: a 100-base query vs the same sequence with 1 base deleted.
    Alignment = 100 columns, 99 matches, 1 gap: id1 = id2 = 99.0, mid (gaps ignored) = 100."""
    rng = np.random.default_rng(0)
    seq = "".join(rng.choice(list("ACGT"), 100))
    ref = seq[:50] + seq[51:]
    write_fasta([("q", seq)], tmp_path / "q.fa")
    write_fasta([("r", ref)], tmp_path / "r.fa")
    search(tmp_path / "q.fa", tmp_path / "r.fa", tmp_path / "h.tsv", iddef=1, min_id=0.5,
           maxaccepts=0, maxrejects=0, threads=1)
    h = read_hits(tmp_path / "h.tsv").iloc[0]
    assert h["id1"] == pytest.approx(99.0)
    assert h["mid"] == pytest.approx(100.0)


def test_chunk_records_sizes():
    """5 records in chunks of 2: [2, 2, 1], order kept."""
    chunks = chunk_records([(str(i), "A") for i in range(5)], 2)
    assert [len(c) for c in chunks] == [2, 2, 1]
    assert [r[0] for c in chunks for r in c] == ["0", "1", "2", "3", "4"]


def test_search_chunked_equals_one_search_and_resumes(tmp_path):
    """Five queries in chunks of 2 give the same hits as one search; a chunk already on disk
    is reused, not searched again (its marker hit survives into the joined output)."""
    rng = np.random.default_rng(1)
    refs = [("r%d" % i, "".join(rng.choice(list("ACGT"), 120))) for i in range(4)]
    queries = [("q%d" % i, refs[i % 4][1][:60] + "T" + refs[i % 4][1][61:]) for i in range(5)]
    write_fasta(refs, tmp_path / "r.fa")
    write_fasta(queries, tmp_path / "q.fa")
    args = dict(iddef=1, min_id=0.75, maxaccepts=0, maxrejects=0, threads=1,
                top_hits_only=True)
    search(tmp_path / "q.fa", tmp_path / "r.fa", tmp_path / "one.tsv", **args)
    search_chunked(tmp_path / "q.fa", tmp_path / "r.fa", tmp_path / "ch.tsv", chunk_size=2,
                   **args)
    key = ["query", "target"]
    one = read_hits(tmp_path / "one.tsv").sort_values(key).reset_index(drop=True)
    ch = read_hits(tmp_path / "ch.tsv").sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(one, ch)
    assert not (tmp_path / "ch_chunks2").exists()
    work = tmp_path / "re_chunks2"
    work.mkdir()
    (work / "hits_0000.tsv").write_text("marker\tr0\t1\t1\t1\t1\t1\t0\t0\t1\t1\t1\n")
    search_chunked(tmp_path / "q.fa", tmp_path / "r.fa", tmp_path / "re.tsv", chunk_size=2,
                   **args)
    re_hits = read_hits(tmp_path / "re.tsv")
    assert "marker" in set(re_hits["query"]) and "q0" not in set(re_hits["query"])
