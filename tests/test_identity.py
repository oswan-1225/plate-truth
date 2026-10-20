import numpy as np
import pandas as pd
import pytest

from plate_truth.identity import best_hits, nearest_relative_table, read_hits, search
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
