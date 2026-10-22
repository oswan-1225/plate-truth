import pandas as pd

from plate_truth.amplicon import cut_window, fastq_urls, reverse_complement

FWD = "GTGCCAGCMGCCGCGGTAA"      # 515F
REV = "GGACTACHVGGGTWTCTAAT"     # 806R
INSERT = "TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGG"


def concrete(primer: str) -> str:
    """One concrete sequence matching a degenerate primer (first base of each code)."""
    pick = {"M": "A", "H": "A", "V": "A", "W": "A", "Y": "C", "R": "A", "N": "A", "B": "C",
            "K": "G", "S": "C", "D": "A"}
    return "".join(pick.get(c, c) for c in primer)


def write_fastq(path, seqs):
    with open(path, "w") as h:
        for i, s in enumerate(seqs):
            h.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


def read_fastq_seqs(path):
    lines = open(path).read().splitlines()
    return lines[1::4]


def test_cut_window_handles_barcode_and_reverse_orientation(tmp_path):
    forward = "ACGTACGTAC" + concrete(FWD) + INSERT + reverse_complement(concrete(REV)) + "TTGA"
    reverse = reverse_complement(forward)
    no_primer = INSERT
    write_fastq(tmp_path / "in.fq", [forward, reverse, no_primer])
    cut_window(tmp_path / "in.fq", tmp_path / "out.fq", FWD, REV, 0.1, 50, 200,
               tmp_path / "rep.json")
    # both orientations give exactly the insert, forward strand; the read without primers is dropped
    assert read_fastq_seqs(tmp_path / "out.fq") == [INSERT, INSERT]


def test_fastq_urls_pair_or_single():
    pair = "ftp.sra.ebi.ac.uk/a/X.fastq.gz;ftp.sra.ebi.ac.uk/a/X_1.fastq.gz;ftp.sra.ebi.ac.uk/a/X_2.fastq.gz"
    assert fastq_urls(pair) == ["https://ftp.sra.ebi.ac.uk/a/X_1.fastq.gz",
                                "https://ftp.sra.ebi.ac.uk/a/X_2.fastq.gz"]
    assert fastq_urls("ftp.sra.ebi.ac.uk/b/Y_subreads.fastq.gz") == [
        "https://ftp.sra.ebi.ac.uk/b/Y_subreads.fastq.gz"]
    assert fastq_urls(pd.NA) == []


def test_asv_end_groups_start_and_end():
    from plate_truth.amplicon import asv_end_groups
    asvs = [("a", "ACGTAA"), ("b", "ACGTCC"), ("c", "TTTTAA")]
    assert asv_end_groups(asvs, 4, "start")["group_id"].tolist() == ["s1", "s1", "s2"]
    # last two bases: a and c share 'AA'
    assert asv_end_groups(asvs, 2, "end")["group_id"].tolist() == ["e1", "e2", "e1"]


def test_split_group_counts_proportional_equal_and_exact_totals():
    from plate_truth.amplicon import split_group_counts
    groups = pd.DataFrame({"asv_id": ["a", "b", "c"], "group_id": ["g1", "g1", "g2"]})
    gc = pd.DataFrame({"run": ["r", "r", "s"], "group_id": ["g1", "g2", "g1"],
                       "count": [10, 4, 5]})
    w = pd.DataFrame({"run": ["r", "r"], "asv_id": ["a", "b"], "count": [3, 1]})
    out = split_group_counts(gc, groups, w).set_index(["run", "asv_id"])["count"]
    # run r: g1 = 10 split 3:1 -> 7.5/2.5 -> largest remainder gives 8 and 2 (sum kept at 10)
    assert out[("r", "a")] + out[("r", "b")] == 10 and out[("r", "a")] in (7, 8)
    assert out[("r", "c")] == 4
    # run s has no weights: equal split of 5 -> 3 and 2
    assert sorted([out[("s", "a")], out[("s", "b")]]) == [2, 3]
