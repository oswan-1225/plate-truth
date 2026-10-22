import pandas as pd

from plate_truth.reference_build import flag_contaminated_genomes, primer_sets

V4 = {"region": "V4", "primers": "515F/806R", "primer_fwd": "GTGCCAGCMGCCGCGGTAA",
      "primer_rev": "GGACTACHVGGGTWTCTAAT"}


def study(sid, region, primers, fwd, rev):
    return {"id": sid, "region": region, "primers": primers, "primer_fwd": fwd, "primer_rev": rev}


def test_primer_sets_groups_identical_pairs_and_maps_full_length_to_v4():
    studies = [
        {**V4, "id": "a"},
        {**V4, "id": "b"},
        study("c", "full-length", "27F/1492R", "AGRGTTTGATYNTGGCTCAG", "TASGGHTACCTTGTTASGACTT"),
        study("d", "V4", "515F/806R", "GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"),
        study("e", "V3-V4", "341F/806R", "CCTAYGGGRBGCASCAG", "GGACTACNNGGGTATCTAAT"),
    ]
    out = primer_sets(studies, V4)
    assert len(out) == 3                                      # a, b and c share the V4 pair
    assert out.loc[0, "studies"] == "a;b;c"
    # same region and primer names but different sequences get distinct reference names
    assert out["ref_name"].tolist() == ["V4_515F_806R", "V4_515F_806R_2", "V3-V4_341F_806R"]


def test_primer_sets_window_of_shares_the_named_study_set():
    studies = [{**V4, "id": "a"},
               {"id": "r", "region": "V4", "primers": "515F/806R", "primer_fwd": None,
                "primer_rev": None, "window_of": "a"}]
    out = primer_sets(studies, V4)
    assert len(out) == 1 and out.loc[0, "studies"] == "a;r"


def test_flag_contaminated_genomes_drops_whole_genome():
    members = pd.DataFrame({
        "ref_id": ["r1", "r2", "r3", "r4"],
        "member_id": ["gtdb_isolate|G1|0", "gtdb_isolate|G1|1", "gtdb_isolate|G2|2",
                      "refseq_type|NR_1"],
        "source": ["gtdb_isolate", "gtdb_isolate", "gtdb_isolate", "refseq_type"]})
    ssu = pd.Series({"r1": "Bacillota", "r2": "Pseudomonadota", "r3": None,   # parse_sintax names
                     "r4": "Pseudomonadota"})
    genome = pd.Series({"G1": "p__Bacillota", "G2": "p__Bacillota"})
    out = flag_contaminated_genomes(members, ssu, genome).set_index("genome")
    assert out.loc["G1", "n_ssu"] == 2 and out.loc["G1", "n_inconsistent"] == 1
    assert out.loc["G1", "dropped"]                           # one foreign copy drops G1
    assert not out.loc["G2", "dropped"]                       # unclassified copy: not judged
    assert "NR_1" not in out.index                            # only GTDB genomes are checked
