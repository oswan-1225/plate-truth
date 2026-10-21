from plate_truth.census import bioproject_accessions, dedupe_candidates, keyword_score


def test_bioproject_accessions():
    text = "Reads are in PRJNA1029248 and PRJEB55863; see also SRP123456 and PRJNA1029248."
    assert bioproject_accessions(text) == ["PRJEB55863", "PRJNA1029248", "SRP123456"]
    assert bioproject_accessions("no accessions, PRJ only") == []


def test_keyword_score():
    s = keyword_score("Colonies were scraped from agar plates; culture-independent DNA too.")
    assert s["plate_terms"] == 2  # 'scraped', 'agar plate'
    assert s["env_terms"] == 1
    assert s["score"] == 5


def test_dedupe_candidates():
    rows = [{"id": "MED:1", "title": "a", "query": "q1"},
            {"id": "MED:1", "title": "a", "query": "q2"},
            {"id": "MED:2", "title": "b", "query": "q1"}]
    out = dedupe_candidates(rows).set_index("id")
    assert out.loc["MED:1", "query"] == "q1;q2" and len(out) == 2


def test_label_runs_uses_sample_field_and_precedence():
    import importlib.util

    import pandas as pd
    spec = importlib.util.spec_from_file_location("s", "scripts/12_census_verify.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    runs = pd.DataFrame({"run_accession": ["r1", "r2", "r3"],
                         "sample_accession": ["s1", "s2", "s3"],
                         "library_strategy": ["AMPLICON"] * 3})
    samples = pd.DataFrame({"sample_accession": ["s1", "s2", "s3"],
                            "description": ["sponges", "scraping colonies on plates", "water"]})
    study = {"id": "x", "match_field": "description", "env_pattern": "^sponges$",
             "plate_pattern": "scraping colonies on plates"}
    out = mod.label_runs(runs, samples, study)
    assert out["role"].tolist() == ["environment", "plate", "other"]
