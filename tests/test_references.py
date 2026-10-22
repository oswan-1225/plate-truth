from plate_truth.references import check_study, citation_parts, name_key

CROSSREF = {"author": [{"family": "Krüger"}, {"family": "Aamand"}],
            "issued": {"date-parts": [[2019, 3, 6]]}, "title": ["Bacterial dispersers"]}
STUDY = {"id": "kruger2019_claytill", "study_group": "kruger2019", "doi": "10.1128/x",
         "citation": "Krüger et al. 2019, Appl Environ Microbiol", "accession": "PRJNA483533"}


def test_name_key_folds_accents_and_hyphens():
    assert name_key("Ortiz-Severín") == "ortizseverin"
    assert name_key("Krüger") == "kruger"


def test_citation_parts():
    assert citation_parts("Krüger et al. 2019, Appl Environ Microbiol") == ("Krüger", 2019)
    assert citation_parts("Martiny 2019, ISME J") == ("Martiny", 2019)
    assert citation_parts("2026, Microbiol Spectrum") == ("", None)


def test_check_study_passes_when_everything_matches():
    out = check_study(STUDY, CROSSREF, "reads are in BioProject PRJNA483533.")
    assert out["ok"] and out["problems"] == ""


def test_check_study_flags_wrong_author_id_and_missing_accession():
    bad = {**STUDY, "id": "lennon2019_claytill", "study_group": "nielsen2019",
           "citation": "Nielsen et al. 2019, Appl Environ Microbiol"}
    out = check_study(bad, CROSSREF, "no accession here")
    assert not out["ok"]
    assert out["problems"].count(";") == 3  # four problems: citation, id, group, accession
    assert "first author is 'Krüger'" in out["problems"]
