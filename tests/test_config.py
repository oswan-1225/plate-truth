from plate_truth.config import cli_option, data_path, load_config, project_root, study_dir


def test_project_root_contains_pyproject():
    assert (project_root() / "pyproject.toml").exists()


def test_grew_defaults_match_project_definition():
    grew = load_config("thresholds")["grew"]
    assert grew["min_fold_enrichment"] == 2.0
    assert grew["min_reads"] == 5


def test_detected_is_looser_than_grew():
    cfg = load_config("thresholds")
    assert cfg["detected"]["min_reads"] < cfg["grew"]["min_reads"]


def test_identity_bin_edges_are_increasing():
    edges = load_config("thresholds")["identity"]["bin_edges_percent"]
    assert edges == sorted(edges)
    assert load_config("thresholds")["identity"]["database_proxy_percent"] in edges


def test_data_path_resolves_under_root():
    assert data_path("raw") == project_root() / "data" / "raw"


def test_cli_option():
    argv = ["x.py", "study", "--tag=maxee2", "--max-ee=2", "--per-run-only"]
    assert cli_option(argv, "tag") == "maxee2"
    assert cli_option(argv, "max-ee") == "2"
    assert cli_option(argv, "missing") is None


def test_study_dir_tag_gets_own_folder():
    assert study_dir("s").name == "s"
    assert study_dir("s", "maxee2").name == "s__maxee2"
    assert study_dir("s", "maxee2").parent == data_path("interim") / "phase2"
