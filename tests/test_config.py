from plate_truth.config import data_path, load_config, project_root


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
