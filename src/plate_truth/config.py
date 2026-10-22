"""Locate the repository and load the YAML config files in config/."""

from pathlib import Path

import yaml


def project_root() -> Path:
    """Return the repository root (the folder containing pyproject.toml)."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").exists():
            return parent
    raise FileNotFoundError("pyproject.toml not found above " + str(Path(__file__)))


def load_config(name: str) -> dict:
    """Load config/<name>.yaml as a dict."""
    path = project_root() / "config" / f"{name}.yaml"
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def study_dir(study_id: str, tag: str | None = None) -> Path:
    """A study's pipeline folder in data/interim/phase2. A tag (a reprocessing for a
    sensitivity analysis) gets its own folder, {study}__{tag}, so primary outputs are untouched."""
    return data_path("interim") / "phase2" / (f"{study_id}__{tag}" if tag else study_id)


def cli_option(argv: list[str], name: str) -> str | None:
    """Value of a '--name=value' command-line option, or None."""
    prefix = f"--{name}="
    return next((a[len(prefix):] for a in argv if a.startswith(prefix)), None)


def data_path(key: str) -> Path:
    """Resolve a key from config/paths.yaml to an absolute path."""
    return project_root() / load_config("paths")[key]
