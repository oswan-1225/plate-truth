"""Locate and run the external command-line tools (VSEARCH, cutadapt)."""

import shutil
import subprocess
import sys
from pathlib import Path

from plate_truth.config import data_path


def vsearch_exe() -> str:
    """Path to VSEARCH: the binary fetched into tools/ (Windows), else one on PATH."""
    found = sorted(data_path("tools").rglob("vsearch.exe"))
    if found:
        return str(found[-1])
    on_path = shutil.which("vsearch")
    if on_path:
        return on_path
    raise FileNotFoundError("VSEARCH not found; run scripts/00_get_vsearch.py")


def cutadapt_exe() -> str:
    """Path to the cutadapt executable installed in the active Python environment."""
    scripts = Path(sys.executable).parent
    for cand in (scripts / "Scripts" / "cutadapt.exe", scripts / "cutadapt"):
        if cand.exists():
            return str(cand)
    on_path = shutil.which("cutadapt")
    if on_path:
        return on_path
    raise FileNotFoundError("cutadapt not found; pip install cutadapt")


def vsearch(*args, quiet: bool = True) -> None:
    """Run VSEARCH with the given arguments; raise if it fails."""
    cmd = [vsearch_exe(), *map(str, args)]
    if quiet:
        cmd.append("--quiet")
    subprocess.run(cmd, check=True)


def cutadapt(*args) -> str:
    """Run cutadapt with the given arguments; return its report text (stdout)."""
    res = subprocess.run([cutadapt_exe(), *map(str, args)], check=True,
                         capture_output=True, text=True)
    return res.stdout
