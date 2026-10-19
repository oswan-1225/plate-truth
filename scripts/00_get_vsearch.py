"""Download the official VSEARCH Windows binary into tools/.

VSEARCH is not packaged for win-64 on conda. Ported from the soil pilot, where this exact
asset was used. On Linux/macOS install it with conda instead: bioconda::vsearch.
Release page: https://github.com/torognes/vsearch/releases
Usage: python scripts/00_get_vsearch.py
"""

import io
import platform
import sys
import zipfile

import requests

from plate_truth.config import data_path

VERSION = "2.32.0"
URL = (
    f"https://github.com/torognes/vsearch/releases/download/v{VERSION}/"
    f"vsearch-{VERSION}-win-x86_64.zip"
)


def find_binary(tools_dir):
    """Return the first vsearch.exe under tools_dir, or None."""
    return next(tools_dir.rglob("vsearch.exe"), None)


def main() -> int:
    if platform.system() != "Windows":
        print("Not Windows: install with `conda install -c bioconda vsearch`.")
        return 1
    tools_dir = data_path("tools")
    exe = find_binary(tools_dir)
    if exe is None:
        tools_dir.mkdir(exist_ok=True)
        response = requests.get(URL, timeout=120)
        response.raise_for_status()
        zipfile.ZipFile(io.BytesIO(response.content)).extractall(tools_dir)
        exe = find_binary(tools_dir)
    print("vsearch at", exe)
    return 0


if __name__ == "__main__":
    sys.exit(main())
