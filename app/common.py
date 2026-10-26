"""Cached loaders shared by the app pages."""

from pathlib import Path

import pandas as pd
import streamlit as st

from plate_truth.app_model import load_bundle

APP = Path(__file__).resolve().parent
CAVEATS = """
- These are **predictions**, not observations: the share of reads in taxa a model expects to be
  enriched on plates, for plating conditions like those of ten published studies (soil,
  activated sludge, a marine sponge, oak, broiler carcasses, rumen).
- The model ranks taxa well in a new study, but it **cannot predict a new study's overall level**
  of growth, which differed more than tenfold between the training studies. That is why every
  number is a range from a low-recovery to a high-recovery study, not a single value.
- "Any of 8 plates" assumes plates are independent beyond their shared study and sample, which
  overstates it; read it as a generous upper reference.
- Nothing here describes media designed for a particular organism or habitat, or isolation
  success; plate scrapings measure enrichment of DNA on a plate.
- Identity must be computed the way the model was trained: global alignment to cultured-strain
  16S cut to the same amplified region, **all gaps counted** (VSEARCH `--iddef 1`).
"""


@st.cache_resource
def bundle() -> dict:
    return load_bundle(APP / "model")


@st.cache_data
def habitats(arm: str) -> pd.DataFrame:
    return pd.read_csv(APP / "data" / f"habitats_{arm}.csv")


@st.cache_data
def example() -> pd.DataFrame:
    return pd.read_csv(APP / "data" / "example_soil_sample.csv")


def pct(x: float, digits: int = 1) -> str:
    return f"{100 * x:.{digits}f}%"
