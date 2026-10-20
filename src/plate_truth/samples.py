"""Sample sheets: turn ENA run tables into (source sample, sample kind, condition) rows."""

import pandas as pd

# Zheng et al. 2024 library names, verified against ENA and Supplementary Data 1:
#   {soil}_{n}_ORI_ORI        original sample (PCR replicates)
#   {soil}_{n}_FIL_FIL        5 um filtrate used as the plating inoculum
#   {soil}_{n}_{AE|AN|MI}_{TSA|inSITU}   one plate scraping
ZHENG_PATTERN = r"^(?P<source>[A-Z]+)_(?P<rep>\d+)_(?P<oxygen>[A-Za-z]+)_(?P<medium>[A-Za-z]+)$"
ZHENG_MEDIUM = {"TSA": "1% TSA", "inSITU": "soil extract"}
ZHENG_OXYGEN = {"AE": "aerobic", "AN": "anaerobic", "MI": "microaerophilic"}


def parse_zheng_library(names: pd.Series) -> pd.DataFrame:
    """Split Zheng et al. library names into source, replicate, kind, medium and oxygen.

    kind is 'environment' (ORI), 'inoculum' (FIL) or 'plate'. Names that do not fit the
    pattern get NaN in every column.
    """
    parts = names.str.extract(ZHENG_PATTERN)
    kind = parts["medium"].map({"ORI": "environment", "FIL": "inoculum"}).fillna("plate")
    kind = kind.where(parts["medium"].notna())
    is_plate = kind == "plate"
    return pd.DataFrame({
        "source_sample": parts["source"],
        "replicate": pd.to_numeric(parts["rep"]),
        "kind": kind,
        "medium": parts["medium"].map(ZHENG_MEDIUM).where(is_plate),
        "oxygen": parts["oxygen"].map(ZHENG_OXYGEN).where(is_plate),
    }, index=names.index)


def condition_id(df: pd.DataFrame) -> pd.Series:
    """A readable plating-condition id: source|medium|oxygen (NaN for non-plate rows)."""
    cid = df["source_sample"] + "|" + df["medium"] + "|" + df["oxygen"]
    return cid.where(df["kind"] == "plate")
