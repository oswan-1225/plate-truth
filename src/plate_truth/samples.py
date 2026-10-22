"""Sample sheets: turn ENA run tables into (source sample, sample kind, condition) rows."""

import re

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


# --- Phase 2: one parser per study ------------------------------------------------------------
# Every parser takes the study's ENA metadata (one row per run: library_name, sample_title,
# sample_alias, description, isolation_source) and returns SHEET_FIELDS for each row, on the same
# index. kind is environment | inoculum | plate | control | other. Fields a parser cannot read
# from the metadata are left NaN; study-wide values (e.g. a single temperature) are filled from
# config/census_studies.yaml by the calling script. Each parser's naming rules were read from the
# paper and its supplement (sources in notes/sources.md).

SHEET_FIELDS = ["source_sample", "kind", "replicate", "medium", "oxygen", "temperature_c",
                "condition_detail"]
OXYGEN_WORDS = {"AE": "aerobic", "AN": "anaerobic", "MI": "microaerophilic",
                "Aer": "aerobic", "Anaer": "anaerobic"}


def _sheet(meta: pd.DataFrame, row_fn) -> pd.DataFrame:
    """Apply a per-run parser to every row and return SHEET_FIELDS on the input's index."""
    rows = [row_fn(r) for r in meta.itertuples(index=False)]
    out = pd.DataFrame([{f: r.get(f) for f in SHEET_FIELDS} for r in rows], index=meta.index)
    out["kind"] = out["kind"].fillna("other")
    out["replicate"] = pd.to_numeric(out["replicate"])
    out["temperature_c"] = pd.to_numeric(out["temperature_c"])
    return out


def _str(value) -> str:
    """A metadata cell as a stripped string ('' for missing)."""
    return "" if pd.isna(value) else str(value).strip()


# Zheng et al. 2024, both habitats. Sludge sources include digits and a tank suffix
# (XM1, XM2_AE, XM2_AN); media codes differ by habitat (Supplementary Data 1 annotation).
ZHENG2_PATTERN = (r"^(?P<source>[A-Z]+\d?(?:_A[EN])?)_(?P<rep>\d+)_(?P<oxygen>[A-Za-z]+)_"
                  r"(?P<medium>[A-Za-z0-9_]+)$")
ZHENG2_MEDIUM = {
    "soil": {"TSA": "1% TSA", "inSITU": "soil extract"},
    "sludge": {"R2A": "R2A", "inSITU": "extraction medium (filtered wastewater)",
               "influent_inSITU": "extraction medium (filtered influent)",
               "Artificial": "artificial wastewater", "Vb12": "R2A + vitamin B12"},
}
# EBPR (lab reactor) is left out: ENA has no original-sludge (ORI) amplicon library for it, and
# its plate labels mix influent and sludge inocula that the supplement does not explain.
ZHENG2_SOURCES = {"soil": ("TA", "SM", "HN", "YN", "NX", "GD"),
                  "sludge": ("XM1", "XM2", "SZ", "NJ", "SH")}


def parse_zheng2024(meta: pd.DataFrame, habitat: str) -> pd.DataFrame:
    """Zheng library names for the soil or sludge rows: {source}_{rep}_{oxygen}_{medium}.

    ORI_ORI is the original sample, FIL_FIL the plating inoculum, anything else one plate.
    """
    media = ZHENG2_MEDIUM[habitat]

    def row(r):
        m = re.match(ZHENG2_PATTERN, _str(r.library_name))
        if not m or m["source"].split("_")[0] not in ZHENG2_SOURCES[habitat]:
            return {}
        # XM2_AE / XM2_AN are one sludge sample (Supplementary Data 1) whose libraries carry the
        # oxygen of the plate set they were sequenced with; the label is kept as detail.
        base, _, label = m["source"].partition("_")
        out = {"source_sample": base, "replicate": m["rep"], "condition_detail": label or None}
        if m["medium"] == "ORI":
            return {**out, "kind": "environment"}
        if m["medium"] == "FIL":
            return {**out, "kind": "inoculum"}
        return {**out, "kind": "plate", "medium": media.get(m["medium"], m["medium"]),
                "oxygen": OXYGEN_WORDS.get(m["oxygen"])}
    return _sheet(meta, row)


# Krüger et al. 2019 clay till. The compartment is in isolation_source; the sample type is in
# description. PSM samples labelled with a distance (11.5-25 mm) are dispersal sectors, not
# whole plates; "NH" PSM samples carry no distance and are taken to be full plates (D-021).
KRUGER_COMPARTMENTS = [("plow layer", "plow_layer"), ("biopore", "biopore"),
                        ("matrix 80-120", "matrix_80_120"), ("matrix 300-350", "matrix_300_350"),
                        ("red fracture", "red_fracture")]


def kruger_compartment(isolation_source: str) -> str | None:
    """Map an isolation_source text to one of the five soil/sediment compartments."""
    text = isolation_source.lower()
    return next((code for key, code in KRUGER_COMPARTMENTS if key in text), None)


def parse_kruger2019(meta: pd.DataFrame) -> pd.DataFrame:
    """Clay-till compartments: soil DNA, Nycodenz extracts, reference plates and PSM plates."""
    def row(r):
        desc = _str(r.description).removeprefix("type of sample:")
        out = {"source_sample": kruger_compartment(_str(r.isolation_source))}
        if desc.startswith("soil sample"):
            return {**out, "kind": "environment"}
        if desc == "Nycodenz extraction":
            return {**out, "kind": "inoculum"}
        if desc == "reference plate":
            return {**out, "kind": "plate", "medium": "25% R2A",
                    "condition_detail": "reference plate"}
        if desc == "PSM_incubation:NH":
            return {**out, "kind": "plate", "medium": "25% R2A",
                    "condition_detail": "PSM full plate"}
        return {**out, "kind": "other", "condition_detail": desc}
    return _sheet(meta, row)


# Li et al. 2023 desert soil: BS/RS are the two source soils; plates are
# {B|R}{L|M|H}{medium number}, L/M/H = 15/30/45 C (sample titles BCL ... RCH).
LI2023_SOURCE = {"B": "bulk", "R": "rhizosphere"}
LI2023_TEMPERATURE = {"L": 15, "M": 30, "H": 45}
# Medium abbreviations M1-M60 in order, from the list above Supplementary Table 8.
LI2023_MEDIA = ("BEN BHI BRA CSA CHIA CZA GN1 GYM GLA GGA HA HVA KB LNMS LB MCA MYP MIC MM PYG "
              "PYS PPS PDA RA RH RZ1 RZ2 SA SSM SSTY CMCA SLA SEAB SEAR SCA SGCY SNA STA STBG T5 "
              "TP TGY TCN YGCP YMG 1xTSA 0.1xTSA 0.01xTSA 1xNA 0.1xNA 0.01xNA 1xAIA 0.1xAIA "
              "0.01xAIA 1xR2A 0.1xR2A 0.01xR2A 1xMA 0.1xMA 0.01xMA").split()


def parse_li2023(meta: pd.DataFrame) -> pd.DataFrame:
    """Desert soil culturomics: 2 source soils x 60 media x 3 temperatures."""
    def row(r):
        name = _str(r.library_name)
        if name in ("BS", "RS"):
            return {"source_sample": LI2023_SOURCE[name[0]], "kind": "environment"}
        m = re.match(r"^([BR])([LMH])(\d+)$", name)
        if not m:
            return {}
        return {"source_sample": LI2023_SOURCE[m[1]], "kind": "plate",
                "medium": f"M{m[3]} ({LI2023_MEDIA[int(m[3]) - 1]})",
                "temperature_c": LI2023_TEMPERATURE[m[2]]}
    return _sheet(meta, row)


def li2023_incubation_days(table9_text: str) -> dict[str, int]:
    """Plate library name -> incubation days, from the text of Li et al. Supplementary
    Table 9 (rows: medium no., abbreviation, then six (library, days) pairs)."""
    days = {}
    for line in table9_text.splitlines():
        if re.match(r"^\s*M\d+\s", line):
            for lib, d in re.findall(r"\b([BR][LMH]\d+)\s+(\d+)\b", line):
                days[lib] = int(d)
    return days


# Verhoeven et al. 2023: two experiments, each with its own sludge sample, told apart by the
# sequencing batch in library_name (MQ210129 oxic; MQ210611 mostly anoxic). The medium of each
# numbered plate is not published, so medium is the plate code.
VERHOEVEN23_BATCH = {"MQ210129": "AS_exp1", "MQ210611": "AS_exp2"}


def parse_verhoeven2023(meta: pd.DataFrame) -> pd.DataFrame:
    """Activated sludge single-cell plating: sludge, strained suspension, plates, controls."""
    def row(r):
        title = _str(r.sample_title)
        source = VERHOEVEN23_BATCH.get(_str(r.library_name)[:8])
        if m := re.match(r"^MV-PHN-(?:AS-|PL-AS)([A-D1-3])$", title):
            return {"source_sample": source, "kind": "environment", "replicate": _rep(m[1])}
        if m := re.match(r"^MV-PHN-(?:AS-str-|PL-STR)([A-D1-3])$", title):
            return {"source_sample": source, "kind": "inoculum", "replicate": _rep(m[1])}
        if re.match(r"^MV-PHN-PL-(?:NEG|POS)$", title):
            return {"source_sample": source, "kind": "control", "condition_detail": title[-3:]}
        if m := re.match(r"^MV-PHN-PL-(\d+)$", title):
            return {"source_sample": source, "kind": "plate", "medium": "unresolved",
                    "oxygen": "aerobic", "condition_detail": f"plate {m[1]}"}
        if m := re.match(r"^MV-PHN-PL(\d+)-(Aer|Anaer)-([A-D])$", title):
            return {"source_sample": source, "kind": "plate", "medium": f"PL{m[1]}",
                    "oxygen": OXYGEN_WORDS[m[2]], "replicate": _rep(m[3])}
        return {"source_sample": source, "kind": "other", "condition_detail": title}
    return _sheet(meta, row)


def _rep(code: str) -> int:
    """Replicate letters A-D or digits as 1-based integers."""
    return int(code) if code.isdigit() else "ABCD".index(code) + 1


# Villanelo 2026: RAW_n = source sludge; CTRL_n = no antibiotic; {ABX}_{LC|WC|HC}_n = one
# antibiotic at one of three concentrations (codes kept as written).
VILLANELO26_ABX = ("AMP", "CAR", "CHL", "CIP", "ERI", "GEN", "KAN", "NEO", "RIF", "STR", "TET")


def parse_villanelo2026(meta: pd.DataFrame) -> pd.DataFrame:
    """Activated sludge on ASF agarose with and without antibiotics."""
    def row(r):
        name = _str(r.library_name)
        if m := re.match(r"^RAW_(\d+)$", name):
            return {"source_sample": "raw_AS", "kind": "environment", "replicate": m[1]}
        if m := re.match(r"^CTRL_(\d+)$", name):
            return {"source_sample": "raw_AS", "kind": "plate", "medium": "ASF agarose",
                    "replicate": m[1]}
        if (m := re.match(r"^([A-Z]{3})_([LWH]C)_(\d+)$", name)) and m[1] in VILLANELO26_ABX:
            return {"source_sample": "raw_AS", "kind": "plate",
                    "medium": f"ASF agarose + {m[1]}", "condition_detail": m[2],
                    "replicate": m[3]}
        if m := re.match(r"^PCR(NEG|POS)\d$", name):
            return {"source_sample": "raw_AS", "kind": "control", "condition_detail": m[1]}
        return {}
    return _sheet(meta, row)


# Indraningrat et al. 2019 sponges: environment libraries Indra_Curacao_{AS|XM}{n}; plate
# sample titles {medium}.{AS|XM}{n}, so the individual sponge is the shared AS/XM number.
# Media codes from the paper (Table 1). "Co"/"Cryo" is not one of the six media: plates were
# inoculated from cryopreserved (glycerol) sponge homogenate, so those libraries are most likely
# the inoculum despite their ENA description; being ambiguous, they are not used (D-022).
SPONGE_MEDIUM = {"MA": "marine agar 1/10", "M3": "M3 agar", "OL": "OLIGO agar", "GP": "GP agar",
                 "Mu": "mucin agar", "CR": "Crenarchaeota agar"}


def parse_indraningrat2019(meta: pd.DataFrame) -> pd.DataFrame:
    """Sponge individuals and the plate scrapings made from each one."""
    def row(r):
        if m := re.match(r"^Indra_Curacao_((?:AS|XM)\d+)$", _str(r.library_name)):
            return {"source_sample": m[1], "kind": "environment"}
        if (m := re.match(r"^([A-Za-z0-9]+)\.((?:AS|XM)\d+)$", _str(r.sample_title))) \
                and "scraping" in _str(r.description) and m[1] in SPONGE_MEDIUM:
            return {"source_sample": m[2], "kind": "plate", "medium": SPONGE_MEDIUM[m[1]]}
        if m and "sponge_inoculum" in _str(r.description):
            return {"source_sample": m[2], "kind": "inoculum", "condition_detail": m[1]}
        return {}
    return _sheet(meta, row)


# Oak 2025: 16S and ITS libraries share the project; only 16S is used. Tissue is the source.
OAK_TISSUE = {"leaf": "leaf", "stem": "stem", "root": "rhizosphere", "rhizosphere": "rhizosphere"}


def parse_ordonez2025(meta: pd.DataFrame) -> pd.DataFrame:
    """Oak tissue microbiota and pooled agar-plate isolates, 16S only."""
    def row(r):
        alias = _str(r.sample_alias)
        if m := re.match(r"^16S_Ale_(Leaf|Stem|Root)_(\d+)$", alias):
            return {"source_sample": OAK_TISSUE[m[1].lower()], "kind": "environment",
                    "replicate": m[2]}
        if m := re.match(r"^16s-isolation_(\d+)$", alias):
            tissue = _str(r.description).split(" ")[0].lower()
            return {"source_sample": OAK_TISSUE.get(tissue), "kind": "plate",
                    "medium": "unresolved", "condition_detail": f"isolation {m[1]}"}
        return {}
    return _sheet(meta, row)


# Ortiz-Severín 2024 seawater: aliases {site}-A = seawater, {site}-B = culturable pool of 7 media;
# D-I are subcultures of the pool (D-020) and metaG libraries are shotgun.
def parse_ortizseverin2024(meta: pd.DataFrame) -> pd.DataFrame:
    """Two seawater sites with direct community and culturable-pool libraries."""
    def row(r):
        m = re.match(r"^(PSh|AFS)-([A-Z])$", _str(r.sample_alias))
        if not m:
            return {}
        rep = re.search(r"(\d+)$", _str(r.library_name))
        out = {"source_sample": m[1], "replicate": rep[1] if rep else None}
        if m[2] == "A":
            return {**out, "kind": "environment"}
        if m[2] == "B":
            return {**out, "kind": "plate", "medium": "pool of 7 media"}
        return {**out, "kind": "other", "condition_detail": "subculture of the pool"}
    return _sheet(meta, row)


# Yu et al. 2020 broiler carcasses: {area}_{carcass}_SH = skin homogenate;
# {area}_{carcass}_{BA|TSA|PCA}{7|30} = plate wash at 7 or 30 C.
BROILER_MEDIUM = {"BA": "blood agar", "TSA": "TSA", "PCA": "plate count agar"}


def parse_yu2020(meta: pd.DataFrame) -> pd.DataFrame:
    """Skin homogenates and plate washes, linked per carcass and skin area."""
    def row(r):
        m = re.match(r"^(Neck|Breast|Back)_([A-E])_(SH|BA|TSA|PCA)(\d*)$", _str(r.library_name))
        if not m:
            return {}
        source = f"{m[2]}_{m[1].lower()}"
        if m[3] == "SH":
            return {"source_sample": source, "kind": "environment"}
        return {"source_sample": source, "kind": "plate", "medium": BROILER_MEDIUM[m[3]],
                "oxygen": "aerobic", "temperature_c": m[4]}
    return _sheet(meta, row)


# Hinsu et al. 2021 groundnut: the five rhizosphere samples of each farm (R31-R35, R41-R45)
# were pooled and plated; plates R3-{medium}, R4-{medium}. Other soil libraries in the project
# belong to a different study and are not linked. Temperatures (27 and 37 C) were pooled.
def parse_hinsu2021(meta: pd.DataFrame) -> pd.DataFrame:
    """Groundnut rhizosphere: per-farm soil replicates and per-medium pooled plates."""
    def row(r):
        name = _str(r.library_name)
        if m := re.match(r"^R([34])([1-5])$", name):
            return {"source_sample": f"R{m[1]}", "kind": "environment", "replicate": m[2]}
        if m := re.match(r"^R([34])-([A-Z0-9]+)$", name):
            medium = re.search(r"on media (.+)$", _str(r.sample_title))
            return {"source_sample": f"R{m[1]}", "kind": "plate",
                    "medium": medium[1] if medium else m[2],
                    "condition_detail": "pooled 27 and 37 C"}
        return {}
    return _sheet(meta, row)


# Zehavi et al. 2018 rumen: one cow's glycerol-stored rumen sample ("Rumen inoculum source
# Dilution factor 0") was diluted 10-1 to 10-6 ("Rumen inoculum Dilution factor 10-k") and each
# dilution plated in duplicate on M10 with or without clarified rumen fluid. The 37 "Bovine
# rumen sample" libraries are other cows and are not linked (D-028).
RUMEN_MEDIUM = {"with": "M10 + clarified rumen fluid", "without": "M10 defined"}


def parse_zehavi2018(meta: pd.DataFrame) -> pd.DataFrame:
    """The plated cow's rumen sample, its dilutions, and plates per medium and dilution."""
    def row(r):
        title = _str(r.sample_title)
        if title == "Rumen inoculum source Dilution factor 0":
            return {"source_sample": "plated_cow", "kind": "environment"}
        if m := re.match(r"^Rumen inoculum Dilution factor (10-\d)$", title):
            return {"source_sample": "plated_cow", "kind": "inoculum",
                    "condition_detail": f"dilution {m[1]}"}
        if m := re.match(r"^Agar plate media M10 (with|without) rumen fluid Dilution factor "
                         r"(10-\d) rep(\d)$", title):
            return {"source_sample": "plated_cow", "kind": "plate",
                    "medium": RUMEN_MEDIUM[m[1]], "condition_detail": f"dilution {m[2]}",
                    "replicate": m[3]}
        return {}
    return _sheet(meta, row)


PARSERS = {
    "zheng2024_soil": lambda meta: parse_zheng2024(meta, "soil"),
    "zheng2024_sludge": lambda meta: parse_zheng2024(meta, "sludge"),
    "kruger2019_claytill": parse_kruger2019,
    "li2023_desertsoil": parse_li2023,
    "verhoeven2023_sludge": parse_verhoeven2023,
    "villanelo2026_sludge_antibiotics": parse_villanelo2026,
    "indraningrat2019_sponges": parse_indraningrat2019,
    "ordonez2025_oak": parse_ordonez2025,
    "ortizseverin2024_seawater": parse_ortizseverin2024,
    "yu2020_broiler": parse_yu2020,
    "hinsu2021_groundnut": parse_hinsu2021,
    "zehavi2018_rumen": parse_zehavi2018,
}


def link_summary(sheet: pd.DataFrame) -> pd.DataFrame:
    """Per study: environment and plate runs, and how many plates have a linked source.

    A plate is linked when its study has at least one environment run with the same
    source_sample. Returns one row per study.
    """
    env = sheet[sheet["kind"] == "environment"]
    plates = sheet[sheet["kind"] == "plate"]
    env_keys = set(zip(env["study_id"], env["source_sample"]))
    linked = [k in env_keys for k in zip(plates["study_id"], plates["source_sample"])]
    plates = plates.assign(linked=linked)
    out = pd.DataFrame({
        "environment_runs": env.groupby("study_id").size(),
        "plate_runs": plates.groupby("study_id").size(),
        "linked_plate_runs": plates.groupby("study_id")["linked"].sum(),
        "linked_sources": plates[plates["linked"]].groupby("study_id")["source_sample"].nunique(),
    }).fillna(0).astype(int)
    return out.rename_axis("study_id").reset_index()
