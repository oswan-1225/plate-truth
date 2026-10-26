"""Upload one sample's ASV table and get per-taxon predictions and the shortcut comparison."""

import io
import tempfile
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st
from common import CAVEATS, bundle, example, pct

from plate_truth.app_model import community_summary, family_composition, predict, prepare_table

b = bundle()
s = b["settings"]
K = s["conditions_k"]
REFERENCES = Path(__file__).resolve().parents[2] / "data" / "interim" / "reference" / "regions"

st.title("Your community")
st.write(
    "Upload one sample's ASV table to see how likely each taxon is to be enriched on plates, "
    "and how the usual shortcut (\"a cultured relative at "
    f"{s['shortcut_percent']:g}% identity or closer means culturable\") compares with what one "
    "plating is predicted to recover.")

with st.expander("Table format", icon=":material/table:"):
    st.markdown(
        "One row per ASV of **one sample**, including every prokaryotic ASV (remove chloroplast "
        "and mitochondrial sequences first). Columns:\n\n"
        "- `asv_id` and `count` (required)\n"
        "- `identity`: percent identity to the nearest cultured relative, computed over the "
        "amplified region with all gaps counted (VSEARCH `--iddef 1`); leave empty when there "
        "is no hit at 75% or above\n"
        "- `phylum` and `family`: SILVA 138 names (optional; without `phylum` the phylum effect "
        "is set to its average)\n\n"
        "Instead of `identity` you can upload the ASV sequences as FASTA when the app runs on a "
        "computer with the project pipeline installed (VSEARCH and the cultured reference).")

source = st.segmented_control("Data", ["Upload a table", "Use the example"],
                              default="Use the example", key="source")
raw, fasta = None, None
if source == "Upload a table":
    up = st.file_uploader("ASV table (CSV or TSV)", type=["csv", "tsv", "txt"])
    fasta = st.file_uploader("ASV sequences (FASTA, optional)", type=["fasta", "fa", "fna"])
    if up is not None:
        sep = "\t" if up.name.endswith((".tsv", ".txt")) else ","
        raw = pd.read_csv(up, sep=sep)
else:
    raw = example()
    st.caption("Example: one soil sample (ENA run ERR11027514) from a public MGnify amplicon "
               "survey of the V3-V4 region, with identities computed by the project pipeline.")

if raw is None:
    st.info("Upload a table to start, or choose the example.", icon=":material/upload:")
    st.stop()


def identities_from_fasta(fasta_bytes: bytes, reference: Path) -> pd.Series:
    """Identity to the nearest cultured relative for each uploaded sequence (local only)."""
    from plate_truth.identity import nearest_relative_table, read_hits, search
    from plate_truth.reference import read_fasta
    with tempfile.TemporaryDirectory() as tmp:
        q = Path(tmp) / "asvs.fasta"
        q.write_bytes(fasta_bytes)
        hits = Path(tmp) / "hits.tsv"
        search(q, reference, hits, 1, 0.75, 0, 0, top_hits_only=True)
        ids = [h.split(";")[0].split()[0] for h, _ in read_fasta(q)]
        return nearest_relative_table(ids, read_hits(hits), id_cols=("id1",))["identity_id1"]


if fasta is not None and "identity" not in {c.lower() for c in raw.columns}:
    refs = sorted(p for p in REFERENCES.glob("cultured_*.fasta") if "unfiltered" not in p.name)
    if not refs:
        st.error("Computing identity needs the cultured reference built by the project "
                 "pipeline (scripts 03 to 15), which is not available here. Add an `identity` "
                 "column instead.", icon=":material/error:")
        st.stop()
    ref = st.selectbox("Cultured reference cut to your primer pair", refs,
                       format_func=lambda p: p.stem.removeprefix("cultured_"))
    with st.spinner("Searching the cultured reference with VSEARCH"):
        try:
            idmap = identities_from_fasta(fasta.getvalue(), ref)
        except FileNotFoundError as err:
            st.error(f"VSEARCH is not available: {err}", icon=":material/error:")
            st.stop()
    raw = raw.copy()
    raw["identity"] = raw[next(c for c in raw.columns if c.lower() == "asv_id")].map(idmap)
    if idmap.isna().mean() > 0.5:
        st.warning(f"{idmap.isna().mean():.0%} of the sequences have no cultured relative at "
                   "75% identity or above. This usually means the reference was cut to a "
                   "different region than your amplicons; check the primer pair.",
                   icon=":material/warning:")

try:
    table, notes = prepare_table(raw, s["present_min_reads"], s["no_hit_value"])
except ValueError as err:
    st.error(f"The table could not be read: {err}", icon=":material/error:")
    st.stop()

total_reads = pd.to_numeric(raw[next(c for c in raw.columns if c.lower() == "count")],
                            errors="coerce").sum()
if total_reads < s["min_prokaryotic_reads"]:
    notes.append(f"The sample has {total_reads:,.0f} reads, fewer than the "
                 f"{s['min_prokaryotic_reads']:,} required in the training and survey data.")
for n in notes:
    st.warning(n, icon=":material/warning:")


@st.cache_data(max_entries=8, show_spinner="Predicting")
def run_predictions(t: pd.DataFrame, phylum_given: bool) -> pd.DataFrame:
    t = t.copy()
    t.attrs["phylum_given"] = phylum_given
    return predict(t, bundle())


pred = run_predictions(table, table.attrs["phylum_given"])
summ = community_summary(pred, K, s["shortcut_percent"])

st.subheader("The sample as a whole")
with st.container(horizontal=True):
    st.metric("Shortcut: reads with a cultured relative", pct(summ["shortcut"]), border=True,
              help=f"Share of the sample's reads in ASVs with identity of at least "
                   f"{s['shortcut_percent']:g}% to a cultured strain.")
    st.metric("Predicted, one plate", f"{pct(summ['one_plate_p10'])} to "
              f"{pct(summ['one_plate_p90'])}", border=True,
              help="Share of reads in taxa predicted to be enriched on one plate, from a "
                   "low-recovery study (10th percentile) to a high-recovery one (90th).")
    st.metric(f"Predicted, any of {K} plates (high-recovery study)",
              pct(summ[f"any_of_{K}_p90"]), border=True,
              help="A generous reference: plates are assumed independent beyond their shared "
                   "study and sample.")

if summ["shortcut"] > summ["one_plate_p90"]:
    st.write(
        f"The shortcut counts **{pct(summ['shortcut'], 0)}** of this sample as culturable. "
        f"One plating is predicted to recover at most about **{pct(summ['one_plate_p90'], 0)}** "
        "even in a study at the high end of what the training studies achieved, so the "
        f"shortcut overstates it about {summ['shortcut'] / summ['one_plate_p90']:.0f}-fold.")
else:
    st.write("Here the shortcut does not exceed what one plating is predicted to recover in a "
             "high-recovery study.")

st.subheader("Each taxon")
plot = pred.assign(
    relative=lambda d: d["has_cultured_relative"].map(
        {True: f"cultured relative ≥{s['shortcut_percent']:g}%",
         False: f"no cultured relative ≥{s['shortcut_percent']:g}%"}),
    identity_plot=lambda d: d["identity"].fillna(s["no_hit_value"]),
    percent=lambda d: 100 * d["predicted_one_plate_p50"])
chart = alt.Chart(plot).mark_circle(opacity=0.6).encode(
    x=alt.X("identity_plot:Q", title="Identity to nearest cultured relative (%; no hit at 75)",
            scale=alt.Scale(domain=[s["no_hit_value"], 100])),
    y=alt.Y("percent:Q", title="Predicted enrichment on one plate, median study (%)"),
    size=alt.Size("a:Q", title="Relative abundance", legend=None),
    color=alt.Color("relative:N", title=None, legend=alt.Legend(orient="top")),
    tooltip=["asv_id", "phylum", "family", alt.Tooltip("identity:Q", format=".1f"),
             alt.Tooltip("a:Q", title="relative abundance", format=".4f"),
             alt.Tooltip("percent:Q", title="predicted (%)", format=".2f")])
st.altair_chart(chart, width="stretch")
st.caption("Each point is an ASV; larger points are more abundant. Predictions average over "
           "the plating conditions of the training studies.")

comp = family_composition(pred, b["family_calls"])
with st.container(border=True):
    st.markdown("**Families that grew less often than relatedness predicts**")
    st.write(
        f"{pct(comp['under'])} of reads are in families that, in the training studies, were "
        "enriched on plates systematically less often than their relatedness to cultured strains "
        f"predicts; {pct(comp['over'])} are in families recovered more often; "
        f"{pct(comp['not evaluated'])} are in families the training data could not evaluate.")

cols = {"asv_id": "ASV", "count": "reads", "a": "relative abundance", "identity": "identity (%)",
        "phylum": "phylum", "family": "family",
        "has_cultured_relative": "shortcut: has cultured relative",
        "predicted_one_plate_p10": "predicted, one plate, low study",
        "predicted_one_plate_p50": "predicted, one plate, median study",
        "predicted_one_plate_p50_lo": "median study, 95% lower",
        "predicted_one_plate_p50_hi": "median study, 95% upper",
        "predicted_one_plate_p90": "predicted, one plate, high study",
        f"predicted_any_of_{K}_p90": f"predicted, any of {K} plates, high study"}
shown = pred[list(cols)].rename(columns=cols).sort_values("reads", ascending=False)
st.dataframe(shown, hide_index=True, height=360)
buf = io.StringIO()
pred[list(cols)].to_csv(buf, index=False)
st.download_button("Download predictions (CSV)", buf.getvalue(), "plate_truth_predictions.csv",
                   "text/csv", icon=":material/download:")

with st.expander("What these numbers mean, and what they do not", icon=":material/help:"):
    st.markdown(CAVEATS)
