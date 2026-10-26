"""Phase 5 results by habitat: the shortcut against the predicted range of recovery."""

import altair as alt
import pandas as pd
import streamlit as st
from common import CAVEATS, habitats, pct

NAMES = {"soil": "Soil", "activated_sludge": "Activated sludge", "plant": "Plant",
         "marine_water": "Marine water", "marine_sediment": "Marine sediment",
         "freshwater": "Freshwater", "human_gut": "Human gut"}

st.title("Public surveys")
st.write(
    "The model applied to public MGnify surveys, selected by rules registered before any "
    "sample was processed. Because a new study's overall level cannot be predicted, each "
    "habitat gets a range, from a low-recovery to a high-recovery study, and the test is "
    "whether the shortcut exceeds even the top of that range.")

arm = st.segmented_control("Survey type", ["amplicon", "shotgun"], default="amplicon",
                           format_func=lambda a: {"amplicon": "16S amplicon",
                                                  "shotgun": "Shotgun (16S fragments, V4)"}[a],
                           key="arm", required=True)
h = habitats(arm).copy()
h["name"] = h["habitat"].map(NAMES) + h["training_habitat"].map({True: "", False: " *"})

bars = pd.DataFrame({
    "habitat": h["name"].repeat(2).to_numpy(),
    "what": ["any of 8 plates", "one plate"] * len(h),
    "low": h[["r8_p10", "r1_p10"]].to_numpy().ravel(),
    "high": h[["r8_p90", "r1_p90"]].to_numpy().ravel()})
order = list(h["name"])
y = alt.Y("habitat:N", sort=order, title=None)
range_chart = alt.Chart(bars).mark_bar(height=14).encode(
    y=y, x=alt.X("low:Q", title="Share of reads", scale=alt.Scale(domain=[0, 1])), x2="high:Q",
    color=alt.Color("what:N", title="Predicted recovery (10th-90th percentile study)",
                    scale=alt.Scale(domain=["one plate", "any of 8 plates"],
                                    range=["#2a78d6", "#a9c9ef"]),
                    legend=alt.Legend(orient="top")),
    tooltip=["habitat", "what", alt.Tooltip("low:Q", format=".3f"),
             alt.Tooltip("high:Q", format=".3f")])
shortcut = alt.Chart(h).mark_point(shape="diamond", size=140, filled=True, color="#e34948"
                                   ).encode(
    y=alt.Y("name:N", sort=order), x="shortcut:Q",
    tooltip=[alt.Tooltip("name:N", title="habitat"),
             alt.Tooltip("shortcut:Q", title="shortcut share", format=".2f")])
st.altair_chart(range_chart + shortcut, width="stretch")
st.caption("Red diamond: share of reads with a cultured relative at 94.5% identity or closer "
           "(the shortcut). * = no training study from this habitat (extrapolation).")

table = pd.DataFrame({
    "habitat": h["name"],
    "studies / samples": h["studies"].astype(str) + " / " + h["samples"].astype(str),
    "shortcut": [f"{pct(r.shortcut, 0)} ({pct(r.shortcut_lo, 0)}-{pct(r.shortcut_hi, 0)})"
                 for r in h.itertuples()],
    "one plate, 10th-90th": [f"{pct(r.r1_p10)}-{pct(r.r1_p90)}" for r in h.itertuples()],
    "shortcut minus one-plate upper end": [
        f"{r.d_mean_r1:.2f} ({r.d_lo_r1:.2f} to {r.d_hi_r1:.2f})" for r in h.itertuples()],
    "exceeds one plate": h["shortcut_exceeds_r1"].map({True: "yes", False: "no"}),
    "exceeds any of 8": h["shortcut_exceeds_r8"].map({True: "yes", False: "no"}),
    "reads in under-recovered families": [pct(v, 0) for v in h["family_under"]],
})
st.dataframe(table, hide_index=True)
st.caption("Intervals: 95% over studies and fixed-effect draws. A habitat with one study has "
           "intervals from the model draws only.")

with st.expander("What these numbers mean, and what they do not", icon=":material/help:"):
    st.markdown(CAVEATS)
