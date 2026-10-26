"""What the model is, how it was validated, and its limits."""

import altair as alt
import pandas as pd
import streamlit as st
from common import APP, CAVEATS

st.title("About the model")
st.markdown(
    "Estimates of how much of the microbial world can be cultured often count a 16S sequence "
    "as culturable when it closely matches a cultured strain. This project tested that "
    "shortcut against what actually grew on plates in ten published studies that sequenced "
    "both an environmental sample and the communities grown from it.\n\n"
    "**Outcome.** An ASV is *enriched on plates* when it has at least 5 reads on a plate and at "
    "least twice its relative abundance in the source sample. DNA carried over from the "
    "inoculum tracks its source abundance; growth raises it.\n\n"
    "**Model.** A mixed-effects logistic regression of enrichment on identity to the nearest "
    "cultured relative (natural spline), source abundance, plate depth, medium class, "
    "temperature, oxygen and incubation time, with random effects for study, source sample, "
    "plating condition and phylum. It was registered before fitting and validated by "
    "holding out one study at a time.\n\n"
    "**What it can and cannot do.** In a held-out study it ranks taxa well, but the overall "
    "level of growth differed between studies more than tenfold and could not be predicted. "
    "The app therefore always reports a range across the between-study distribution.")

curve = pd.read_csv(APP / "data" / "curve.csv")
curve = curve[curve["no_hit"] == 0].assign(
    **{c: lambda d, c=c: 100 * d[c] for c in ("p", "p_lo", "p_hi", "study_p10", "study_p90")})
base = alt.Chart(curve).encode(x=alt.X("identity:Q", title="Identity to nearest cultured "
                                                            "relative (%)",
                                       scale=alt.Scale(domain=[75, 100])))
band = base.mark_area(opacity=0.2, color="#2a78d6").encode(
    y=alt.Y("p_lo:Q", title="Enriched on plates, new study (%)"), y2="p_hi:Q")
spread = base.mark_line(strokeDash=[4, 3], color="#2a78d6").encode(y="study_p10:Q") + \
    base.mark_line(strokeDash=[4, 3], color="#2a78d6").encode(y="study_p90:Q")
line = base.mark_line(color="#2a78d6", strokeWidth=2.5).encode(
    y="p:Q", tooltip=[alt.Tooltip("identity:Q", format=".2f"),
                      alt.Tooltip("p:Q", title="predicted (%)", format=".2f")])
st.altair_chart(band + spread + line, width="stretch")
st.caption("Probability that an ASV is enriched on plates in a study the model has not seen, "
           "averaged over the analysed communities and conditions. Band: 95% interval; dashed: "
           "10th and 90th percentile study.")

st.subheader("Caveats")
st.markdown(CAVEATS)
st.markdown("Code, analysis plan and data: "
            "[github.com/oswan-1225/plate-truth](https://github.com/oswan-1225/plate-truth).")
