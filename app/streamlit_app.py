"""plate-truth: what a cultured relative predicts about growth on plates.

Run from the repository root:  streamlit run app/streamlit_app.py
"""

import sys
from pathlib import Path

import streamlit as st

# The app uses the project package (src/plate_truth) without requiring an install, and the
# pages import app/common.py.
APP = Path(__file__).resolve().parent
for path in (APP, APP.parent / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

st.set_page_config(page_title="plate-truth", page_icon=":material/biotech:", layout="wide")

pages = st.navigation([
    st.Page("app_pages/community.py", title="Your community", icon=":material/upload_file:",
            default=True),
    st.Page("app_pages/surveys.py", title="Public surveys", icon=":material/public:"),
    st.Page("app_pages/about.py", title="About the model", icon=":material/info:"),
])
pages.run()
