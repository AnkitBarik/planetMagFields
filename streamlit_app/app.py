# app.py
"""
The planetMagFields web app has moved to GitHub Pages. This Streamlit app only
points visitors there, for old links and bookmarks.
"""
import streamlit as st

NEW_URL = "https://ankitbarik.github.io/planetMagFields/app/"

st.set_page_config(page_title="planetMagFields has moved", page_icon="🪐", layout="centered")

st.title("planetMagFields has moved")
st.markdown(
    f"The interactive 3D explorer now runs directly in your browser at "
    f"**[{NEW_URL}]({NEW_URL})**. It is faster, never sleeps, and includes "
    f"everything this app did, including magnetic fields along spacecraft trajectories."
)
st.link_button("Open the new app", NEW_URL, type="primary")
