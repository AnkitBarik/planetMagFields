# app.py
import streamlit as st

st.set_page_config(
    page_title="planetMagFields App",
    page_icon="🚀",
    layout="wide"
)

pg = st.navigation([
    st.Page("pages/3d_explorer.py", title="3D Explorer", icon="🌍"),
    st.Page("pages/orbit_trajectory.py", title="Orbit Trajectory", icon="🛰️"),
])
pg.run()