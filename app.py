#!/usr/bin/env python3
"""
Streamlit webapp for interactive 3D visualization of planetary magnetic fields.
Run locally:  streamlit run app.py
"""

import numpy as np
import streamlit as st
import pyvista as pv
from stpyvista import stpyvista

from planetmagfields import Planet, get_models
from planetmagfields.utils import planetlist
from planetmagfields.lib3d import plot_surface, render_field_lines

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="PlanetMagFields 3D",
    page_icon="🪐",
    layout="wide",
    initial_sidebar_state="auto"
)

# ---------------------------------------------------------------------------
# Cached helpers
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def fetch_models(planet_name: str) -> list[str]:
    return sorted(get_models(planet_name).tolist())

@st.cache_data(show_spinner="Computing field...")
def load_planet(name: str, model: str, r: float, units: str, nphi: int) -> Planet:
    return Planet(name=name, model=model, r=r, units=units, nphi=nphi, info=False)

def plot_field_surf(res: int, fieldlines: bool, bgcolor: str, rout: float, colorlim_fac: float, units: str, cmap: str):
    st.session_state["render_params"] = dict(
        res=res, fieldlines=fieldlines, bgcolor=bgcolor, rout=rout, colorlim_fac=colorlim_fac, units=units, cmap=cmap
    )

with st.sidebar:
    st.title("Controls")

    planet_name = st.selectbox(
        "Planet",
        options=planetlist,
        index=planetlist.index("earth"),
        format_func=str.capitalize,
    )

    models = fetch_models(planet_name)
    model = st.selectbox("Model", options=models, index=len(models) - 1)

    r = st.slider(
        "Radius  (r / R_planet)",
        min_value=0.5, max_value=5.0, value=1.0, step=0.05,
    )

    units = st.selectbox("Field unit", ["muT", "nT", "Gauss"])

    resolution = st.select_slider(
        "Resolution  (nphi)",
        options=[64, 128, 256, 512],
        value=128,
    )

    bgcolor = st.color_picker("Background color", value="#000000")

    fieldlines = st.checkbox("Show field lines", value=False)

    routmax = st.slider(
        "Outer radius  (rout / R_planet)",
        min_value=1.0, max_value=10.0, value=2.0, step=0.1,
    )

    nr = st.slider(
        "Radial resolution  (nr) for field lines",
        min_value=16, max_value=256, value=32, step=4,
    )

    rout = np.linspace(1, routmax, nr)

    colorlim_fac = st.slider(
        "Color limit factor  (clim_fac)",
        min_value=0.01, max_value=1.0, value=0.05, step=0.01,
    )

    # Only diverging colormaps
    cmap = st.selectbox("Colormap", ["seismic", "RdBu_r", "Spectral_r"])

    st.button("Render", on_click=plot_field_surf,args=[resolution, fieldlines,
                                                       bgcolor, rout,
                                                       colorlim_fac, units, cmap])

st.title("Planetary Magnetic Fields — 3D Explorer")
planet = load_planet(planet_name, model, r, units, resolution)

# Metrics row
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Planet",      planet.name.capitalize())
c2.metric("Model",       planet.model)
c3.metric("l_max",       planet.lmax)
c4.metric("Dipole tilt", f"{planet.dipTheta:.2f}°")
c5.metric("r / Rp",      f"{r:.2f}")

if "render_params" in st.session_state:
    p = st.session_state["render_params"]
    with st.spinner("Plotting field...", show_time=True):
        if not p["fieldlines"]:
            pl, _ = plot_surface(
                planet.theta, planet.phi, planet.Br,
                clim_fac=p["colorlim_fac"], bgcolor=p["bgcolor"],
                cmap=p["cmap"]
            )
        else:
            pl = render_field_lines(
                planet.name, planet.glm.copy(), planet.hlm.copy(), planet.idx, planet.lmax, planet.mmax,
                1, p["rout"], nphi=p["res"], surf=True,
                clim_fac=p["colorlim_fac"], bgcolor=p["bgcolor"],
                units=planet.units, cmap=p["cmap"]
            )
        stpyvista(pl, key=f"sphere_{p['res']}", backend="panel")
