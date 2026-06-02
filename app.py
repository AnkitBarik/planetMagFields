#!/usr/bin/env python3
"""
Streamlit webapp for interactive 3D visualization of planetary magnetic fields.
Run locally:  streamlit run app.py
"""

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy.integrate import solve_ivp

from planetmagfields import Planet, get_models
from planetmagfields.libgauss import getB
from planetmagfields.potextra import get_field_along_path_scipy, get_pol_from_Gauss
from planetmagfields.utils import planetlist

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="PlanetMagFields 3D",
    page_icon="🪐",
    layout="wide",
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


@st.cache_data(show_spinner=False)
def compute_surface_br(
    name: str, model: str, r: float, units: str, nphi: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute Br on a uniform (regular) grid to avoid pole holes in the 3D sphere."""
    planet = load_planet(name, model, r, units, nphi)
    ntheta = nphi // 2
    phi   = np.linspace(0, 2 * np.pi, nphi + 1)
    theta = np.linspace(0, np.pi, ntheta + 1)
    p2D, th2D = np.meshgrid(phi, theta, indexing='ij')
    br = getB(
        planet.lmax, planet.mmax,
        planet.glm, planet.hlm, planet.idx,
        r, p2D, th2D,
        planetname=planet.name,
    )
    br *= planet.unitfac
    return p2D, th2D, br


# ---------------------------------------------------------------------------
# Seeding
# ---------------------------------------------------------------------------

def get_seed_points(
    planet: Planet,
    n_strong: int,
    n_random: int,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Mix top-|Br| surface locations with uniform-random seeds on the sphere."""
    br_abs = np.abs(planet.Br)          # shape (nphi, ntheta)

    # Strongest-field seeds
    flat_idx = np.argsort(br_abs.ravel())[-n_strong:]
    iphi, itheta = np.unravel_index(flat_idx, br_abs.shape)
    strong_theta = planet.theta[itheta]
    strong_phi   = planet.phi[iphi]

    # Uniform random seeds on the sphere
    rng = np.random.default_rng(rng_seed)
    rand_theta = np.arccos(rng.uniform(-1.0, 1.0, n_random))
    rand_phi   = rng.uniform(0.0, 2 * np.pi, n_random)

    return (
        np.concatenate([strong_theta, rand_theta]),
        np.concatenate([strong_phi,   rand_phi]),
    )


# ---------------------------------------------------------------------------
# Field-line tracing  (shtns if available, scipy fallback)
# ---------------------------------------------------------------------------

def compute_field_lines(
    planet: Planet,
    seed_thetas: np.ndarray,
    seed_phis: np.ndarray,
    rmax: float,
) -> list[tuple]:
    """Trace field lines from surface seed points; return (x, y, z) lists."""
    glm, hlm, idx, lmax, mmax = planet.glm, planet.hlm, planet.idx, planet.lmax, planet.mmax

    # --- one-time shtns setup -------------------------------------------
    _sh = None
    _L = _bpolcmb = None
    try:
        import shtns
        bpol = get_pol_from_Gauss(planet.name, glm, hlm, lmax, mmax, idx)
        _sh = shtns.sht(lmax, mmax=lmax, norm=shtns.sht_orthonormal)
        _L = _sh.l * (_sh.l + 1)
        _bpolcmb = _sh.spec_array()
        for l in range(1, lmax + 1):
            for m in range(l + 1):
                _bpolcmb[_sh.idx(l, m)] = bpol[idx[l, m]]
    except ImportError:
        pass  # fall back to scipy

    def _get_field(r, theta, phi):
        if _sh is not None:
            bp_r = _bpolcmb * (1.0 / r) ** _sh.l
            qlm  = bp_r * _L / r ** 2
            slm  = -_sh.l / r ** 2 * bp_r
            tlm  = np.zeros_like(qlm)
            br, bt, bp = _sh.SHqst_to_point(qlm, slm, tlm, np.cos(theta), phi)
            return float(br), float(bt), float(bp)
        br, bt, bp = get_field_along_path_scipy(
            glm, hlm, idx, lmax,
            np.array([r]), np.array([theta]), np.array([phi]),
        )
        return float(br[0]), float(bt[0]), float(bp[0])

    def _ode(t, state, sign):
        r, theta, phi = state
        phi   = phi % (2 * np.pi)
        theta = np.clip(theta, 1e-5, np.pi - 1e-5)
        br, bt, bp = _get_field(r, theta, phi)
        mag  = np.sqrt(br ** 2 + bt ** 2 + bp ** 2) + 1e-12
        sin_t = np.sin(theta)
        return [sign * br / mag,
                sign * bt / (r * mag),
                sign * bp / (r * sin_t * mag)]

    def hit_surface(t, s, sign): return s[0] - 1.0
    def hit_rmax(t, s, sign):    return s[0] - rmax
    hit_surface.terminal = True;  hit_surface.direction = -1
    hit_rmax.terminal    = True;  hit_rmax.direction    =  1

    lines = []
    for theta0, phi0 in zip(seed_thetas, seed_phis):
        fwd = solve_ivp(_ode, [0, 40], [1.01, theta0, phi0], args=(1,),
                        events=[hit_surface, hit_rmax], max_step=0.1)
        bwd = solve_ivp(_ode, [0, 40], [1.01, theta0, phi0], args=(-1,),
                        events=[hit_surface, hit_rmax], max_step=0.1)
        r     = np.concatenate([bwd.y[0, ::-1], fwd.y[0, 1:]])
        theta = np.concatenate([bwd.y[1, ::-1], fwd.y[1, 1:]])
        phi   = np.concatenate([bwd.y[2, ::-1], fwd.y[2, 1:]])
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)
        lines.append((x, y, z))
    return lines


# ---------------------------------------------------------------------------
# 3D figure builder
# ---------------------------------------------------------------------------

# Plotly diverging colorscales (RdBu_r → RdBu + reversescale)
_CMAPS = {
    "RdBu_r":   ("RdBu",      True),
    "RdBu":     ("RdBu",      False),
    "Spectral":  ("Spectral",  False),
    "Spectral_r":("Spectral",  True),
    "RdYlBu_r": ("RdYlBu",    True),
    "RdYlBu":   ("RdYlBu",    False),
    "Picnic":    ("Picnic",    False),
    "Bluered":   ("Bluered",   False),
}


def build_figure(
    planet: Planet,
    p2D: np.ndarray,
    th2D: np.ndarray,
    br: np.ndarray,
    cmap: str,
    vmin: float | None,
    vmax: float | None,
    field_lines: list | None = None,
    line_color: str = "white",
    line_width: int = 2,
) -> go.Figure:
    x = np.sin(th2D) * np.cos(p2D)
    y = np.sin(th2D) * np.sin(p2D)
    z = np.cos(th2D)

    abs_max = float(np.abs(br).max())
    cmin = vmin if vmin is not None else -abs_max
    cmax = vmax if vmax is not None else  abs_max

    colorscale, reversescale = _CMAPS.get(cmap, ("RdBu", True))

    surface = go.Surface(
        x=x, y=y, z=z,
        surfacecolor=br,
        colorscale=colorscale,
        reversescale=reversescale,
        cmin=cmin,
        cmax=cmax,
        colorbar=dict(
            title=dict(text=f"Br ({planet.units})", side="right"),
            thickness=18,
        ),
        lighting=dict(ambient=0.7, diffuse=0.6, specular=0.1),
        hovertemplate=(
            "Br: %{surfacecolor:.3f} " + planet.units + "<extra></extra>"
        ),
    )

    fig = go.Figure(data=[surface])

    if field_lines:
        colors = px.colors.sample_colorscale("plasma", len(field_lines))
        for i, (lx, ly, lz) in enumerate(field_lines):
            fig.add_trace(go.Scatter3d(
                x=lx, y=ly, z=lz,
                mode="lines",
                line=dict(color=line_color if line_color != "rainbow" else colors[i],
                          width=line_width),
                hoverinfo="skip",
                showlegend=False,
            ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor="black",
            aspectmode="data",
            camera=dict(eye=dict(x=1.4, y=1.4, z=0.8)),
        ),
        paper_bgcolor="#0e1117",
        margin=dict(l=0, r=0, t=0, b=0),
        height=650,
    )
    return fig


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
_DEFAULT_MODELS = {
    "earth":    "igrf14",
    "mercury":  "wardinski2019",
    "jupiter":  "jrm33",
    "saturn":   "cassini11+",
    "uranus":   "holme1996",
    "neptune":  "connerny1991",
    "ganymede": "kivelson2002",
}

with st.sidebar:
    st.title("Controls")

    planet_name = st.selectbox(
        "Planet",
        options=planetlist,
        index=planetlist.index("earth"),
        format_func=str.capitalize,
    )

    models = fetch_models(planet_name)
    default_model = _DEFAULT_MODELS.get(planet_name, models[-1])
    default_idx = models.index(default_model) if default_model in models else len(models) - 1
    model = st.selectbox("Model", options=models, index=default_idx)

    r_reset_col, r_col = st.columns([1, 4])
    with r_reset_col:
        st.write("")  # vertical align
        if st.button("↺", key="reset_r", help="Reset radius to 1.0"):
            st.session_state["r_val"] = 1.0
    with r_col:
        r = st.slider(
            "Radius  (r / R_planet)",
            min_value=0.5, max_value=5.0,
            value=1.0, step=0.05,
            key="r_val",
        )

    units = st.selectbox("Field units", ["muT", "nT", "Gauss"])

    res_reset_col, res_col = st.columns([1, 4])
    with res_reset_col:
        st.write("")  # vertical align
        if st.button("↺", key="reset_res", help="Reset resolution to 128"):
            st.session_state["res_val"] = 128
    with res_col:
        resolution = st.select_slider(
            "Resolution  (nphi)",
            options=[64, 128, 256, 512],
            value=128,
            key="res_val",
        )

    st.divider()
    st.subheader("Colour scale")

    cmap = st.selectbox("Colormap", list(_CMAPS.keys()))

    sym = st.checkbox("Symmetric  (auto ±max)", value=True)
    if not sym:
        col1, col2 = st.columns(2)
        vmin = col1.number_input("vmin", value=0.0, format="%.3f")
        vmax = col2.number_input("vmax", value=0.0, format="%.3f")
    else:
        vmin = vmax = None

    st.divider()
    st.subheader("Field lines")

    show_lines = st.checkbox("Show field lines", value=False)
    fl_settings = {}
    if show_lines:
        fl_settings["n_strong"] = 1  # fixed, not exposed
        fl_settings["n_random"] = st.slider(
            "Number of field lines", min_value=0, max_value=30, value=5,
        )
        fl_settings["rmax"] = st.slider(
            "Max radius (Rp)", min_value=1.5, max_value=8.0, value=3.0, step=0.5
        )
        fl_settings["line_color"] = st.selectbox(
            "Line colour", ["white", "rainbow", "yellow", "cyan", "black"],
        )
        fl_settings["line_width"] = st.slider(
            "Line width", min_value=1, max_value=5, value=2
        )

# ---------------------------------------------------------------------------
# Main panel
# ---------------------------------------------------------------------------
st.title("Planetary Magnetic Fields — 3D Explorer")

planet = load_planet(planet_name, model, r, units, resolution)
p2D, th2D, br_surf = compute_surface_br(planet_name, model, r, units, resolution)

st.markdown("""
<style>
[data-testid="stMetricLabel"] { font-size: 1.5rem !important; }
[data-testid="stMetricValue"] { font-size: 1.6rem !important; }
</style>
""", unsafe_allow_html=True)

# Metrics row
c1, c2, c3, c4, c5 = st.columns(5)

def _metric_html(label: str, value: str) -> str:
    return (
        f'<div style="text-align:center">'
        f'<div style="font-size:2.0rem;color:#aaa;margin-bottom:2px">{label}</div>'
        f'<div style="font-size:1.6rem;font-weight:600">{value}</div>'
        f'</div>'
    )

c1.markdown(_metric_html("Planet",      planet.name.capitalize()),        unsafe_allow_html=True)
c2.markdown(_metric_html("Model",       planet.model),                     unsafe_allow_html=True)
c3.markdown(_metric_html("l_max",       str(planet.lmax)),                 unsafe_allow_html=True)
c4.markdown(_metric_html("Dipole tilt", f"{planet.dipTheta:.2f}°"),       unsafe_allow_html=True)
c5.markdown(_metric_html("r / Rp",      f"{r:.2f}"),                       unsafe_allow_html=True)

# Optionally compute and trace field lines
field_lines = None
if show_lines and fl_settings:
    with st.spinner("Tracing field lines..."):
        seed_thetas, seed_phis = get_seed_points(
            planet,
            n_strong=fl_settings["n_strong"],
            n_random=fl_settings["n_random"],
        )
        field_lines = compute_field_lines(
            planet, seed_thetas, seed_phis,
            rmax=fl_settings["rmax"],
        )

st.plotly_chart(
    build_figure(
        planet, p2D, th2D, br_surf,
        cmap, vmin, vmax,
        field_lines=field_lines,
        line_color=fl_settings.get("line_color", "white") if fl_settings else "white",
        line_width=fl_settings.get("line_width", 2) if fl_settings else 2,
    ),
    width='stretch',
)

st.caption(
    "Drag to rotate · Scroll to zoom · "
    "Source: [planetMagFields](https://github.com/AnkitBarik/planetMagFields)"
)
