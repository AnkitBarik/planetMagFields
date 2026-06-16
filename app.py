#!/usr/bin/env python3
"""
Streamlit webapp for interactive 3D visualization of planetary magnetic fields.
Uses planetmagfields' built-in potential extrapolation for fast field line tracing.
"""

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy.integrate import solve_ivp
from scipy.interpolate import RegularGridInterpolator
import warnings

from planetmagfields import Planet, get_models
from planetmagfields.libgauss import getB
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


@st.cache_resource(show_spinner="Loading planet...")
def load_planet(name: str, model: str, r: float, units: str, nphi: int) -> Planet:
    return Planet(name=name, model=model, r=r, units=units, nphi=nphi, info=False)


@st.cache_data(show_spinner=False)
def compute_surface_br(
    name: str, model: str, r: float, units: str, nphi: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute Br on a uniform (regular) grid to avoid pole holes in the 3D sphere."""
    planet = load_planet(name, model, r, units, nphi)
    ntheta = nphi // 2
    phi = np.linspace(0, 2 * np.pi, nphi + 1)
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
# Field extrapolation and interpolation using Planet.extrapolate()
# ---------------------------------------------------------------------------

@st.cache_resource(show_spinner="Extrapolating field to 3D grid...")
def get_extrapolated_field(
    name: str, model: str, units: str, nphi: int, rmax: float, nr: int = 40
) -> tuple[Planet, np.ndarray]:
    """
    Use Planet.extrapolate() to compute B-field on radial shells.
    Returns the planet object (with extrapolated field) and r array.

    Field arrays have shape (nphi, ntheta, nr).
    """
    planet = Planet(name=name, model=model, r=1.0, units=units, nphi=nphi, info=False)
    rout = np.linspace(1.0, rmax, nr)
    planet.extrapolate(rout)
    return planet, rout


def create_field_interpolators(
    planet: Planet, rout: np.ndarray
) -> tuple[RegularGridInterpolator, RegularGridInterpolator, RegularGridInterpolator]:
    """
    Create interpolators from the extrapolated field arrays.

    Field arrays have shape (nphi, ntheta, nr).
    Grid order for interpolator: (phi, theta, r).
    Query points: [phi, theta, r].
    """
    grid = (planet.phi, planet.theta, rout)

    interp_br = RegularGridInterpolator(
        grid, planet.br_ex, method='linear', bounds_error=False, fill_value=None
    )
    interp_bt = RegularGridInterpolator(
        grid, planet.btheta_ex, method='linear', bounds_error=False, fill_value=None
    )
    interp_bp = RegularGridInterpolator(
        grid, planet.bphi_ex, method='linear', bounds_error=False, fill_value=None
    )

    return interp_br, interp_bt, interp_bp


# ---------------------------------------------------------------------------
# Automatic seeding (similar to PyVista streamlines)
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def get_seed_points(
    _planet_name: str,
    _planet_model: str,
    _planet_units: str,
    _nphi: int,
    n_lines: int,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Automatically generate seed points for field line tracing.

    Strategy (similar to PyVista's streamlines):
    1. Use |Br| as a probability distribution to bias seeds toward strong field regions
    2. Ensure coverage of both positive and negative Br regions (for closed field lines)
    3. Add some uniform random seeds for overall coverage
    """
    planet = load_planet(_planet_name, _planet_model, 1.0, _planet_units, _nphi)

    rng = np.random.default_rng(rng_seed)
    br = planet.Br  # shape: (nphi, ntheta)
    br_abs = np.abs(br)

    n_positive = max(1, int(0.4 * n_lines))
    n_negative = max(1, int(0.4 * n_lines))
    n_uniform = n_lines - n_positive - n_negative

    thetas = []
    phis = []

    # Seeds in strong POSITIVE Br regions
    positive_mask = br > 0
    if positive_mask.any():
        weights_pos = np.where(positive_mask, br_abs, 0).ravel()
        weights_pos = weights_pos / weights_pos.sum()
        indices = rng.choice(len(weights_pos), size=n_positive, replace=False, p=weights_pos)
        iphi, itheta = np.unravel_index(indices, br.shape)
        thetas.extend(planet.theta[itheta])
        phis.extend(planet.phi[iphi])

    # Seeds in strong NEGATIVE Br regions
    negative_mask = br < 0
    if negative_mask.any():
        weights_neg = np.where(negative_mask, br_abs, 0).ravel()
        weights_neg = weights_neg / weights_neg.sum()
        indices = rng.choice(len(weights_neg), size=n_negative, replace=False, p=weights_neg)
        iphi, itheta = np.unravel_index(indices, br.shape)
        thetas.extend(planet.theta[itheta])
        phis.extend(planet.phi[iphi])

    # Uniform random seeds
    if n_uniform > 0:
        rand_theta = np.arccos(rng.uniform(-1.0, 1.0, n_uniform))
        rand_phi = rng.uniform(0.0, 2 * np.pi, n_uniform)
        thetas.extend(rand_theta)
        phis.extend(rand_phi)

    return np.array(thetas), np.array(phis)


# ---------------------------------------------------------------------------
# Smooth adaptive field-line tracing with dense_output (CACHED)
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner="Tracing field lines...")
def compute_field_lines(
    _planet_name: str,
    _planet_model: str,
    _planet_units: str,
    _nphi: int,
    rmax: float,
    nr: int,
    seed_thetas: np.ndarray,
    seed_phis: np.ndarray,
    max_arc_length: float = 100.0,
    points_per_unit_length: int = 50,
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Trace field lines using adaptive ODE solver with dense_output for smooth curves.

    Parameters
    ----------
    _planet_name, _planet_model, _planet_units, _nphi : str/int
        Planet parameters for cache key (prefixed with _ to exclude from hashing).
    rmax : float
        Maximum radius for tracing.
    nr : int
        Number of radial grid points.
    seed_thetas, seed_phis : ndarray
        Seed point coordinates.
    max_arc_length : float
        Maximum arc length to integrate.
    points_per_unit_length : int
        Number of points per unit arc length for smooth rendering.

    State vector: [phi, theta, r] to match (nphi, ntheta, nr) convention.
    """
    # Load extrapolated field inside cached function
    planet, rout = get_extrapolated_field(
        _planet_name, _planet_model, _planet_units, _nphi, rmax, nr
    )

    interp_br, interp_bt, interp_bp = create_field_interpolators(planet, rout)

    theta_min, theta_max = planet.theta.min(), planet.theta.max()
    r_min, r_max_grid = rout.min(), rout.max()

    def ode_rhs(t, state, direction):
        """ODE for field line tracing. State: [phi, theta, r]."""
        phi, theta, r = state

        phi_wrapped = phi % (2 * np.pi)
        theta_clamped = np.clip(theta, theta_min + 1e-6, theta_max - 1e-6)
        r_clamped = np.clip(r, r_min, r_max_grid)

        point = np.array([[phi_wrapped, theta_clamped, r_clamped]])

        br = float(interp_br(point)[0])
        bt = float(interp_bt(point)[0])
        bp = float(interp_bp(point)[0])

        mag = np.sqrt(br**2 + bt**2 + bp**2) + 1e-12
        sin_theta = np.sin(theta_clamped) + 1e-10

        # Derivatives: [dphi/ds, dtheta/ds, dr/ds]
        return [
            direction * bp / (r_clamped * sin_theta * mag),
            direction * bt / (r_clamped * mag),
            direction * br / mag,
        ]

    def hit_surface(t, state, direction):
        return state[2] - 1.0  # r is state[2]

    def hit_rmax(t, state, direction):
        return state[2] - rmax  # r is state[2]

    hit_surface.terminal = True
    hit_surface.direction = -1
    hit_rmax.terminal = True
    hit_rmax.direction = 1

    lines = []
    r_start = 1.02

    for theta0, phi0 in zip(seed_thetas, seed_phis):
        points_forward = []
        points_backward = []

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            # Forward integration: state = [phi, theta, r]
            fwd = solve_ivp(
                ode_rhs,
                [0, max_arc_length],
                [phi0, theta0, r_start],
                args=(1,),
                method='RK45',
                events=[hit_surface, hit_rmax],
                dense_output=True,
                rtol=1e-6,
                atol=1e-8,
            )

            # Backward integration
            bwd = solve_ivp(
                ode_rhs,
                [0, max_arc_length],
                [phi0, theta0, r_start],
                args=(-1,),
                method='RK45',
                events=[hit_surface, hit_rmax],
                dense_output=True,
                rtol=1e-6,
                atol=1e-8,
            )

        # Resample forward path at uniform arc-length intervals
        if fwd.success and fwd.t[-1] > 1e-6:
            n_pts_fwd = max(20, int(points_per_unit_length * fwd.t[-1]))
            s_fwd = np.linspace(0, fwd.t[-1], n_pts_fwd)
            states_fwd = fwd.sol(s_fwd)  # shape: (3, n_pts)
            for i in range(len(s_fwd)):
                phi, theta, r = states_fwd[0, i], states_fwd[1, i], states_fwd[2, i]
                if 1.0 - 0.01 <= r <= rmax + 0.01:
                    points_forward.append([phi, theta, r])

        # Resample backward path at uniform arc-length intervals
        if bwd.success and bwd.t[-1] > 1e-6:
            n_pts_bwd = max(20, int(points_per_unit_length * bwd.t[-1]))
            s_bwd = np.linspace(0, bwd.t[-1], n_pts_bwd)
            states_bwd = bwd.sol(s_bwd)
            for i in range(len(s_bwd)):
                phi, theta, r = states_bwd[0, i], states_bwd[1, i], states_bwd[2, i]
                if 1.0 - 0.01 <= r <= rmax + 0.01:
                    points_backward.append([phi, theta, r])

        # Combine: backward (reversed) + forward (skip duplicate at seed)
        points_backward.reverse()
        if points_forward and points_backward:
            all_points = points_backward + points_forward[1:]
        elif points_forward:
            all_points = points_forward
        elif points_backward:
            all_points = points_backward
        else:
            continue

        if len(all_points) < 5:
            continue

        # Convert to Cartesian: points are [phi, theta, r]
        all_points = np.array(all_points)
        phi = all_points[:, 0]
        theta = all_points[:, 1]
        r = all_points[:, 2]

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        lines.append((x, y, z))

    return lines


# ---------------------------------------------------------------------------
# 3D figure builder - separated into cached traces and layout
# ---------------------------------------------------------------------------
def get_diverging_colormaps():
    """Get list of diverging colormap names from plotly, excluding reversed variants."""
    return sorted([
        name for name in dir(px.colors.diverging)
        if not name.startswith('_')
        and name[0].isupper()
        and not name.endswith('_r')
    ])

# Colormap selection
colormaps = get_diverging_colormaps()
default_cmap = "RdBu" if "RdBu" in colormaps else colormaps[0]


@st.cache_data(show_spinner=False)
def compute_surface_trace(
    _planet_name: str,
    _planet_model: str,
    _planet_units: str,
    p2D: np.ndarray,
    th2D: np.ndarray,
    br: np.ndarray,
    selected_cmap: str,
    use_reversescale: bool,
    vmin: float | None,
    vmax: float | None,
) -> dict:
    """Compute surface trace data (cached). Returns dict for go.Surface."""
    x = np.sin(th2D) * np.cos(p2D)
    y = np.sin(th2D) * np.sin(p2D)
    z = np.cos(th2D)

    abs_max = float(np.abs(br).max())
    cmin = vmin if vmin is not None else -abs_max
    cmax = vmax if vmax is not None else abs_max

    # Vectorized hover text creation
    use_scientific = (abs_max >= 1000) | ((abs_max < 0.01) & (abs_max > 0))

    hover_text = np.where(
        use_scientific,
        np.char.add(np.char.add("Br: ", np.char.mod("%.2e", br)), f" {_planet_units}"),
        np.char.add(np.char.add("Br: ", np.char.mod("%.3f", br)), f" {_planet_units}"),
    )

    if _planet_units == "nT":
        cbar_title = "Bᵣ (nT)"
    elif _planet_units == "muT":
        cbar_title = "Bᵣ (μT)"
    elif _planet_units == "Gauss":
        cbar_title = "Bᵣ (G)"
    else:
        cbar_title = f"Bᵣ ({_planet_units})"

    return {
        "x": x,
        "y": y,
        "z": z,
        "surfacecolor": br,
        "colorscale": selected_cmap,
        "reversescale": use_reversescale,
        "cmin": cmin,
        "cmax": cmax,
        "cbar_title": cbar_title,
        "hover_text": hover_text,
    }


@st.cache_data(show_spinner=False)
def compute_field_line_traces(
    field_lines: list[tuple[np.ndarray, np.ndarray, np.ndarray]] | None,
    line_color: str,
    line_width: int,
) -> list[dict]:
    """Compute field line trace data (cached). Returns list of dicts for go.Scatter3d."""
    if not field_lines:
        return []

    traces = []
    n_lines = len(field_lines)

    if line_color == "rainbow":
        colors = px.colors.sample_colorscale("plasma", n_lines)
    else:
        colors = [line_color] * n_lines

    for i, (lx, ly, lz) in enumerate(field_lines):
        traces.append({
            "x": lx,
            "y": ly,
            "z": lz,
            "color": colors[i],
            "width": line_width,
        })

    return traces


def build_figure(
    surface_data: dict,
    field_line_data: list[dict],
    show_axes: bool,
    rmax: float,
) -> go.Figure:
    """Build the 3D Plotly figure from cached trace data. Only layout depends on show_axes."""

    surface = go.Surface(
        x=surface_data["x"],
        y=surface_data["y"],
        z=surface_data["z"],
        surfacecolor=surface_data["surfacecolor"],
        colorscale=surface_data["colorscale"],
        reversescale=surface_data["reversescale"],
        cmin=surface_data["cmin"],
        cmax=surface_data["cmax"],
        colorbar=dict(
            title=dict(text=surface_data["cbar_title"], side="right", font=dict(size=20)),
            thickness=18,
            tickfont=dict(size=16),
        ),
        lighting=dict(ambient=0.7, diffuse=0.6, specular=0.1),
        text=surface_data["hover_text"],
        hovertemplate="%{text}<extra></extra>",
    )

    fig = go.Figure(data=[surface])

    for fl in field_line_data:
        fig.add_trace(go.Scatter3d(
            x=fl["x"],
            y=fl["y"],
            z=fl["z"],
            mode="lines",
            line=dict(color=fl["color"], width=fl["width"]),
            hoverinfo="none",
            showlegend=False,
        ))

    # Set symmetric axis ranges to keep rotation centered on origin
    axis_range = [-rmax, rmax]

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=show_axes, range=axis_range),
            yaxis=dict(visible=show_axes, range=axis_range),
            zaxis=dict(visible=show_axes, range=axis_range),
            bgcolor="black",
            aspectmode="cube",
            camera=dict(eye=dict(x=1.1, y=1.1, z=0.8)),
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
    "earth": "igrf14",
    "mercury": "wardinski2019",
    "jupiter": "jrm33",
    "saturn": "cassini11+",
    "uranus": "holme1996",
    "neptune": "connerny1991",
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
        st.write("")
        if st.button("↺", key="reset_r", help="Reset radius to 1.0"):
            st.session_state["r_val"] = 1.0
    with r_col:
        r = st.slider(
            "Radius (r / R_planet)",
            min_value=0.5, max_value=5.0,
            value=1.0, step=0.05,
            key="r_val",
        )

    units = st.selectbox("Field units", ["muT", "nT", "Gauss"])

    res_reset_col, res_col = st.columns([1, 4])
    with res_reset_col:
        st.write("")
        if st.button("↺", key="reset_res", help="Reset resolution to 128"):
            st.session_state["res_val"] = 128
    with res_col:
        resolution = st.select_slider(
            "Resolution (nphi)",
            options=[64, 128, 256, 512],
            value=128,
            key="res_val",
        )

    show_axes = st.checkbox("Show axes", value=False)

    st.divider()
    st.subheader("Colour scale")

    selected_cmap = st.selectbox("Colormap", colormaps, index=colormaps.index(default_cmap))

    # Invert checkbox - default unchecked means blue=negative, red=positive
    invert_cmap = st.checkbox("Invert colormap", value=False,
                            help="Flip color mapping (default: blue=negative, red=positive)")

    use_reversescale = not invert_cmap  # Default True, becomes False when inverted

    sym = st.checkbox("Symmetric (auto ±max)", value=True)
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
        fl_settings["n_lines"] = st.slider(
            "Number of field lines",
            min_value=1, max_value=50, value=20,
        )
        fl_settings["rmax"] = st.slider(
            "Max radius (Rp)",
            min_value=1.5, max_value=10.0, value=4.0, step=0.5
        )
        fl_settings["line_color"] = st.selectbox(
            "Line colour",
            ["white", "rainbow", "yellow", "cyan", "lime", "black"],
        )
        fl_settings["line_width"] = st.slider(
            "Line width",
            min_value=1, max_value=6, value=3
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

c1, c2, c3, c4, c5 = st.columns(5)

def _metric_html(label: str, value: str) -> str:
    return (
        f'<div style="text-align:center">'
        f'<div style="font-size:2.0rem;color:#aaa;margin-bottom:2px">{label}</div>'
        f'<div style="font-size:1.6rem;font-weight:600">{value}</div>'
        f'</div>'
    )

c1.markdown(_metric_html("Planet", planet.name.capitalize()), unsafe_allow_html=True)
c2.markdown(_metric_html("Model", planet.model), unsafe_allow_html=True)
c3.markdown(_metric_html("l_max", str(planet.lmax)), unsafe_allow_html=True)
c4.markdown(_metric_html("Dipole tilt", f"{planet.dipTheta:.2f}°"), unsafe_allow_html=True)
c5.markdown(_metric_html("r / Rp", f"{r:.2f}"), unsafe_allow_html=True)

# Compute field lines if enabled (cached)
field_lines = None
if show_lines and fl_settings:
    seed_thetas, seed_phis = get_seed_points(
        planet_name,
        model,
        units,
        resolution,
        n_lines=fl_settings["n_lines"],
    )

    field_lines = compute_field_lines(
        planet_name,
        model,
        units,
        resolution,
        rmax=fl_settings["rmax"],
        nr=50,
        seed_thetas=seed_thetas,
        seed_phis=seed_phis,
        max_arc_length=150.0,
        points_per_unit_length=50,
    )

# Compute cached trace data (does NOT depend on show_axes)
surface_data = compute_surface_trace(
    planet_name,
    model,
    units,
    p2D,
    th2D,
    br_surf,
    selected_cmap,
    use_reversescale,
    vmin,
    vmax,
)

field_line_data = compute_field_line_traces(
    field_lines,
    fl_settings.get("line_color", "white") if fl_settings else "white",
    fl_settings.get("line_width", 3) if fl_settings else 3,
)

# Build figure with layout (show_axes only affects layout, not cached data)
rmax_for_layout = fl_settings.get("rmax", 1.0) if fl_settings else 1.0

fig = build_figure(
    surface_data,
    field_line_data,
    show_axes,
    rmax_for_layout,
)

st.plotly_chart(fig, width="stretch")

st.caption(
    "Drag to rotate · Scroll to zoom · "
    "Documentation: [planetMagFields](https://github.com/AnkitBarik/planetMagFields) · "
    "Source: [planetMagFields](https://github.com/AnkitBarik/planetMagFields)"
)