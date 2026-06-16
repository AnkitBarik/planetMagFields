#!/usr/bin/env python3
"""
Streamlit page for computing magnetic field along an orbit trajectory.
Users can upload coordinates or paste them directly, then download field values.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import io

from planetmagfields import Planet, get_models
from planetmagfields.libgauss import getB
from planetmagfields.utils import planetlist

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="Orbit Trajectory",
    page_icon="🛰️",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Cached helpers
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def fetch_models(planet_name: str) -> list[str]:
    return sorted(get_models(planet_name).tolist())


@st.cache_resource(show_spinner="Loading planet...")
def load_planet(name: str, model: str, units: str) -> Planet:
    return Planet(name=name, model=model, r=1.0, units=units, info=False)


@st.cache_data(show_spinner="Computing field along trajectory...")
def compute_orbit_field(
    planet_name: str,
    model: str,
    units: str,
    r: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute magnetic field components along orbit trajectory.
    """
    planet = load_planet(planet_name, model, units)
    planet.orbit_path(r, theta, phi)
    return planet.br_orb.copy(), planet.btheta_orb.copy(), planet.bphi_orb.copy()


@st.cache_data(show_spinner=False)
def compute_surface_br(
    planet_name: str, model: str, units: str, nphi: int = 128
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute Br on planet surface for visualization."""
    planet = load_planet(planet_name, model, units)
    ntheta = nphi // 2
    phi = np.linspace(0, 2 * np.pi, nphi + 1)
    theta = np.linspace(0, np.pi, ntheta + 1)
    p2D, th2D = np.meshgrid(phi, theta, indexing='ij')
    br = getB(
        planet.lmax, planet.mmax,
        planet.glm, planet.hlm, planet.idx,
        1.0, p2D, th2D,
        planetname=planet.name,
    )
    br *= planet.unitfac
    return p2D, th2D, br


def parse_csv_data(data: str | io.BytesIO, source: str = "input") -> pd.DataFrame | None:
    """Parse CSV data from string or file and validate columns."""
    try:
        if isinstance(data, str):
            if not data.strip():
                return None
            df = pd.read_csv(io.StringIO(data))
        else:
            df = pd.read_csv(data)

        df.columns = df.columns.str.lower().str.strip()

        required = {'r', 'theta', 'phi'}
        if not required.issubset(set(df.columns)):
            st.error(f"Data must contain columns: r, theta, phi. Found: {list(df.columns)}")
            return None

        for col in ['r', 'theta', 'phi']:
            df[col] = pd.to_numeric(df[col], errors='coerce')

        n_before = len(df)
        df = df.dropna(subset=['r', 'theta', 'phi'])
        n_after = len(df)

        if n_after == 0:
            st.error("No valid numeric data found.")
            return None

        if n_after < n_before:
            st.warning(f"Dropped {n_before - n_after} rows with invalid/missing values.")

        return df

    except Exception as e:
        st.error(f"Error parsing {source}: {e}")
        return None


def convert_angles(df: pd.DataFrame, theta_unit: str, phi_unit: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert coordinates to required units (r in Rp, angles in radians)."""
    r = df['r'].values.copy()
    theta = df['theta'].values.copy()
    phi = df['phi'].values.copy()

    if theta_unit == "degrees":
        theta = np.deg2rad(theta)

    if phi_unit == "degrees":
        phi = np.deg2rad(phi)

    return r, theta, phi


def build_trajectory_figure(
    r: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    br: np.ndarray,
    btheta: np.ndarray,
    bphi: np.ndarray,
    units: str,
    color_by: str,
    planet_name: str,
    model: str,
) -> go.Figure:
    """Build 3D figure showing trajectory colored by field component with Br surface."""

    # Convert trajectory to Cartesian
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    # Select color data for trajectory
    if color_by == "Br":
        color_data = br
        cbar_title = "Bᵣ"
    elif color_by == "Bθ":
        color_data = btheta
        cbar_title = "Bθ"
    elif color_by == "Bφ":
        color_data = bphi
        cbar_title = "Bφ"
    else:  # |B|
        color_data = np.sqrt(br**2 + btheta**2 + bphi**2)
        cbar_title = "|B|"

    # Unit label
    if units == "nT":
        unit_label = "nT"
    elif units == "muT":
        unit_label = "μT"
    else:
        unit_label = "G"

    cbar_title_traj = f"{cbar_title} ({unit_label})"

    # Create hover text for trajectory
    hover_text = [
        f"r: {ri:.3f} Rp<br>"
        f"θ: {np.rad2deg(ti):.2f}°<br>"
        f"φ: {np.rad2deg(pi):.2f}°<br>"
        f"Br: {bri:.4g} {unit_label}<br>"
        f"Bθ: {bti:.4g} {unit_label}<br>"
        f"Bφ: {bpi:.4g} {unit_label}"
        for ri, ti, pi, bri, bti, bpi in zip(r, theta, phi, br, btheta, bphi)
    ]

    # Get cached surface Br
    p2D, th2D, br_surf = compute_surface_br(planet_name, model, units)

    # Surface coordinates
    xs = np.sin(th2D) * np.cos(p2D)
    ys = np.sin(th2D) * np.sin(p2D)
    zs = np.cos(th2D)

    # Surface hover text
    abs_max_surf = float(np.abs(br_surf).max())
    use_scientific = (abs_max_surf >= 1000) | ((abs_max_surf < 0.01) & (abs_max_surf > 0))

    hover_text_surf = np.where(
        use_scientific,
        np.char.add(np.char.add("Br: ", np.char.mod("%.2e", br_surf)), f" {unit_label}"),
        np.char.add(np.char.add("Br: ", np.char.mod("%.3f", br_surf)), f" {unit_label}"),
    )

    # Planet surface with Br
    planet_surface = go.Surface(
        x=xs, y=ys, z=zs,
        surfacecolor=br_surf,
        colorscale="RdBu",
        reversescale=True,
        cmin=-abs_max_surf,
        cmax=abs_max_surf,
        colorbar=dict(
            title=dict(text=f"Surface Bᵣ ({unit_label})", side="right", font=dict(size=14)),
            thickness=12,
            tickfont=dict(size=10),
            x=1.0,
            len=0.4,
            y=0.25,
        ),
        lighting=dict(ambient=0.7, diffuse=0.6, specular=0.1),
        text=hover_text_surf,
        hovertemplate="%{text}<extra></extra>",
        name="Surface Br",
    )

    # Trajectory trace
    trajectory = go.Scatter3d(
        x=x, y=y, z=z,
        mode="lines+markers",
        marker=dict(
            size=4,
            color=color_data,
            colorscale="Viridis",
            colorbar=dict(
                title=dict(text=cbar_title_traj, side="right", font=dict(size=14)),
                thickness=12,
                tickfont=dict(size=10),
                x=1.0,
                len=0.4,
                y=0.75,
            ),
            cmin=color_data.min(),
            cmax=color_data.max(),
        ),
        line=dict(color="white", width=2),
        text=hover_text,
        hoverinfo="text",
        name="Trajectory",
    )

    fig = go.Figure(data=[planet_surface, trajectory])

    # Set axis ranges
    rmax = max(r.max() * 1.1, 1.5)
    axis_range = [-rmax, rmax]

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False, range=axis_range),
            yaxis=dict(visible=False, range=axis_range),
            zaxis=dict(visible=False, range=axis_range),
            bgcolor="black",
            aspectmode="cube",
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.0)),
        ),
        paper_bgcolor="#0e1117",
        margin=dict(l=0, r=0, t=0, b=0),
        height=550,
    )

    return fig


def build_timeseries_figure(
    br: np.ndarray,
    btheta: np.ndarray,
    bphi: np.ndarray,
    units: str,
) -> go.Figure:
    """Build figure showing field components vs point index."""

    if units == "nT":
        unit_label = "nT"
    elif units == "muT":
        unit_label = "μT"
    else:
        unit_label = "G"

    n_points = len(br)
    x = np.arange(n_points)

    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=(f"Bᵣ ({unit_label})", f"Bθ ({unit_label})", f"Bφ ({unit_label})"),
    )

    fig.add_trace(
        go.Scatter(x=x, y=br, mode="lines", name="Bᵣ", line=dict(color="#EF553B")),
        row=1, col=1,
    )
    fig.add_trace(
        go.Scatter(x=x, y=btheta, mode="lines", name="Bθ", line=dict(color="#00CC96")),
        row=2, col=1,
    )
    fig.add_trace(
        go.Scatter(x=x, y=bphi, mode="lines", name="Bφ", line=dict(color="#636EFA")),
        row=3, col=1,
    )

    fig.update_xaxes(title_text="Point index", row=3, col=1)
    fig.update_layout(
        height=450,
        showlegend=False,
        paper_bgcolor="#0e1117",
        plot_bgcolor="#1a1a2e",
        font=dict(color="white"),
        margin=dict(l=60, r=20, t=40, b=40),
    )

    for annotation in fig['layout']['annotations']:
        annotation['font'] = dict(size=14, color='white')

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
    st.title("Orbit Settings")

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

    units = st.selectbox("Field units", ["muT", "nT", "Gauss"])

    st.divider()
    st.subheader("Input format")

    theta_unit = st.radio("θ (colatitude) unit", ["radians", "degrees"], horizontal=True)
    phi_unit = st.radio("φ (longitude) unit", ["radians", "degrees"], horizontal=True)

    st.divider()
    st.subheader("Visualization")

    color_by = st.selectbox("Color trajectory by", ["Br", "Bθ", "Bφ", "|B|"])


# ---------------------------------------------------------------------------
# Main panel
# ---------------------------------------------------------------------------

st.title("🛰️ Orbit Trajectory Field Calculator")

st.markdown("""
Compute the magnetic field along a spacecraft/satellite trajectory.
Provide coordinates via **file upload** or **paste directly** below.

**Required columns:** `r`, `theta`, `phi`
- **r**: Radial distance in planetary radii (Rp)
- **theta**: Colatitude (0 = north pole, π = south pole)
- **phi**: Longitude (azimuthal angle)
""")

# Initialize session state for trajectory data
if "trajectory_df" not in st.session_state:
    st.session_state.trajectory_df = None

# Input method tabs
input_tab1, input_tab2, input_tab3 = st.tabs(["📋 Paste Data", "📁 Upload File", "📝 Generate Example"])

with input_tab1:
    st.markdown("Paste CSV data with columns `r`, `theta`, `phi`:")

    example_text = """r,theta,phi
1.5,1.57,0.0
1.6,1.50,0.3
1.8,1.40,0.6
2.0,1.30,0.9
2.3,1.20,1.2
2.6,1.10,1.5
3.0,1.00,1.8
3.4,0.90,2.1
3.8,0.85,2.4
4.0,0.80,2.7"""

    pasted_data = st.text_area(
        "CSV data",
        height=200,
        placeholder=example_text,
        help="Paste comma-separated values with header row",
        label_visibility="collapsed",
    )

    col1, col2, _, _, _ = st.columns(5)

    with col1:
        if st.button("Compute field along trajectory", type="primary", width='content'):
            if pasted_data.strip():
                parsed = parse_csv_data(pasted_data, source="pasted data")
                if parsed is not None:
                    st.session_state.trajectory_df = parsed
                    st.rerun()
            else:
                st.warning("Please paste some data first.")

    with col2:
        if st.button("Load example", width='content'):
            st.session_state.trajectory_df = parse_csv_data(example_text, source="example")
            st.rerun()

with input_tab2:
    uploaded_file = st.file_uploader(
        "Upload trajectory coordinates (CSV)",
        type=["csv"],
        help="CSV file with columns: r, theta, phi",
    )

    if uploaded_file is not None:
        parsed = parse_csv_data(uploaded_file, source="uploaded file")
        if parsed is not None:
            st.session_state.trajectory_df = parsed

with input_tab3:
    st.markdown("Generate a sample elliptical orbit for testing:")

    col1, col2, col3 = st.columns(3)
    with col1:
        n_points = st.number_input("Number of points", min_value=10, max_value=1000, value=100)
    with col2:
        r_min = st.number_input("Periapsis (Rp)", min_value=1.0, max_value=5.0, value=1.5, step=0.1)
    with col3:
        r_max = st.number_input("Apoapsis (Rp)", min_value=1.5, max_value=10.0, value=4.0, step=0.5)

    inclination = st.slider("Orbital inclination (°)", min_value=0, max_value=90, value=45)

    gen_col1, gen_col2, _, _, _ = st.columns(5)

    with gen_col1:
        if st.button("🎲 Generate & compute", type="primary", width='content'):
            t = np.linspace(0, 2 * np.pi, int(n_points))
            a = (r_min + r_max) / 2
            e = (r_max - r_min) / (r_max + r_min)
            r_orbit = a * (1 - e**2) / (1 + e * np.cos(t))

            inc_rad = np.deg2rad(inclination)
            theta_orbit = np.pi/2 - np.arcsin(np.sin(inc_rad) * np.sin(t))
            phi_orbit = np.arctan2(np.cos(inc_rad) * np.sin(t), np.cos(t))
            phi_orbit = phi_orbit % (2 * np.pi)

            st.session_state.trajectory_df = pd.DataFrame({
                'r': r_orbit,
                'theta': np.rad2deg(theta_orbit) if theta_unit == "degrees" else theta_orbit,
                'phi': np.rad2deg(phi_orbit) if phi_unit == "degrees" else phi_orbit,
            })
            st.rerun()

    with gen_col2:
        t = np.linspace(0, 2 * np.pi, int(n_points))
        a = (r_min + r_max) / 2
        e = (r_max - r_min) / (r_max + r_min)
        r_orbit = a * (1 - e**2) / (1 + e * np.cos(t))
        inc_rad = np.deg2rad(inclination)
        theta_orbit = np.pi/2 - np.arcsin(np.sin(inc_rad) * np.sin(t))
        phi_orbit = np.arctan2(np.cos(inc_rad) * np.sin(t), np.cos(t)) % (2 * np.pi)

        example_df = pd.DataFrame({
            'r': r_orbit,
            'theta': np.rad2deg(theta_orbit) if theta_unit == "degrees" else theta_orbit,
            'phi': np.rad2deg(phi_orbit) if phi_unit == "degrees" else phi_orbit,
        })

        st.download_button(
            "⬇️ Download example CSV",
            data=example_df.to_csv(index=False),
            file_name="example_orbit.csv",
            mime="text/csv",
            width='content',
        )

# Clear data button
if st.session_state.trajectory_df is not None:
    if st.button("🗑️ Clear trajectory data"):
        st.session_state.trajectory_df = None
        st.rerun()

# ---------------------------------------------------------------------------
# Process and display results
# ---------------------------------------------------------------------------

df = st.session_state.trajectory_df

if df is not None:
    st.divider()
    st.success(f"✅ Loaded {len(df)} trajectory points")

    with st.expander("Preview input data", expanded=False):
        st.dataframe(df.head(20), width='content')

    # Convert coordinates
    r, theta, phi = convert_angles(df, theta_unit, phi_unit)

    # Validate ranges
    if (r < 0.5).any():
        st.warning("⚠️ Some r values are < 0.5 Rp. Results may be unreliable inside the planet.")
    if (theta < 0).any() or (theta > np.pi).any():
        st.warning("⚠️ θ values should be in [0, π] radians (colatitude).")

    # Compute field (cached)
    br, btheta, bphi = compute_orbit_field(planet_name, model, units, r, theta, phi)

    # Create output dataframe
    output_df = df.copy()
    output_df['Br'] = br
    output_df['Btheta'] = btheta
    output_df['Bphi'] = bphi
    output_df['Bmag'] = np.sqrt(br**2 + btheta**2 + bphi**2)

    # Unit label
    if units == "muT":
        unit_label = "μT"
    elif units == "nT":
        unit_label = "nT"
    else:
        unit_label = "G"

    # Results header
    res_col1, res_col2 = st.columns([1, 3])

    with res_col1:
        st.download_button(
            "⬇️ Download results CSV",
            data=output_df.to_csv(index=False),
            file_name=f"{planet_name}_{model}_orbit_field.csv",
            mime="text/csv",
            type="primary",
            width='content',
        )

    with res_col2:
        bmag = output_df['Bmag']
        st.markdown(
            f"**Field magnitude range:** {bmag.min():.4g} – {bmag.max():.4g} {unit_label} &nbsp;|&nbsp; "
            f"**Planet:** {planet_name.capitalize()} &nbsp;|&nbsp; **Model:** {model}"
        )

    # Visualizations
    viz_tab1, viz_tab2, viz_tab3 = st.tabs(["🌐 3D Trajectory", "📈 Field Components", "📊 Data Table"])

    with viz_tab1:
        fig_3d = build_trajectory_figure(
            r, theta, phi, br, btheta, bphi, units, color_by,
            planet_name, model,
        )
        st.plotly_chart(fig_3d, width="stretch")
        st.caption("Drag to rotate · Scroll to zoom · Surface shows Bᵣ at r = 1 Rp")

    with viz_tab2:
        fig_ts = build_timeseries_figure(br, btheta, bphi, units)
        st.plotly_chart(fig_ts, width="stretch")

    with viz_tab3:
        st.dataframe(output_df, width="stretch", height=400)

else:
    st.info("👆 Paste data, upload a file, or generate an example trajectory to get started.")