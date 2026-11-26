#!/usr/bin/env python3
"""
Coordinate-convention and getB regression tests.

phi is east longitude in the body-fixed frame of each field model, and Gauss
coefficients use Schmidt semi-normalised Legendre functions without the
Condon-Shortley phase, for every planet.
"""

import numpy as np
import pytest
from pathlib import Path
from scipy.special import sph_harm_y

from planetmagfields import Planet, planetlist, get_models
from planetmagfields.libgauss import getB

DATA = Path(__file__).parent
EARTH_RADIUS_KM = 6371.2


def getB_reference(lmax, mmax, glm, hlm, idx, r, phi, theta):
    """Direct, slow evaluation of Br with scipy spherical harmonics, removing
    scipy's Condon-Shortley phase."""
    br = np.zeros(np.broadcast(phi, theta).shape)
    for l in range(1, lmax + 1):
        for m in range(0, (l if mmax > 0 else 0) + 1):
            ylm = sph_harm_y(l, m, theta, phi)
            norm = (-1) ** m * np.sqrt(4 * np.pi / (2 * l + 1)) * (np.sqrt(2) if m else 1)
            br += (l + 1) * r ** (-l - 2) * norm * (
                glm[idx[l, m]] * ylm.real + hlm[idx[l, m]] * ylm.imag)
    return br


class TestEarthLongitude:
    def test_external_reference_values(self):
        """tests/earth/ref_val.dat holds IGRF vertical field Z (positive down,
        so Z = -Br) at geographic co-latitude/longitude, from an external
        calculator. phi must be the geographic east longitude."""
        for year, colat, lon, r_km, z_ref in np.loadtxt(DATA / 'earth/ref_val.dat'):
            p = Planet('earth', year=year, info=False, units='nT')
            r = np.array([r_km / EARTH_RADIUS_KM])
            theta = np.array([np.deg2rad(colat)])
            phi = np.array([np.deg2rad(lon)])
            p.orbit_path(r, theta, phi)
            # Tolerance covers small differences between IGRF generations.
            np.testing.assert_allclose(-p.br_orb[0], z_ref, rtol=1e-4, atol=1.0)
            br = getB(p.lmax, p.mmax, p.glm, p.hlm, p.idx, r[0], phi, theta)
            np.testing.assert_allclose(-br[0], z_ref, rtol=1e-4, atol=1.0)

    def test_south_atlantic_anomaly(self):
        """Total field near 25S, 50W is ~22,600 nT (IGRF 2020); the antipodal
        meridian (25S, 130E, northern Australia) is ~54,000 nT."""
        p = Planet('earth', year=2020, info=False, units='nT')
        p.orbit_path(np.array([1.0]), np.array([np.deg2rad(115.)]), np.array([np.deg2rad(-50.)]))
        F = np.sqrt(p.br_orb ** 2 + p.btheta_orb ** 2 + p.bphi_orb ** 2)[0]
        assert 21000 < F < 24000

    def test_dipole_pole_longitude(self):
        """The north geomagnetic pole of IGRF 2020 is at about 72.7 W."""
        p = Planet('earth', year=2020, info=False)
        g11, h11 = p.glm[p.idx[1, 1]], p.hlm[p.idx[1, 1]]
        lon = np.rad2deg(np.arctan2(-h11, -g11))
        assert abs(lon - (-72.7)) < 0.5

    def test_plot_longitude_matches_phi(self):
        """Map longitude is phi itself, wrapped to [-180, 180]."""
        from planetmagfields.plotlib import lonlat_grid
        p = Planet('earth', year=2020, nphi=64, info=False)
        lon2D, lat2D, B = lonlat_grid(p.p2D, p.th2D, p.Br)
        assert lon2D[0, 0] == -np.pi and lon2D[-1, 0] == np.pi
        assert np.all(np.diff(lon2D[:, 0]) > 0)
        # Interior columns are the original phi columns, unchanged
        for i in range(1, lon2D.shape[0] - 1):
            phi = lon2D[i, 0] % (2 * np.pi)
            j = np.argmin(np.abs(p.phi - phi))
            np.testing.assert_allclose(B[i], p.Br[j])


class TestJupiterLongitude:
    """Connerney et al. (2022, JGR Planets, doi:10.1029/2021JE007055), JRM33.
    System III (1965) is a west-longitude system: lambda_W = 360 - phi."""

    def test_dipole_tilt_and_longitude(self):
        p = Planet('jupiter', model='jrm33', info=False)
        g10, g11, h11 = (float(p.glm[p.idx[1, 0]]), float(p.glm[p.idx[1, 1]]),
                         float(p.hlm[p.idx[1, 1]]))
        tilt = np.rad2deg(np.arctan(np.hypot(g11, h11) / g10))
        lam_w = (-np.rad2deg(np.arctan2(h11, g11))) % 360
        assert abs(tilt - 10.25) < 0.01      # paper: 10.25 deg
        assert abs(lam_w - 196.38) < 0.01    # paper: towards 196.38 deg System III

    def test_great_blue_spot(self):
        """Isolated equatorial negative flux patch near 90 deg W (Fig. 7),
        about -12 G at the surface; on maps it must appear at that longitude."""
        from planetmagfields.plotlib import lonlat_grid
        p = Planet('jupiter', model='jrm33', nphi=360, info=False)
        lon2D, lat2D, B = lonlat_grid(p.p2D, p.th2D, p.Br)
        band = np.abs(np.rad2deg(lat2D)) < 30
        i, j = np.unravel_index(np.argmin(np.where(band, B, np.inf)), B.shape)
        assert abs(np.rad2deg(lon2D[i, j]) - (-87)) < 10   # ~84-90 W on the map
        assert -10 < np.rad2deg(lat2D[i, j]) < 0
        assert -16e2 < B[i, j] < -10e2                    # muT: about -12 G


class TestFastGetB:
    @pytest.mark.parametrize("planet,model", [(p, str(m)) for p in planetlist for m in get_models(p)])
    def test_matches_direct_evaluation(self, planet, model):
        p = Planet(planet, model=model, nphi=32, info=False, units='nT')
        rng = np.random.default_rng(1)
        theta = rng.uniform(0, np.pi, 50)
        phi = rng.uniform(0, 2 * np.pi, 50)
        for r in (0.8, 1.0, 2.5):
            ref = getB_reference(p.lmax, p.mmax, p.glm, p.hlm, p.idx, r, phi, theta)
            br = getB(p.lmax, p.mmax, p.glm, p.hlm, p.idx, r, phi, theta)
            np.testing.assert_allclose(br, ref, rtol=1e-10, atol=1e-10 * np.abs(ref).max())

    def test_grid_shape_and_values(self):
        p = Planet('jupiter', nphi=64, info=False, units='nT')
        assert p.Br.shape == p.p2D.shape
        ref = getB_reference(p.lmax, p.mmax, p.glm, p.hlm, p.idx, 1.0, p.p2D, p.th2D)
        np.testing.assert_allclose(p.Br, ref, rtol=1e-10, atol=1e-10 * np.abs(ref).max())

    def test_low_degree_closed_forms(self):
        """Single modes against closed-form Schmidt functions at r = 1."""
        from planetmagfields.libgauss import gen_idx
        lmax = 2
        idx = gen_idx(lmax)
        n = idx[lmax, lmax] + 1
        th = np.linspace(0.1, 3.0, 7)
        ph = np.linspace(0.2, 6.0, 7)
        c, s_ = np.cos(th), np.sin(th)
        cases = {(1, 0): c, (1, 1): s_, (2, 0): 1.5 * c ** 2 - 0.5,
                 (2, 1): np.sqrt(3) * c * s_, (2, 2): np.sqrt(3) / 2 * s_ ** 2}
        for (l, m), P in cases.items():
            g = np.zeros(n); h = np.zeros(n)
            g[idx[l, m]] = 1.0; h[idx[l, m]] = 0.5
            br = getB(lmax, lmax, g, h, idx, 1.0, ph, th)
            expect = (l + 1) * P * (np.cos(m * ph) + 0.5 * np.sin(m * ph))
            np.testing.assert_allclose(br, expect, atol=1e-13)

    def test_integer_radius_and_float32_coefficients(self):
        """r passed as an int and float32 coefficient files must not lose precision."""
        p = Planet('jupiter', model='jrm09', nphi=32, info=False, units='nT')
        assert p.glm.dtype == np.float32
        br_int = getB(p.lmax, p.mmax, p.glm, p.hlm, p.idx, 1, p.p2D, p.th2D)
        ref = getB_reference(p.lmax, p.mmax, p.glm.astype(float), p.hlm.astype(float),
                             p.idx, 1.0, p.p2D, p.th2D)
        np.testing.assert_allclose(br_int, ref, rtol=1e-10, atol=1e-10 * np.abs(ref).max())


class TestMapCentre:
    def test_default_central_lon(self):
        from planetmagfields.plotlib import default_central_lon
        assert default_central_lon('earth') == 0
        assert default_central_lon('Earth') == 0
        for name in ('mercury', 'jupiter', 'saturn', 'uranus', 'neptune', 'ganymede'):
            assert default_central_lon(name) == 180
        assert default_central_lon('jupiter', 0) == 0
        assert default_central_lon('earth', 90) == 90

    @pytest.mark.parametrize("clon", [0, 180, 90, -60])
    def test_lonlat_grid_centre(self, clon):
        """Longitudes span [clon-180, clon+180] and data stay with their phi."""
        from planetmagfields.plotlib import lonlat_grid
        p = Planet('jupiter', model='jrm33', nphi=64, info=False)
        lon2D, lat2D, B = lonlat_grid(p.p2D, p.th2D, p.Br, central_lon=clon)
        lon = np.rad2deg(lon2D[:, 0])
        np.testing.assert_allclose([lon[0], lon[-1]], [clon - 180, clon + 180])
        assert np.all(np.diff(lon) > 0)
        for i in range(1, len(lon) - 1):
            j = np.argmin(np.abs((p.phi - np.deg2rad(lon[i]) + np.pi) % (2*np.pi) - np.pi))
            np.testing.assert_allclose(B[i], p.Br[j])
