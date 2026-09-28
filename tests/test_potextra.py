import numpy as np
import pytest
from pathlib import Path
from planetmagfields import Planet
from planetmagfields.libgauss import gen_idx, get_grid
from planetmagfields.potextra import (
    get_pol_from_Gauss,
    extrapot_scipy,
    get_field_along_path_scipy
)

DATA = Path(__file__).parent


def _make_single_mode(lmax, l1, m1, g_val=1.0, h_val=0.0):
    ncoeff = (lmax + 1) * (lmax + 2) // 2
    idx    = gen_idx(lmax)
    glm    = np.zeros(ncoeff)
    hlm    = np.zeros(ncoeff)
    glm[idx[l1, m1]] = g_val
    hlm[idx[l1, m1]] = h_val
    return glm, hlm, idx


class TestJupiterBr:
    def test_surface_br(self):
        p = Planet(name='jupiter', r=0.85, nphi=256, info=False, model='jrm33')

        br_ref = np.loadtxt(DATA / 'jupiter/Br_reference085.dat')
        percent_err = np.abs((br_ref - p.Br) / br_ref) * 100

        np.testing.assert_allclose(percent_err, 0, rtol=0.1, atol=0.1)


class TestPotentialExtrapolation:
    def test_internal_consistency_jupiter(self):
        p = Planet(name='jupiter', r=10, nphi=256, info=False, model='jrm33')
        p.extrapolate(np.array([10]))
        err = np.abs(np.squeeze(p.br_ex) - p.Br)

        np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)

    def test_jupiter_extrapolated_field(self):
        p = Planet(name='jupiter', nphi=512, info=False, model='jrm33')
        p.extrapolate([2])

        br_ref = np.loadtxt(DATA / 'jupiter/Br_reference.dat')
        bt_ref = np.loadtxt(DATA / 'jupiter/Bt_reference.dat')
        bp_ref = np.loadtxt(DATA / 'jupiter/Bp_reference.dat')

        percent_err = (
              np.mean(np.abs(p.br_ex[..., 0] * 1e3 - br_ref) / br_ref)
            + np.mean(np.abs(p.btheta_ex[..., 0] * 1e3 - bt_ref) / bt_ref)
            + np.mean(np.abs(p.bphi_ex[..., 0] * 1e3 - bp_ref) / bp_ref)
        ) * 100

        np.testing.assert_allclose(percent_err, 0, rtol=1, atol=1)

    def test_internal_consistency_saturn_m0(self):
        p = Planet(name='saturn', r=10, nphi=256, info=False, model='cassini11+')
        p.extrapolate(np.array([10]))
        err = np.abs(np.squeeze(p.br_ex) - p.Br)

        np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)


class TestOrbitPath:
    def test_orbit_earth(self):
        p = Planet(name='earth', r=1, nphi=256, info=False, model='igrf14')
        p.extrapolate([2])
        p.orbit_path([2], [p.theta[10]], [p.phi[10]])

        br_ref  = p.br_ex[10, 10, 0]
        bt_ref  = p.btheta_ex[10, 10, 0]
        bp_ref  = p.bphi_ex[10, 10, 0]
        err = np.sqrt(
            (p.br_orb - br_ref) ** 2
            + (p.btheta_orb - bt_ref) ** 2
            + (p.bphi_orb - bp_ref) ** 2
        )

        np.testing.assert_allclose(err, 0, rtol=1e-3, atol=1e-3)

    def test_orbit_saturn_m0(self):
        p = Planet(name='saturn', r=1, nphi=256, info=False, model='cassini11+')
        p.extrapolate([2])
        p.orbit_path([2], [p.theta[10]], [p.phi[10]])

        br_ref  = p.br_ex[10, 10, 0]
        bt_ref  = p.btheta_ex[10, 10, 0]
        bp_ref  = p.bphi_ex[10, 10, 0]
        err = np.sqrt(
            (p.br_orb - br_ref) ** 2
            + (p.btheta_orb - bt_ref) ** 2
            + (p.bphi_orb - bp_ref) ** 2
        )

        np.testing.assert_allclose(err, 0, rtol=1e-4, atol=1e-4)


class TestScipyImplementations:
    def test_extrapot_scipy_vs_field_along_path(self):
        """extrapot_scipy (full 3-D grid) and get_field_along_path_scipy
        (arbitrary points) must agree to floating-point precision."""
        lmax = 3
        glm, hlm, idx = _make_single_mode(lmax, 2, 1, g_val=1.0, h_val=0.5)
        r_test = 1.5
        nphi   = 32
        ntheta = nphi // 2

        br_grid, bt_grid, bp_grid = extrapot_scipy(
            glm, hlm, idx, lmax, lmax, 1.0, np.array([r_test]), nphi=nphi
        )

        _, _, phi, theta = get_grid(nphi, ntheta)

        for i_phi, i_theta in [(5, 4), (10, 7), (20, 12)]:
            br_pt, bt_pt, bp_pt = get_field_along_path_scipy(
                glm, hlm, idx, lmax,
                np.array([r_test]),
                np.array([theta[i_theta]]),
                np.array([phi[i_phi]]),
            )
            np.testing.assert_allclose(
                br_grid[i_phi, i_theta, 0], br_pt[0], rtol=1e-10, atol=1e-10
            )
            np.testing.assert_allclose(
                bt_grid[i_phi, i_theta, 0], bt_pt[0], rtol=1e-10, atol=1e-10
            )
            np.testing.assert_allclose(
                bp_grid[i_phi, i_theta, 0], bp_pt[0], rtol=1e-10, atol=1e-10
            )


class TestGetPolFromGauss:
    def test_m0_value(self):
        """For g(1,0)=1 the poloidal coefficient at (l=1,m=0) must equal
        sqrt(4*pi/3)/l; all other modes must be zero."""
        lmax = 2
        glm, hlm, idx = _make_single_mode(lmax, 1, 0)
        bpol = get_pol_from_Gauss('earth', glm, hlm, lmax, lmax, idx)

        expected = np.sqrt(4 * np.pi / 3) / 1
        np.testing.assert_allclose(bpol[idx[1, 0]], expected, rtol=1e-14)
        bpol[idx[1, 0]] = 0.0
        np.testing.assert_allclose(bpol, 0.0, atol=1e-14)

    def test_fac_m_earth_vs_nonearth(self):
        """For m>0, earth uses fac_m=1 while other planets use (-1)^m."""
        lmax = 2
        glm, hlm, idx = _make_single_mode(lmax, 1, 1)
        norm_m1 = np.sqrt(2 * np.pi / 3) / 1

        bpol_earth = get_pol_from_Gauss('earth',   glm, hlm, lmax, lmax, idx)
        bpol_jup   = get_pol_from_Gauss('jupiter', glm, hlm, lmax, lmax, idx)

        np.testing.assert_allclose(bpol_earth[idx[1, 1]],  norm_m1, rtol=1e-14)
        np.testing.assert_allclose(bpol_jup[idx[1, 1]],   -norm_m1, rtol=1e-14)

    def test_mmax0_ignores_nonzero_m(self):
        """When mmax=0, modes with m>0 must not contribute."""
        lmax = 3
        glm, hlm, idx = _make_single_mode(lmax, 2, 1, g_val=5.0, h_val=3.0)
        bpol = get_pol_from_Gauss('earth', glm, hlm, lmax, 0, idx)

        np.testing.assert_allclose(bpol, 0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# extrapot_scipy – radial power-law scaling
# ---------------------------------------------------------------------------

def test_extrapot_scipy_radial_scaling():
    """For a pure l=1 dipole the radial field must scale as (r1/r2)^(l+2) = (r1/r2)^3
    between any two radii.  This is exact (not statistical) so the tolerance is tight."""
    lmax = 1
    glm, hlm, idx = _make_single_mode(lmax, 1, 0)
    rplanet = 1.0
    r1, r2 = 1.0, 2.0
    nphi = 16

    br1, bt1, _ = extrapot_scipy(glm, hlm, idx, lmax, lmax, rplanet,
                                  np.array([r1]), nphi=nphi)
    br2, bt2, _ = extrapot_scipy(glm, hlm, idx, lmax, lmax, rplanet,
                                  np.array([r2]), nphi=nphi)

    expected_ratio = (r1 / r2) ** 3
    np.testing.assert_allclose(br2[..., 0], br1[..., 0] * expected_ratio, rtol=1e-12)
    np.testing.assert_allclose(bt2[..., 0], bt1[..., 0] * expected_ratio, rtol=1e-12)


# ---------------------------------------------------------------------------
# extrapot_scipy – axisymmetric field (mmax=0) has no phi-component
# ---------------------------------------------------------------------------

def test_extrapot_scipy_bphi_zero_for_mmax0():
    """An axisymmetric field (mmax=0) must produce bphi=0 everywhere."""
    lmax = 3
    glm, hlm, idx = _make_single_mode(lmax, 1, 0)
    _, _, bp = extrapot_scipy(glm, hlm, idx, lmax, 0, 1.0,
                               np.array([1.5]), nphi=16)
    np.testing.assert_allclose(bp, 0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# get_field_along_path_scipy – shape-mismatch guard
# ---------------------------------------------------------------------------

def test_get_field_along_path_scipy_shape_mismatch():
    """Mismatched r/theta/phi arrays must raise ValueError."""
    lmax = 1
    glm, hlm, idx = _make_single_mode(lmax, 1, 0)
    with pytest.raises(ValueError):
        get_field_along_path_scipy(
            glm, hlm, idx, lmax,
            np.array([1.0, 1.5]),
            np.array([np.pi / 2]),
            np.array([0.0]),
        )


@pytest.mark.parametrize("planet", ["earth", "jupiter", "saturn"])
def test_extrapot_scipy_matches_surface_field(planet):
    """At r=1 the extrapolated radial field must equal getB, including the
    Condon-Shortley phase convention for Earth."""
    p = Planet(name=planet, nphi=32, info=False, units='nT')
    br, _, _ = extrapot_scipy(p.glm, p.hlm, p.idx, p.lmax, p.mmax, 1.0,
                              np.array([1.0]), nphi=32, planetname=planet)
    np.testing.assert_allclose(br[..., 0], p.Br, atol=1e-8 * np.abs(p.Br).max())


@pytest.mark.parametrize("planet", ["earth", "jupiter", "saturn"])
def test_orbit_matches_grid(planet):
    p = Planet(name=planet, nphi=64, info=False)
    p.extrapolate([2.0])
    p.orbit_path([2.0], [p.theta[7]], [p.phi[5]])
    np.testing.assert_allclose(p.br_orb[0], p.br_ex[5, 7, 0], rtol=1e-6)
    np.testing.assert_allclose(p.btheta_orb[0], p.btheta_ex[5, 7, 0], rtol=1e-6)
    np.testing.assert_allclose(p.bphi_orb[0], p.bphi_ex[5, 7, 0], rtol=1e-6, atol=1e-9)
