import numpy as np
import os
import sys
import pytest
sys.path.insert(0, os.path.abspath('../'))
from planetmagfields import Planet
from planetmagfields.libgauss import gen_idx, get_grid
from planetmagfields.potextra import (
    get_pol_from_Gauss,
    extrapot_scipy,
    get_field_along_path_scipy
)


class TestJupiterBr:
    def test_surface_br(self):
        p = Planet(name='jupiter', r=0.85, nphi=256, info=False, model='jrm33')

        br_ref = np.loadtxt('jupiter/Br_reference085.dat')
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

        br_ref = np.loadtxt('jupiter/Br_reference.dat')
        bt_ref = np.loadtxt('jupiter/Bt_reference.dat')
        bp_ref = np.loadtxt('jupiter/Bp_reference.dat')

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
    @staticmethod
    def _make_single_mode(lmax, l1, m1, g_val=1.0, h_val=0.0):
        ncoeff = (lmax + 1) * (lmax + 2) // 2
        idx    = gen_idx(lmax)
        glm    = np.zeros(ncoeff)
        hlm    = np.zeros(ncoeff)
        glm[idx[l1, m1]] = g_val
        hlm[idx[l1, m1]] = h_val
        return glm, hlm, idx

    def test_extrapot_scipy_vs_field_along_path(self):
        """extrapot_scipy (full 3-D grid) and get_field_along_path_scipy
        (arbitrary points) must agree to floating-point precision."""
        lmax = 3
        glm, hlm, idx = self._make_single_mode(lmax, 2, 1, g_val=1.0, h_val=0.5)
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

    def test_get_cart_spherical_basis_vectors(self):
        """At (theta=pi/2, phi=0): r-hat -> +x, theta-hat -> -z, phi-hat -> +y."""
        from planetmagfields.lib3d import get_cart

        th    = np.array([[[np.pi / 2]]])
        p     = np.array([[[0.0]]])
        ones  = np.ones((1, 1, 1))
        zeros = np.zeros((1, 1, 1))

        vx, vy, vz = get_cart(ones, zeros, zeros, th, p)
        np.testing.assert_allclose(
            [vx[0, 0, 0], vy[0, 0, 0], vz[0, 0, 0]], [1, 0, 0], atol=1e-14
        )

        vx, vy, vz = get_cart(zeros, ones, zeros, th, p)
        np.testing.assert_allclose(
            [vx[0, 0, 0], vy[0, 0, 0], vz[0, 0, 0]], [0, 0, -1], atol=1e-14
        )

        vx, vy, vz = get_cart(zeros, zeros, ones, th, p)
        np.testing.assert_allclose(
            [vx[0, 0, 0], vy[0, 0, 0], vz[0, 0, 0]], [0, 1, 0], atol=1e-14
        )


class TestGetPolFromGauss:
    @staticmethod
    def _make_single_mode(lmax, l1, m1, g_val=1.0, h_val=0.0):
        ncoeff = (lmax + 1) * (lmax + 2) // 2
        idx    = gen_idx(lmax)
        glm    = np.zeros(ncoeff)
        hlm    = np.zeros(ncoeff)
        glm[idx[l1, m1]] = g_val
        hlm[idx[l1, m1]] = h_val
        return glm, hlm, idx

    def test_m0_value(self):
        """For g(1,0)=1 the poloidal coefficient at (l=1,m=0) must equal
        sqrt(4*pi/3)/l; all other modes must be zero."""
        lmax = 2
        glm, hlm, idx = self._make_single_mode(lmax, 1, 0)
        bpol = get_pol_from_Gauss('earth', glm, hlm, lmax, lmax, idx)

        expected = np.sqrt(4 * np.pi / 3) / 1
        np.testing.assert_allclose(bpol[idx[1, 0]], expected, rtol=1e-14)
        bpol[idx[1, 0]] = 0.0
        np.testing.assert_allclose(bpol, 0.0, atol=1e-14)

    def test_fac_m_earth_vs_nonearth(self):
        """For m>0, earth uses fac_m=1 while other planets use (-1)^m."""
        lmax = 2
        glm, hlm, idx = self._make_single_mode(lmax, 1, 1)
        norm_m1 = np.sqrt(2 * np.pi / 3) / 1

        bpol_earth = get_pol_from_Gauss('earth',   glm, hlm, lmax, lmax, idx)
        bpol_jup   = get_pol_from_Gauss('jupiter', glm, hlm, lmax, lmax, idx)

        np.testing.assert_allclose(bpol_earth[idx[1, 1]],  norm_m1, rtol=1e-14)
        np.testing.assert_allclose(bpol_jup[idx[1, 1]],   -norm_m1, rtol=1e-14)

    def test_mmax0_ignores_nonzero_m(self):
        """When mmax=0, modes with m>0 must not contribute."""
        lmax = 3
        glm, hlm, idx = self._make_single_mode(lmax, 2, 1, g_val=5.0, h_val=3.0)
        bpol = get_pol_from_Gauss('earth', glm, hlm, lmax, 0, idx)

        np.testing.assert_allclose(bpol, 0.0, atol=1e-14)

    p = Planet(name='jupiter',r=0.85,nphi=256,info=False,model='jrm33')

    br_ref = np.loadtxt('jupiter/Br_reference085.dat')
    br = p.Br

    percent_err = np.abs( (br_ref - br)/br_ref ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=0.1,
                                atol=0.1)

def test_potextra_internal():
    p = Planet(name='jupiter',r=10,nphi=256,info=False,model='jrm33')

    p.extrapolate(np.array([10]))
    br_ref = np.squeeze(p.br_ex)
    err = np.abs( (br_ref - p.Br) )

    np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)

def test_potextra_jupiterMag():

    p = Planet(name='jupiter',nphi=512,info=False,model='jrm33')
    p.extrapolate([2])

    br_ref = np.loadtxt('jupiter/Br_reference.dat')
    bt_ref = np.loadtxt('jupiter/Bt_reference.dat')
    bp_ref = np.loadtxt('jupiter/Bp_reference.dat')

    percent_err = ( np.mean(np.abs(p.br_ex[...,0]*1e3 - br_ref)/br_ref)
                +  np.mean(np.abs(p.btheta_ex[...,0]*1e3 - bt_ref)/bt_ref)
                +  np.mean(np.abs(p.bphi_ex[...,0]*1e3 - bp_ref)/bp_ref)
                ) * 100

    np.testing.assert_allclose(percent_err, 0, rtol=1,
                            atol=1)

def test_potextra_internal_m0():
    p = Planet(name='saturn',r=10,nphi=256,info=False,model='cassini11+')

    p.extrapolate(np.array([10]))
    br_ref = np.squeeze(p.br_ex)
    err = np.abs( (br_ref - p.Br) )

    np.testing.assert_allclose(err, 0, rtol=1e-2, atol=1e-2)

def test_extrapot_scipy_vs_field_along_path():
    """extrapot_scipy (full 3-D grid) and get_field_along_path_scipy (arbitrary
    points) are independent implementations of the same formula.  They must
    agree to floating-point precision at matching (r, theta, phi) coordinates."""
    lmax = 3
    glm, hlm, idx = _make_single_mode(lmax, 2, 1, g_val=1.0, h_val=0.5)
    rplanet = 1.0
    r_test  = 1.5
    nphi   = 32
    ntheta = nphi // 2

    br_grid, bt_grid, bp_grid = extrapot_scipy(
        glm, hlm, idx, lmax, lmax, rplanet, np.array([r_test]), nphi=nphi
    )

    # Reconstruct the exact grid used internally by extrapot_scipy
    _,_,phi,theta = get_grid(nphi, ntheta)

    for i_phi, i_theta in [(5, 4), (10, 7), (20, 12)]:
        br_pt, bt_pt, bp_pt = get_field_along_path_scipy(
            glm, hlm, idx, lmax,
            np.array([r_test]),
            np.array([theta[i_theta]]),
            np.array([phi[i_phi]]),
        )
        np.testing.assert_allclose(br_grid[i_phi, i_theta, 0], br_pt[0], rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(bt_grid[i_phi, i_theta, 0], bt_pt[0], rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(bp_grid[i_phi, i_theta, 0], bp_pt[0], rtol=1e-10, atol=1e-10)


def test_orbit():
   p = Planet(name='earth',r=1,nphi=256,info=False,model='igrf14')
   p.extrapolate([2])
   p.orbit_path([2],[p.theta[10]],[p.phi[10]])

   br_ref,bt_ref,bp_ref = p.br_ex[10,10,0],p.btheta_ex[10,10,0],p.bphi_ex[10,10,0]

   err = np.sqrt((p.br_orb-br_ref)**2
                 +(p.btheta_orb-bt_ref)**2
                 +(p.bphi_orb-bp_ref)**2)

   np.testing.assert_allclose(err, 0, rtol=1e-3, atol=1e-3)


def test_orbit_m0():
   p = Planet(name='saturn',r=1,nphi=256,info=False,model='cassini11+')
   p.extrapolate([2])
   p.orbit_path([2],[p.theta[10]],[p.phi[10]])

   br_ref,bt_ref,bp_ref = p.br_ex[10,10,0],p.btheta_ex[10,10,0],p.bphi_ex[10,10,0]

   err = np.sqrt((p.br_orb-br_ref)**2
                 +(p.btheta_orb-bt_ref)**2
                 +(p.bphi_orb-bp_ref)**2)

   np.testing.assert_allclose(err, 0, rtol=1e-4, atol=1e-4)


# ---------------------------------------------------------------------------
# Helpers shared by the unit tests below
# ---------------------------------------------------------------------------

def _make_single_mode(lmax, l1, m1, g_val=1.0, h_val=0.0):
    """Return (glm, hlm, idx) for a field with exactly one non-zero Gauss mode."""
    ncoeff = (lmax + 1) * (lmax + 2) // 2
    idx    = gen_idx(lmax)
    glm    = np.zeros(ncoeff)
    hlm    = np.zeros(ncoeff)
    glm[idx[l1, m1]] = g_val
    hlm[idx[l1, m1]] = h_val
    return glm, hlm, idx


# ---------------------------------------------------------------------------
# get_pol_from_Gauss
# ---------------------------------------------------------------------------

def test_get_pol_from_Gauss_m0_value():
    """For g(1,0)=1, the poloidal coefficient at (l=1,m=0) must equal
    sqrt(4*pi/3)/l (norm * 1 * (1 - 0j))."""
    lmax = 2
    glm, hlm, idx = _make_single_mode(lmax, 1, 0)
    bpol = get_pol_from_Gauss('earth', glm, hlm, lmax, lmax, idx)
    expected = np.sqrt(4 * np.pi / 3) / 1
    np.testing.assert_allclose(bpol[idx[1, 0]], expected, rtol=1e-14)
    # All other modes must remain zero
    bpol[idx[1, 0]] = 0.0
    np.testing.assert_allclose(bpol, 0.0, atol=1e-14)


def test_get_pol_from_Gauss_fac_m_earth_vs_nonearth():
    """For m>0 modes, earth always uses fac_m=1 while other planets use (-1)^m.
    For l=1, m=1, g=1: bpol(earth) = +norm, bpol(jupiter) = -norm."""
    lmax = 2
    glm, hlm, idx = _make_single_mode(lmax, 1, 1)
    norm_m1 = np.sqrt(2 * np.pi / 3) / 1   # sqrt(2*pi/(2*1+1)) / l with l=1
    bpol_earth = get_pol_from_Gauss('earth',   glm, hlm, lmax, lmax, idx)
    bpol_jup   = get_pol_from_Gauss('jupiter', glm, hlm, lmax, lmax, idx)
    np.testing.assert_allclose(bpol_earth[idx[1, 1]],  norm_m1, rtol=1e-14)
    np.testing.assert_allclose(bpol_jup[idx[1, 1]],   -norm_m1, rtol=1e-14)


def test_get_pol_from_Gauss_mmax0_ignores_nonzero_m():
    """When mmax=0, modes with m>0 must not contribute even if g/h are nonzero."""
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