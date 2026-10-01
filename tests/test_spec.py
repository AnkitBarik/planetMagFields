import numpy as np
import pytest
from planetmagfields.libgauss import gen_idx, get_spec


def _make_single_mode(lmax, l1, m1, g_val=1.0, h_val=0.0):
    ncoeff = (lmax + 1) * (lmax + 2) // 2
    idx    = gen_idx(lmax)
    glm    = np.zeros(ncoeff)
    hlm    = np.zeros(ncoeff)
    glm[idx[l1, m1]] = g_val
    hlm[idx[l1, m1]] = h_val
    return glm, hlm, idx


# Equatorially symmetric modes have l-m even, anti-symmetric modes l-m odd
@pytest.mark.parametrize("mmax", [0, 4])
@pytest.mark.parametrize("l1, m1, symmetric", [
    (1, 0, False),
    (2, 0, True),
    (3, 0, False),
    (4, 0, True),
])
def test_axisymmetric_mode_symmetry(mmax, l1, m1, symmetric):
    lmax = 4
    glm, hlm, idx = _make_single_mode(lmax, l1, m1)
    E, _, E_symm, E_antisymm, _ = get_spec(glm, hlm, idx, lmax, mmax)

    if symmetric:
        assert E_symm == pytest.approx(E[l1])
        assert E_antisymm == 0.
    else:
        assert E_antisymm == pytest.approx(E[l1])
        assert E_symm == 0.


@pytest.mark.parametrize("l1, m1, symmetric", [
    (1, 1, True),
    (2, 1, False),
    (3, 1, True),
    (3, 2, False),
])
def test_nonaxisymmetric_mode_symmetry(l1, m1, symmetric):
    lmax = 4
    glm, hlm, idx = _make_single_mode(lmax, l1, m1, h_val=0.5)
    E, _, E_symm, E_antisymm, _ = get_spec(glm, hlm, idx, lmax, lmax)

    if symmetric:
        assert E_symm == pytest.approx(E[l1])
        assert E_antisymm == 0.
    else:
        assert E_antisymm == pytest.approx(E[l1])
        assert E_symm == 0.


def test_branches_agree():
    """mmax=0 and general branches must give identical results for an
    axisymmetric field."""
    lmax = 6
    ncoeff = (lmax + 1) * (lmax + 2) // 2
    idx = gen_idx(lmax)
    rng = np.random.default_rng(0)
    glm = np.zeros(ncoeff)
    hlm = np.zeros(ncoeff)
    for l in range(1, lmax + 1):
        glm[idx[l, 0]] = rng.normal()

    res0 = get_spec(glm, hlm, idx, lmax, 0, r=0.8)
    resg = get_spec(glm, hlm, idx, lmax, lmax, r=0.8)
    for a, b in zip(res0, resg):
        np.testing.assert_allclose(a, b)
