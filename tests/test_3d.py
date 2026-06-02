#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from planetmagfields.lib3d import get_cart, get_grid3d


def test_get_cart_spherical_basis_vectors():
    """get_cart must map each unit spherical basis vector to the correct
    Cartesian direction.  At (theta=pi/2, phi=0):
      r-hat     -> +x-hat
      theta-hat -> -z-hat
      phi-hat   -> +y-hat
    """

    th    = np.array([[[np.pi / 2]]])
    p     = np.array([[[0.0]]])
    ones  = np.ones((1, 1, 1))
    zeros = np.zeros((1, 1, 1))

    vx, vy, vz = get_cart(ones, zeros, zeros, th, p)
    np.testing.assert_allclose(
        [vx[0,0,0], vy[0,0,0], vz[0,0,0]], [1, 0, 0], atol=1e-14
    )

    vx, vy, vz = get_cart(zeros, ones, zeros, th, p)
    np.testing.assert_allclose(
        [vx[0,0,0], vy[0,0,0], vz[0,0,0]], [0, 0, -1], atol=1e-14
    )

    vx, vy, vz = get_cart(zeros, zeros, ones, th, p)
    np.testing.assert_allclose(
        [vx[0,0,0], vy[0,0,0], vz[0,0,0]], [0, 1, 0], atol=1e-14
    )

# ---------------------------------------------------------------------------
# get_grid – Cartesian coordinates at a known point
# ---------------------------------------------------------------------------

def test_get_grid_cartesian_equator():
    """At r=1, theta=pi/2 (equator), phi=0: x=1, y=0, z=0."""
    r3D, th3D, p3D, x, y, z = get_grid3d(
        np.array([1.0]),
        np.array([np.pi / 2]),
        np.array([0.0]),
    )
    np.testing.assert_allclose(x[0, 0, 0], 1.0, atol=1e-14)
    np.testing.assert_allclose(y[0, 0, 0], 0.0, atol=1e-14)
    np.testing.assert_allclose(z[0, 0, 0], 0.0, atol=1e-14)


def test_get_grid_cartesian_pole():
    """At r=1, theta=0 (north pole), any phi: x=0, y=0, z=1."""
    r3D, th3D, p3D, x, y, z = get_grid3d(
        np.array([1.0]),
        np.array([0.0]),
        np.array([1.23]),   # phi irrelevant at pole
    )
    np.testing.assert_allclose(x[0, 0, 0], 0.0, atol=1e-14)
    np.testing.assert_allclose(y[0, 0, 0], 0.0, atol=1e-14)
    np.testing.assert_allclose(z[0, 0, 0], 1.0, atol=1e-14)