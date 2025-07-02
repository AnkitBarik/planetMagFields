#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from .potextra import

def _cell_bounds(points, bound_position=0.5):
    """
    Calculate coordinate cell boundaries. Copied from:
    https://docs.pyvista.org/examples/02-plot/spherical.html

    Parameters
    ----------
    points: numpy.ndarray
        One-dimensional array of uniformly spaced values of shape (M,).

    bound_position: bool, optional
        The desired position of the bounds relative to the position
        of the points.

    Returns
    -------
    bounds: numpy.ndarray
        Array of shape (M+1,)

    Examples
    --------
    >>> a = np.arange(-1, 2.5, 0.5)
    >>> a
    array([-1. , -0.5,  0. ,  0.5,  1. ,  1.5,  2. ])
    >>> cell_bounds(a)
    array([-1.25, -0.75, -0.25,  0.25,  0.75,  1.25,  1.75,  2.25])

    """
    if points.ndim != 1:
        msg = 'Only 1D points are allowed.'
        raise ValueError(msg)
    diffs = np.diff(points)
    delta = diffs[0] * bound_position
    return np.concatenate([[points[0] - delta], points + delta])

def plot_surface(theta,phi,dat,fieldname='Br',cmap='seismic'):

    th_bounds = _cell_bounds(theta * 180/np.pi)
    ph_bounds = _cell_bounds(phi * 180/np.pi)

    p2d,th2d = np.meshgrid(phi,theta,indexing='ij')

    grid = pv.grid_from_sph_coords(ph_bounds,th_bounds,r=1)

    grid.cell_data[fieldname] = dat.ravel('C')

    pl = pv.Plotter()
    pl.add_mesh(grid, cmap='seismic',smooth_shading=True)
    pl.show()