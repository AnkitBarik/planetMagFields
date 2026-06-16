#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from .libgauss import get_grid, getB
from .potextra import extrapot
from .utils import is_dark_color

def get_grid3d(r,theta,phi):
    """
    Produces a 3D grid for storing values of extrapolated field. Used in producing
    vtk file for 3D visualization.

    Parameters
    ----------
    r : array_like
        Array of radial levels
    theta : array_like
        Array of colatitudes, ranging from 0 to pi
    phi : array_like
        Array of longitudes, ranging from 0 to 2*pi

    Returns
    -------
    r3D : ndarray(float, ndim=3)
        3D array of radius values at each grid point, shape (nphi,ntheta,nr)
    th3D : ndarray(float, ndim=3)
        3D array of colatitude values at each grid point, shape (nphi,ntheta,nr)
    p3D : ndarray(float, ndim=3)
        3D array of longitude values at each grid point, shape (nphi,ntheta,nr)
    x : ndarray(float, ndim=3)
        3D array of Cartesian x values at each grid point, shape (nphi,ntheta,nr)
    y : ndarray(float, ndim=3)
        3D array of Cartesian y values at each grid point, shape (nphi,ntheta,nr)
    z : ndarray(float, ndim=3)
        3D array of Cartesian z values at each grid point, shape (nphi,ntheta,nr)
    s : ndarray(float, ndim=3)
        3D array of cylindrical radius values at each grid point, shape (nphi,ntheta,nr)
    """

    p3D, th3D, r3D = np.meshgrid(phi, theta, r, indexing='ij')

    s = r3D * np.sin(th3D)
    x = s *   np.cos(p3D)
    y = s *   np.sin(p3D)
    z = r3D * np.cos(th3D)

    return r3D,th3D,p3D, x,y,z

def get_cart(vr,vt,vp,th3D,p3D):
    """
    Converts a set of vector components in spherical coordinate system
    to Cartesian. All components, both at input and output have shape
    (nphi,ntheta,nr), number of points in longitude, colatitude and radial
    directions, respectively.

    Parameters
    ----------
    vr : ndarray(float, ndim=3)
        3D array of radial component of vector
    vt : ndarray(float, ndim=3)
        3D array of colatitude component of vector
    vp : ndarray(float, ndim=3)
        3D array of azimuthal component of vector
    th3D : ndarray(float, ndim=3)
        3D array of colatitude values on the grid
    p3D : ndarray(float, ndim=3)
        3D array of longitude values on the grid

    Returns
    -------
    vx : ndarray(float, ndim=3)
        3D array of vector component in Cartesian x direction
    vy : ndarray(float, ndim=3)
        3D array of vector component in Cartesian x direction
    vz : ndarray(float, ndim=3)
        3D array of vector component in Cartesian x direction
    """

    vs = vr * np.sin(th3D) + vt * np.cos(th3D)
    vz = vr * np.cos(th3D) - vt * np.sin(th3D)

    vx = vs * np.cos(p3D) - vp * np.sin(p3D)
    vy = vs * np.sin(p3D) + vp * np.cos(p3D)

    return vx,vy,vz


def writeVts(name,br,btheta,bphi,r,theta,phi,r_planet=1):
    """
    Writes an unstructured vtk file for 3D visualization.

    Parameters
    ----------
    name : str
        Name of the file
    br : ndarray(float, ndim=3)
        3D array of radial component of vector
    btheta : ndarray(float, ndim=3)
        3D array of colatitudinal component of vector
    bphi : ndarray(float, ndim=3)
        3D array of azimuthal component of vector
    r : array_like
        Array of radial levels
    theta : array_like
        Array of colatitudes
    phi : array_like
        Array of longitudes

    Returns
    -------
    None
    """

    r3D,th3D,p3D, x,y,z = get_grid3d(r*r_planet,theta,phi)

    print("grid shape=",th3D.shape)

    bx,by,bz = get_cart(br, btheta, bphi,th3D,p3D)

    try:
        from pyevtk.hl import gridToVTK
    except:
        print("This requires the use of pyevtk library!")
        print("You can get it from https://github.com/paulo-herrera/PyEVTK")

    br = np.asfortranarray(br)
    bx = np.asfortranarray(bx)
    by = np.asfortranarray(by)
    bz = np.asfortranarray(bz)

    gridToVTK("%s"%name,x,y,z,pointData= {"radius":r3D,
                                          "Radial mag field":br,
                                          "Mag Field":(bx, by,bz)})

    print("Output written to %s!" %name)

# PyVista plotting
##################

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

def set_plotter_properties(plotter, bgcolor):
    plotter.set_background(bgcolor)

    if is_dark_color(bgcolor):
        font_color = 'white'
        # Add axes with white labels
        plotter.add_axes(
            line_width=2,
            color='white',           # This sets the axes line color
            x_color='white',         # X-axis color (line and label)
            y_color='white',         # Y-axis color (line and label)
            z_color='white',         # Z-axis color (line and label)
            xlabel='X',
            ylabel='Y',
            zlabel='Z',
            label_size=(0.1, 0.1),   # Label size as fraction of viewport
            labels_off=False,
        )
    else:
        font_color = 'black'
        plotter.add_axes(line_width=2)

    return font_color

def render_tight(plotter, mesh, padding=0.1):
    """Render with tight bounds around the mesh."""

    # Get mesh bounds
    bounds = mesh.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)

    # Calculate center and size
    center = mesh.center
    size = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4]
    )

    # Set camera to fit mesh tightly
    plotter.camera.focal_point = center
    plotter.camera.position = (
        center[0] + size * (1 + padding),
        center[1] + size * (1 + padding),
        center[2] + size * (1 + padding)
    )

    # Reset camera to fit all actors with some padding
    plotter.reset_camera(bounds=bounds)
    plotter.camera.zoom(0.8)

    return plotter

def plot_surface(theta,phi,dat,fieldname='Br',cmap='seismic',clim_fac=1.0, bgcolor='black'):

    th_bounds = _cell_bounds(theta * 180/np.pi)
    ph_bounds = _cell_bounds(phi * 180/np.pi)

    try:
        import pyvista as pv
    except ImportError:
        print("This requires the use of pyvista library!")
        print("You can install it with pip install pyvista")

    grid = pv.grid_from_sph_coords(ph_bounds,th_bounds,r=1)

    grid.cell_data[fieldname] = dat.ravel('C')

    pl = pv.Plotter(window_size=(800, 800))

    font_color = set_plotter_properties(pl, bgcolor)

    datMax = np.max(np.abs(dat))
    clim = [-clim_fac * datMax, clim_fac * datMax]
    pl.add_mesh(grid, cmap=cmap,clim=clim,smooth_shading=True, show_scalar_bar=False)

    pl = render_tight(pl, grid, padding=0.1)

    return pl, font_color

def render_field_lines(planetname, glm, hlm, idx, lmax, mmax, rplanet,
                       rout, nphi=128, surf=False, clim_fac=1.0,
                       units='nT', bgcolor='black', cmap='seismic',
                       lightweight=False):

    """
    Renders field lines using pyvista.

    Parameters
    ----------
    planetname : str
        Name of the planet, used for title of the plot
    glm : array_like
        Array of Gauss coefficients, shape (lmax+1,mmax+1)
    hlm : array_like
        Array of harmonic coefficients, shape (lmax+1,mmax+1)
    idx : int
        Index of the radial level to plot
    lmax : int
        Maximum degree of spherical harmonic expansion
    mmax : int
        Maximum order of spherical harmonic expansion
    rplanet : float
        Radius of the planet, used for scaling the plot
    rout : float
        Outer radius of the grid, used for scaling the plot
    nphi : int, optional
        Number of points in longitude direction. If None, it will be set to 2*mmax+1

    Returns
    -------
    pl : pyvista.Plotter
        PyVista plotter object with the rendered field lines.
    """

    p2D, th2D, phi, theta = get_grid(nphi, nphi//2)

    if units.lower() == 'mut':
        unitfac = 1e-3
        cbar_title = r'$B_r (\mu$T)'
    elif units.lower() == 'nt':
        unitfac = 1.
        cbar_title = r'$B_r$ (nT)'
    elif units.lower() == 'gauss':
        unitfac = 1e-5
        cbar_title = r'$B_r$ (Gauss)'

    glm = glm * unitfac
    hlm = hlm * unitfac

    brout, btout, bpout = extrapot(planetname, glm, hlm, idx, lmax, mmax, rplanet,
                                   np.atleast_1d(rout * rplanet), nphi=nphi)

    r = np.atleast_1d(rout * rplanet)

    try:
        import pyvista as pv
    except ImportError:
        print("This requires the use of pyvista library!")
        print("You can install it with pip install pyvista")

    grid = pv.grid_from_sph_coords(phi * 180/np.pi, theta * 180/np.pi, r)

    clim = [-clim_fac * np.max(np.abs(brout)), clim_fac * np.max(np.abs(brout))]

    if surf:
        br_surf = getB(lmax,mmax,glm,hlm,idx,1,p2D,th2D,planetname=planetname)
        pl, font_color = plot_surface(theta, phi, br_surf, clim_fac=clim_fac, bgcolor=bgcolor)
    else:
        pl = pv.Plotter(window_size=(800, 800))
        font_color = set_plotter_properties(pl, bgcolor)

    # Transform vectors to cartesian coordinates
    _, th3d, ph3d, _, _, _ = get_grid3d(r*rplanet, theta, phi)
    bx, by, bz = get_cart(brout, btout, bpout, th3d, ph3d)

    bx = np.ascontiguousarray(np.transpose(bx,[2,0,1]), dtype=np.float32)
    by = np.ascontiguousarray(np.transpose(by,[2,0,1]), dtype=np.float32)
    bz = np.ascontiguousarray(np.transpose(bz,[2,0,1]), dtype=np.float32)
    br = np.ascontiguousarray(np.transpose(brout,[2,0,1]), dtype=np.float32)

    vectors = np.stack((bx, by, bz), axis=-1)

    grid.point_data['B'] = vectors.reshape(-1, 3, order='C')
    grid.point_data['Br'] = br.ravel(order='C')
    grid.point_data['Energy'] = 0.5 * (bx**2 + by**2 + bz**2).ravel(order='C')

    grid.set_active_vectors('B')

    max_length = rout.max() * rplanet * 2.0
    initial_step = 0.2
    max_steps = 1000

    n_seeds = 30 if lightweight else 100
    n_steps = 500 if lightweight else 1000

    streamlines = grid.streamlines(
            vectors='B',
            n_points=n_seeds,
            source_radius=r.max(),
            source_center=(0, 0, 0),
            integration_direction='both',
            max_length=max_length,
            initial_step_length=initial_step,
            max_steps=n_steps,
            terminal_speed=1e-12,
            interpolator_type='cell'
        )

    if streamlines.n_points > 0:
        # Get and normalize energy for tube radius scaling
        b_mag_sq_vals = streamlines.point_data['Energy']
        max_b_sq = np.max(b_mag_sq_vals)

        if max_b_sq > 0:
            radius_scale = b_mag_sq_vals / max_b_sq
            # np.clip(radius_scale, 0.1, 1.0, out=radius_scale)
        else:
            radius_scale = np.ones_like(b_mag_sq_vals) * 0.5

        streamlines.point_data['radius_scale'] = radius_scale

        scalar_bar_args = {'title': cbar_title,
                           'vertical': True,
                           'position_x': 0.88,
                           'position_y': 0.2,
                           'width': 0.08,
                           'height': 0.6,
                           'title_font_size': 25,
                           'label_font_size': 18,
                           'n_labels': 5,
                           'fmt': '%.1e',
                           'font_family': 'times',
                           'color': font_color}

        if lightweight:
            # Render as lines — much faster than tubes
            pl.add_mesh(
                streamlines,
                scalars='Br',
                cmap=cmap,
                clim=clim,
                line_width=1.5,
                show_scalar_bar=True,
                scalar_bar_args=scalar_bar_args
            )
        else:
            tubes = streamlines.tube(
                scalars='radius_scale',
                radius=rplanet * 0.015,
                radius_factor=8,
                n_sides=8,
                capping=False
            )
            pl.add_mesh(
                tubes,
                scalars='Br',
                cmap=cmap,
                clim=clim,
                smooth_shading=True,
                show_scalar_bar=True,
                scalar_bar_args=scalar_bar_args
            )

    pl = render_tight(pl, grid, padding=0.1)

    return pl
