#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def _import_cartopy():
    try:
        import cartopy.crs as ccrs
        return ccrs
    except ImportError:
        print("cartopy library not available, using Hammer projection")
        return None


def radius_label(r):
    if r == 1:
        return '  Surface'
    return r'  $r/r_{\rm surface}=%.2f$' %r


def _range_label(sym, arr, vmin, vmax, lmax):
    if arr is not None:
        return r', $%s = %s$' %(sym, str(arr))
    if vmin > 0:
        if vmax < lmax:
            return r', $ %d \leq %s \leq %d$' %(vmin, sym, vmax)
        return r', $%s \geq %d$' %(sym, vmin)
    if vmax < lmax:
        return r', $%s \leq %d$' %(sym, vmax)
    return ''


def filter_label(lmax, larr=None, marr=None, lCutMin=0, lCutMax=None,
                 mmin=0, mmax=None):
    """Title suffix describing a spherical harmonic filter."""
    lCutMax = lmax if lCutMax is None else lCutMax
    mmax = lmax if mmax is None else mmax
    return (_range_label('l', larr, lCutMin, lCutMax, lmax)
            + _range_label('m', marr, mmin, mmax, lmax))


def get_color_limits(dat,vmin=None,vmax=None):
    """Computes minimum and maximum of colorbar for plots

    Parameters
    ----------
    dat : ndarray(float, ndim=2)
        Data for plotting
    vmin : float, optional
        Minimum of colorscale, by default None
    vmax : float, optional
        Maximum of colorscale, by default None

    Returns
    -------
    vmin : float, optional
        Minimum of colorscale
    vmax : float, optional
        Minimum of colorscale
    """
    if vmin is not None and vmax is not None:
        return vmin,vmax
    elif vmin is None and vmax is not None:
        vmin = -vmax
    elif vmin is not None and vmax is None:
        vmax = -vmin
    else:
        bmax = np.abs(dat).max()
        digits = int(np.log10(bmax)) + 1

        if digits > 1:
            bmax = np.round(bmax)
        else:
            bmax = np.round(bmax,decimals=1)

        vmin = -bmax
        vmax = bmax
    return vmin, vmax

def lonlat_grid(p2D, th2D, B):
    """
    Convert a (phi, theta) grid to longitude/latitude arrays for map plotting.

    phi is east longitude, so longitude is phi wrapped into [-pi, pi]. Columns
    are reordered so longitude increases monotonically, and exact columns at
    -pi and +pi (interpolated across the seam) are added so filled contours
    cover the whole map without a gap.

    Parameters
    ----------
    p2D : ndarray(float, ndim=2)
        Longitude (phi) on a (nphi, ntheta) grid, ranging over [0, 2*pi]
    th2D : ndarray(float, ndim=2)
        Co-latitude on the same grid
    B : ndarray(float, ndim=2)
        Data on the same grid

    Returns
    -------
    lon2D, lat2D, B : ndarray(float, ndim=2)
        Longitude and latitude in radians and the reordered data
    """
    phi = p2D[:, 0]
    th = th2D[0, :]
    # Drop the duplicated periodic end point (phi = 2 pi) if present
    if np.isclose(phi[-1] - phi[0], 2*np.pi):
        phi, B = phi[:-1], B[:-1]
    lon = (phi + np.pi) % (2*np.pi) - np.pi
    order = np.argsort(lon, kind='stable')
    lon, B = lon[order], B[order]
    # Linear interpolation across the seam to get the values at exactly +-pi
    gap = lon[0] + 2*np.pi - lon[-1]
    w = (np.pi - lon[-1]) / gap if gap > 0 else 0.
    B_seam = (1 - w) * B[-1] + w * B[0]
    lon = np.concatenate([[-np.pi], lon, [np.pi]])
    B = np.concatenate([B_seam[None], B, B_seam[None]])
    lon2D, th2D = np.meshgrid(lon, th, indexing='ij')
    return lon2D, np.pi/2 - th2D, B


def hammer2cart(ttheta, pphi, colat=False):
    """
    This function is used to define the Hammer projection for
    default plotting. Copied from MagIC python plotting script:

    https://github.com/magic-sph/magic/blob/master/python/magic/plotlib.py

    Parameters
    ----------
    ttheta : ndarray(float, ndim=2)
        2D array defining latitude or co-latitude (theta).
        This ranges from 0 to pi and has a shape (nphi,ntheta)
    pphi : ndarray(float, ndim=2)
        2D array defining longitude (phi)
        This ranges from 0 to 2*pi and has a shape (nphi,ntheta)
    colat : Boolean
        Flag defining whether ttheta is colatitude (True) or latitude(False)
    """

    if not colat: # for lat and phi \in [-pi, pi]
        xx = 2.*np.sqrt(2.) * np.cos(ttheta)*np.sin(pphi/2.)\
             /np.sqrt(1.+np.cos(ttheta)*np.cos(pphi/2.))
        yy = np.sqrt(2.) * np.sin(ttheta)\
             /np.sqrt(1.+np.cos(ttheta)*np.cos(pphi/2.))
    else:  # for colat and phi \in [0, 2pi]
        xx = -2.*np.sqrt(2.) * np.sin(ttheta)*np.cos(pphi/2.)\
             /np.sqrt(1.+np.sin(ttheta)*np.sin(pphi/2.))
        yy = np.sqrt(2.) * np.cos(ttheta)\
             /np.sqrt(1.+np.sin(ttheta)*np.sin(pphi/2.))
    return xx, yy

def plotSurf(p2D,th2D,B,levels=60,cmap='RdBu_r',
             proj='Mollweide',vmin=None,vmax=None):
    """
    Plots magnetic field on a surface defined by 2D arrays
    of longitude and co-latitude.

    Parameters
    ----------
    p2D : ndarray(float, ndim=2)
        2D array defining east longitude (phi)
        This ranges from 0 to 2*pi and has a shape (nphi,ntheta)
    th2D : ndarray(float, ndim=2)
        2D array defining co-latitude (theta)
        This ranges from 0 to pi and has a shape (nphi,ntheta)
    B : ndarray(float, ndim=2)
        This defines the data to plot corresponding to the grid defined by p2D
        and th2D, usually radial magnetic field
    levels : int
        Number of contour levels
    cmap : str, optional
        Colormap for plotting contours, by default RdBu_r
    proj : str, optional
        Map projection for plotting. Supports all projections used by cartopy,
        by default Mollweide
    vmin : float, optional
        Minimum of colorscale, by default None
    vmax : float, optional
        Maximum of colorscale, by default None

    Returns
    -------
    ax : matplotlib.axes instance
        This is the axes handle of the plot
    cbar : matplotlib.colorbar.Colorbar instance
        This is the handle of the colorbar
    proj : cartopy.crs projection class
        This is the map projection of the plot
    """

    vmin,vmax = get_color_limits(B,vmin,vmax)

    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cs = np.linspace(vmin,vmax,levels)


    lon2D, lat2D, B = lonlat_grid(p2D, th2D, B)

    ccrs = _import_cartopy()
    if ccrs is None:
        proj = 'hammer'

    if proj.lower() == 'hammer':
        ax = plt.axes()
        xx,yy = hammer2cart(lat2D,lon2D)
        cont = ax.contourf(xx,yy,B,cs,cmap=cmap,norm=divnorm,extend='both')
    else:
        ax = plt.axes(projection=getattr(ccrs, proj)())

        cont = ax.contourf(lon2D*180/np.pi,lat2D*180/np.pi,B,cs,  \
            transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm,extend='both')

    cbar = plt.colorbar(cont,orientation='horizontal',fraction=0.06,
                        pad=0.04,ticks=[vmin,0,vmax])

    ax.axis('equal')
    ax.axis('off')

    return ax, cbar, proj

def plotB_subplot(ax,p2D,th2D,B,planetname="earth",levels=60,cmap='RdBu_r',
                  proj='Mollweide',vmin=None,vmax=None):
    """
    Plots subplot of magnetic field on a surface defined by 2D arrays
    of longitude and co-latitude. The subplot is defined by the axes handle ax.

    Parameters
    ----------
    ax :  matplotlib.axes instance
        Axes of the subplot
    p2D : ndarray(float, ndim=2)
        2D array defining east longitude (phi)
        This ranges from 0 to 2*pi and has a shape (nphi,ntheta)
    th2D : ndarray(float, ndim=2)
        2D array defining co-latitude (theta)
        This ranges from 0 to pi and has a shape (nphi,ntheta)
    B : ndarray(float, ndim=2)
        This defines the data to plot corresponding to the grid defined by p2D
        and th2D, usually radial magnetic field
    planetname : str
        Name of the planet
    levels : int
        Number of contour levels
    cmap : str, optional
        Colormap for plotting contours, by default RdBu_r
    proj : str, optional
        Map projection for plotting. Supports all projections used by cartopy,
        by default Mollweide
    vmin : float, optional
        Minimum of colorscale, by default None
    vmax : float, optional
        Maximum of colorscale, by default None

    Returns
    -------
    None
    """
    planetname = planetname.lower()

    lon2D, lat2D, B = lonlat_grid(p2D, th2D, B)

    vmin,vmax = get_color_limits(B,vmin,vmax)
    divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    cs = np.linspace(vmin,vmax,levels)

    ccrs = _import_cartopy()
    if ccrs is None:
        proj = 'hammer'

    if proj.lower() == 'hammer':
        xx,yy = hammer2cart(lat2D,lon2D)
        cont = ax.contourf(xx,yy,B,cs,cmap=cmap,norm=divnorm,extend='both')
    else:
        if planetname == "earth":
            ax.coastlines()

        cont = ax.contourf(lon2D*180/np.pi,lat2D*180/np.pi,B,cs,  \
            transform=ccrs.PlateCarree(),cmap=cmap,norm=divnorm,extend='both')
    cont.set_edgecolor("face")
    cbar = plt.colorbar(cont,orientation='horizontal',fraction=0.06,
                        pad=0.04,norm=divnorm,ticks=[vmin,0,vmax])
    cbar.ax.tick_params(labelsize=15)

    ax.set_title(planetname.capitalize(),fontsize=20)
    ax.axis('equal')
    ax.axis('off')

def plot_spec(l,E,r,planetname):
    """
    Plots Lowes spectrum of a planet.

    Parameters
    ----------
    l : int array
        Array of spherical harmonic degrees from 0 to planet.lmax
    E : array_like
        Array of magnetic energy in each spherical harmonic degree
    r : float
        Radial level of spectrum, r=1 is the planet's surface
    planetname : str
        Name of the planet

    Returns
    -------
    None
    """

    plt.semilogy(l[1:],E[1:],'-o',lw=1.2,color='#449c99',mfc='#347b79',mec='#347b79',ms=8)
    plt.xlabel(r'$l$',fontsize=30)
    plt.ylabel(r'$R_l$ (nT$^2$)',fontsize=30)
    plt.xticks(l)

    plt.tick_params(which='major',labelsize=15,length=10)
    plt.tick_params(which='minor',length=5)

    plt.grid(True,alpha=0.5)

    plt.title('Lowes spectrum, '+ planetname.capitalize() + radius_label(r),fontsize=20,pad=10)
