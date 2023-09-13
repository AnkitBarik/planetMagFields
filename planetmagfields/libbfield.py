#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from .libgauss import get_grid,getB,getBm0
from .plotlib import *
from .utils import planetlist, stdDatDir


def getBr(planet, r=1.0, nphi=256, ntheta=128, info=True):
    """
    Computes radial magnetic field for a planet on a radial surface.
    Radius is scaled to planetary radius.

    Parameters
    ----------
    planet : Planet class instance
        Class containing Gauss coefficients,
    r : float, optional
        Radial level for radial field computation, by default 1.0
    nphi : int, optional
        Number of points in longitude, by default 256
    ntheta : int, optional
        Number of points in co-latitude, by default 128
    info : bool, optional
        Whether to print information about the planet, by default True

    Returns
    -------
    p2D : (2,) array_like
        Longitude at every point on a (longitude,co-latitude) grid
    th2D : (2,) array_like
        Co-latitude at every point on a (longitude,co-latitude) grid
    Br : (2,) array_like
        Radial magnetic field at every point on a (longitude,co-latitude) grid
    dipTheta : float
        Dipole tilt co-latitude in degrees
    dipPhi : float
        Dipole longitude in degrees
    """

    p2D,th2D = get_grid(nphi=nphi,ntheta=ntheta)

    if planet.name in ["mercury", "saturn"]:
        Br = getBm0(planet.lmax,
                    planet.glm,
                    r,
                    p2D,
                    th2D) * 1e-3
        dipTheta = 0.
        dipPhi = 0.
    else:
        Br = getB(planet.lmax,
                  planet.glm,
                  planet.hlm,
                  planet.idx,
                  r,
                  p2D,
                  th2D,
                  planetname=planet.name) * 1e-3

        dipTheta = np.arctan(np.sqrt(planet.glm[planet.idx[1,1]]**2 + planet.hlm[planet.idx[1,1]]**2)
                                    /planet.glm[planet.idx[1,0]]) * 180./np.pi
        dipPhi = np.arctan(planet.hlm[planet.idx[1,1]]/planet.glm[planet.idx[1,1]]) * 180./np.pi

    if info:
        print(("Planet: %s" %planet.name.capitalize()))
        #print(("Depth (fraction of surface radius) = %.2f" %r))
        print(("l_max = %d" %planet.lmax))
        print(("Dipole tilt (degrees) = %f" %dipTheta))

    return p2D, th2D, Br, dipTheta, dipPhi

def plotAllFields(datDir=stdDatDir,r=1.0,levels=30,cmap='RdBu_r',proj='Mollweide'):
    """
    Plots fields of all the planets for which data is available. It's provided in
    utils.planetlist.

    Parameters
    ----------
    datDir : str, optional
        Data directory, where the Gauss coefficient data is present,
        named as <planetname>.dat, the standard directory is ./data,
        by default stdDatDir
    r : float, optional
        Radial level to compute and plot field on, scaled by the planetary
        radius, by default 1.0
    levels : int, optional
        Number of contour levels, by default 30
    cmap : str, optional
        Colormap for contours, by default 'RdBu_r'
    proj : str, optional
        Map projection, by default 'Mollweide'
    """

    from .planet import Planet

    print("")
    print('|=========|======|=======|')
    print(('|%-8s | %-2s| %-5s |' %('Planet','Theta','Phi')))
    print('|=========|======|=======|')

    plt.figure(figsize=(12,12))

    for k, name in enumerate(planetlist):
        planet = Planet(name=name,datDir=datDir,r=r,info=False)

        if name == "ganymede":
            nplot = 8
        else:
            nplot = k+1

        if proj.lower() == 'hammer':
            ax = plt.subplot(3,3,nplot)
        else:
            import cartopy.crs as ccrs
            projection = eval('ccrs.'+proj+'()')
            ax = plt.subplot(3,3,nplot,projection=projection)

        plotB_subplot(ax,planet.p2D,
                      planet.th2D,
                      planet.Br,
                      planetname=name,
                      levels=levels,
                      cmap=cmap,
                      proj=proj)

        if name in ["mercury","saturn"]:
            print(('|%-8s | %-4.1f | %-5.1f |' %(name.capitalize(),planet.dipTheta, planet.dipPhi)))
        else:
            print(('|%-8s | %-3.1f | %-5.1f |' %(name.capitalize(),planet.dipTheta, planet.dipPhi)))

    print('|---------|------|-------|')

    if r == 1:
        plt.suptitle(r'Radial magnetic field ($\mu$T) at surface', fontsize=20)
    else:
        plt.suptitle(r'Radial magnetic field ($\mu$T) at $r/r_{\rm surface} = %.2f$' %r, fontsize=20)